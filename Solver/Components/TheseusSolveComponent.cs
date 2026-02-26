using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Threading;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino.Geometry;
using Ariadne.FDM;
using Ariadne.Graphs;

namespace Ariadne.Solver;

/// <summary>
/// Grasshopper component that solves FDM networks using the Theseus optimization engine.
/// Uses a manual async state machine (not GH_TaskCapableComponent) so the UI thread
/// is never blocked during long-running optimizations. Forward solves run synchronously.
/// Intermediate optimization results are streamed as real component outputs via
/// coalesced ExpireSolution calls, so downstream components receive live data.
/// </summary>
public class TheseusSolveComponent : GH_Component
{
    private enum SolverState { Idle, Running, Done }

    private SolverState _state = SolverState.Idle;
    private int _solveGeneration;
    private SolveResult? _cachedResult;
    private SolveResult? _pendingResult;
    private Exception? _pendingError;
    private bool _lastWasOptimization;

    private bool _prevRunInput;
    private int _triggerGeneration;
    private int _consumedGeneration;
    private int _lastInputHash;

    private CancellationTokenSource? _cts;
    private SynchronizationContext? _uiContext;

    private readonly object _intermediateLock = new();
    private IntermediateOutput? _intermediateOutput;
    private volatile bool _expirePending;

    private record IntermediateOutput(
        List<Point3d> Nodes,
        List<LineCurve> Edges,
        List<double> Q,
        List<double> Forces,
        List<double> Lengths,
        int EvalCount);

    public TheseusSolveComponent()
        : base("Theseus Solve", "Theseus",
            "Solve FDM network using Theseus optimization engine.",
            "Ariadne", "Design")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Network", "Network", "FDM Network to solve", GH_ParamAccess.item);
        pManager.AddNumberParameter("Force Densities", "q", "Initial force densities", GH_ParamAccess.list, 10.0);
        pManager.AddVectorParameter("Loads", "Loads", "Loads on free nodes", GH_ParamAccess.list, new Vector3d(0, 0, -1));
        pManager.AddGenericParameter("Opt Config", "OptConfig", "Optimization configuration (optional, connect for optimization)", GH_ParamAccess.item);

        pManager[3].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Network", "Network", "Solved network with updated geometry", GH_ParamAccess.item);
        pManager.AddPointParameter("Nodes", "Nodes", "Solved node positions", GH_ParamAccess.list);
        pManager.AddCurveParameter("Edges", "Edges", "Solved edge curves", GH_ParamAccess.list);
        pManager.AddNumberParameter("Force Densities", "Q", "Optimized force densities", GH_ParamAccess.list);
        pManager.AddNumberParameter("Member Forces", "Forces", "Member forces (Q x Length)", GH_ParamAccess.list);
        pManager.AddNumberParameter("Member Lengths", "Lengths", "Member lengths", GH_ParamAccess.list);
        pManager.AddIntegerParameter("Iterations", "Iter", "Number of solver iterations", GH_ParamAccess.item);
        pManager.AddBooleanParameter("Converged", "Conv", "Did the solver converge?", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        _uiContext ??= SynchronizationContext.Current;

        // ── Gather inputs (always, so edge detection stays in sync) ──
        FDM_Network? network = null;
        List<double> q = [];
        List<Vector3d> loads = [];
        OptimizationConfig? config = null;

        if (_state != SolverState.Done)
        {
            if (!DA.GetData(0, ref network)) return;
            DA.GetDataList(1, q);
            DA.GetDataList(2, loads);
            DA.GetData(3, ref config);

            if (network == null || !network.Valid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid or null network");
                return;
            }
        }

        // ── Harvest completed optimization result ────────────────
        // Done state skips input gathering, so edge detection is deferred
        // until the next solve where inputs are actually read.
        if (_state == SolverState.Done)
        {
            HarvestResult(DA);
            return;
        }

        // ── Rising-edge detection ────────────────────────────────
        bool currentRun = config?.Run == true;
        if (currentRun && !_prevRunInput)
            _triggerGeneration++;
        _prevRunInput = currentRun;

        // ── No config -> forward solve ───────────────────────────
        if (config == null)
        {
            if (_state == SolverState.Running)
            {
                CancelAndDisposeCts();
                _state = SolverState.Idle;
                Message = "";
            }
            _lastWasOptimization = false;
            var snap = new SolveSnapshot(network!, q.AsReadOnly(), loads.AsReadOnly(), null);
            var result = RunSolve(snap, null);
            _cachedResult = result;
            OutputResult(DA, result);
            return;
        }

        // ── Optimization is in progress ──────────────────────────
        if (_state == SolverState.Running)
        {
            bool newTrigger = _triggerGeneration != _consumedGeneration;
            int inputHash = ComputeInputHash(network!, q, loads, config);
            bool inputsChanged = inputHash != _lastInputHash;

            if (newTrigger || (currentRun && inputsChanged))
            {
                _consumedGeneration = _triggerGeneration;
                _lastInputHash = inputHash;
                StartOptimization(network!, q, loads, config);
            }

            OutputIntermediateOrCached(DA);
            return;
        }

        // ── State is Idle ────────────────────────────────────────

        if (!currentRun)
        {
            if (_cachedResult != null)
                OutputResult(DA, _cachedResult);
            return;
        }

        if (config.Objectives.Count == 0)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                "No objectives provided. Falling back to forward solve.");
            _lastWasOptimization = false;
            var snap = new SolveSnapshot(network!, q.AsReadOnly(), loads.AsReadOnly(), null);
            var result = RunSolve(snap, null);
            _cachedResult = result;
            OutputResult(DA, result);
            return;
        }

        if (_triggerGeneration == _consumedGeneration)
        {
            int inputHash = ComputeInputHash(network!, q, loads, config);
            if (inputHash == _lastInputHash)
            {
                if (_cachedResult != null)
                    OutputResult(DA, _cachedResult);
                return;
            }
            _lastInputHash = inputHash;
        }
        else
        {
            _consumedGeneration = _triggerGeneration;
            _lastInputHash = ComputeInputHash(network!, q, loads, config);
        }

        // ── Start background optimization ────────────────────────
        StartOptimization(network!, q, loads, config);
        OutputIntermediateOrCached(DA);
    }

    // ── Async optimization lifecycle ─────────────────────────────────

    private void StartOptimization(
        FDM_Network network, List<double> q, List<Vector3d> loads, OptimizationConfig config)
    {
        CancelAndDisposeCts();
        _expirePending = false;
        _lastWasOptimization = true;
        _pendingResult = null;
        _pendingError = null;

        lock (_intermediateLock)
            _intermediateOutput = null;

        _cts = new CancellationTokenSource();
        var token = _cts.Token;
        var generation = ++_solveGeneration;

        var edgeIndices = BuildEdgeIndexPairs(network);
        bool stream = config.StreamPreview;

        Func<int, double, double[], double[], bool> callback = (iter, loss, xyz, qVals) =>
        {
            if (token.IsCancellationRequested) return false;
            if (stream)
                OnProgress(iter, loss, xyz, qVals, edgeIndices);
            return true;
        };

        var snap = new SolveSnapshot(network, q.AsReadOnly(), loads.AsReadOnly(), config);

        _state = SolverState.Running;
        Message = "Solving...";

        Task.Run(() =>
        {
            try
            {
                var result = RunSolve(snap, callback);
                if (token.IsCancellationRequested || generation != _solveGeneration)
                    return;

                _pendingResult = result;
                _state = SolverState.Done;
                _uiContext?.Post(_ => ExpireSolution(true), null);
            }
            catch (Exception ex)
            {
                if (token.IsCancellationRequested || generation != _solveGeneration)
                    return;

                _pendingError = ex;
                _state = SolverState.Done;
                _uiContext?.Post(_ => ExpireSolution(true), null);
            }
        });
    }

    private void HarvestResult(IGH_DataAccess DA)
    {
        lock (_intermediateLock)
            _intermediateOutput = null;

        _state = SolverState.Idle;
        Message = "";

        if (_pendingError != null)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, _pendingError.Message);
            _pendingError = null;
            if (_cachedResult != null)
                OutputResult(DA, _cachedResult);
            return;
        }

        if (_pendingResult != null)
        {
            _cachedResult = _pendingResult;
            _pendingResult = null;
            OutputResult(DA, _cachedResult);

            if (_lastWasOptimization)
            {
                if (_cachedResult.Converged)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Remark,
                        $"Converged in {_cachedResult.Iterations} iterations");
                }
                else
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                        $"Did not converge after {_cachedResult.Iterations} iterations");
                }
            }
            return;
        }

        if (_cachedResult != null)
            OutputResult(DA, _cachedResult);
    }

    // ── Progress callback ────────────────────────────────────────────

    private void OnProgress(int evalCount, double loss, double[] xyz, double[] qVals, (int start, int end)[] edgeIndices)
    {
        int nn = xyz.Length / 3;
        var nodes = new List<Point3d>(nn);
        for (int i = 0; i < nn; i++)
            nodes.Add(new Point3d(xyz[i * 3], xyz[i * 3 + 1], xyz[i * 3 + 2]));

        int ne = edgeIndices.Length;
        var edges = new List<LineCurve>(ne);
        var lengths = new List<double>(ne);
        var forces = new List<double>(ne);
        var qList = new List<double>(qVals);

        for (int i = 0; i < ne; i++)
        {
            var (s, e) = edgeIndices[i];
            var lc = new LineCurve(nodes[s], nodes[e]);
            edges.Add(lc);
            double len = lc.GetLength();
            lengths.Add(len);
            forces.Add(i < qVals.Length ? qVals[i] * len : 0.0);
        }

        var intermediate = new IntermediateOutput(nodes, edges, qList, forces, lengths, evalCount);

        lock (_intermediateLock)
            _intermediateOutput = intermediate;

        if (!_expirePending)
        {
            _expirePending = true;
            _uiContext?.Post(_ =>
            {
                _expirePending = false;
                ExpireSolution(true);
            }, null);
        }
    }

    // ── Output helpers ───────────────────────────────────────────────

    private void OutputIntermediateOrCached(IGH_DataAccess DA)
    {
        IntermediateOutput? intermediate;
        lock (_intermediateLock)
            intermediate = _intermediateOutput;

        if (intermediate != null)
        {
            DA.SetDataList(1, intermediate.Nodes);
            DA.SetDataList(2, intermediate.Edges);
            DA.SetDataList(3, intermediate.Q);
            DA.SetDataList(4, intermediate.Forces);
            DA.SetDataList(5, intermediate.Lengths);
            DA.SetData(6, intermediate.EvalCount);
            DA.SetData(7, false);
        }
        else if (_cachedResult != null)
        {
            OutputResult(DA, _cachedResult);
        }
    }

    private void CancelAndDisposeCts()
    {
        _cts?.Cancel();
        _cts?.Dispose();
        _cts = null;
    }

    private static int ComputeInputHash(
        FDM_Network network, List<double> q, List<Vector3d> loads, OptimizationConfig config)
    {
        var hash = new HashCode();
        hash.Add(RuntimeHelpers.GetHashCode(network));
        foreach (var v in q) hash.Add(v);
        foreach (var l in loads) { hash.Add(l.X); hash.Add(l.Y); hash.Add(l.Z); }
        hash.Add(config.MaxIterations);
        hash.Add(config.AbsTol);
        hash.Add(config.RelTol);
        hash.Add(config.BarrierWeight);
        hash.Add(config.BarrierSharpness);
        foreach (var obj in config.Objectives)
            hash.Add(obj.GetContentHashCode());
        foreach (var lb in config.LowerBounds) hash.Add(lb);
        foreach (var ub in config.UpperBounds) hash.Add(ub);
        return hash.ToHashCode();
    }

    private static (int start, int end)[] BuildEdgeIndexPairs(FDM_Network network)
    {
        var graph = network.Graph;
        var nodeMap = new Dictionary<Node, int>(graph.Nn);
        for (int i = 0; i < graph.Nn; i++)
            nodeMap[graph.Nodes[i]] = i;

        var pairs = new (int, int)[graph.Ne];
        for (int i = 0; i < graph.Ne; i++)
        {
            var edge = graph.Edges[i];
            pairs[i] = (nodeMap[edge.Start], nodeMap[edge.End]);
        }
        return pairs;
    }

    // ── Solve logic ─────────────────────────────────────────────────

    private static SolveResult RunSolve(SolveSnapshot snap, Func<int, double, double[], double[], bool>? callback)
    {
        if (snap.Config is { } cfg && cfg.Objectives.Count > 0)
        {
            var inputs = new SolverInputs
            {
                QInit = snap.Q.ToList(),
                Loads = snap.Loads.ToList(),
                LowerBounds = cfg.LowerBounds.ToList(),
                UpperBounds = cfg.UpperBounds.ToList(),
                Objectives = cfg.Objectives.ToList(),
            };
            var options = new SolverOptions
            {
                MaxIterations = cfg.MaxIterations,
                AbsTol = cfg.AbsTol,
                RelTol = cfg.RelTol,
                BarrierWeight = cfg.BarrierWeight,
                BarrierSharpness = cfg.BarrierSharpness,
                ReportFrequency = cfg.ReportFrequency,
            };
            return TheseusSolverService.Solve(snap.Network, inputs, options, callback);
        }

        var fwdInputs = new SolverInputs
        {
            QInit = snap.Q.ToList(),
            Loads = snap.Loads.ToList(),
        };
        return TheseusSolverService.SolveForward(snap.Network, fwdInputs);
    }

    private static void OutputResult(IGH_DataAccess DA, SolveResult result)
    {
        DA.SetData(0, result.Network);
        DA.SetDataList(1, result.NodePositions);
        DA.SetDataList(2, result.EdgeCurves);
        DA.SetDataList(3, result.ForceDensities);
        DA.SetDataList(4, result.MemberForces);
        DA.SetDataList(5, result.MemberLengths);
        DA.SetData(6, result.Iterations);
        DA.SetData(7, result.Converged);
    }

    // ── Cleanup ─────────────────────────────────────────────────────

    public override void RemovedFromDocument(GH_Document document)
    {
        CancelAndDisposeCts();
        lock (_intermediateLock)
            _intermediateOutput = null;
        _state = SolverState.Idle;
        base.RemovedFromDocument(document);
    }

    protected override Bitmap Icon => Properties.Resources.Create;

    public override Guid ComponentGuid => new("F8A7B2C1-3D4E-5F60-A1B2-C3D4E5F60718");
}

/// <summary>
/// Immutable snapshot of all inputs captured on the main thread before
/// the solve runs on a background thread.
/// </summary>
internal sealed record SolveSnapshot(
    FDM_Network Network,
    IReadOnlyList<double> Q,
    IReadOnlyList<Vector3d> Loads,
    OptimizationConfig? Config);
