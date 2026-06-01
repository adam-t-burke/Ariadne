using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino.Geometry;
using Ariadne.FDM;
using Ariadne.Graphs;
using Theseus.Interop;

namespace Ariadne.Solver;

/// <summary>
/// Grasshopper component that solves FDM networks using the Theseus optimization engine.
/// Uses a manual async state machine (not GH_TaskCapableComponent) so the UI thread
/// is never blocked during long-running optimizations or forward solves.
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
    private bool _prevResetCacheInput;
    private int _triggerGeneration;
    private int _consumedGeneration;
    private int _lastInputHash;
    private bool _hasLastInputHash;
    private int _lastResultInputHash;
    private bool _hasLastResultInputHash;

    private CancellationTokenSource? _cts;
    private SynchronizationContext? _uiContext;
    private readonly Dictionary<EdgeKey, double> _qCache = [];
    private PendingCacheUpdate? _pendingCacheUpdate;

    private readonly object _intermediateLock = new();
    private IntermediateOutput? _intermediateOutput;
    private List<double> _intermediateLossTrace = [];
    private volatile bool _expirePending;
    private long _lastPreviewTicks;
    private int _lastPreviewIteration;
    private const int PreviewMinIntervalMs = 200;

    private record IntermediateOutput(
        List<Point3d> Nodes,
        List<LineCurve> Edges,
        List<double> Q,
        List<double> Forces,
        List<double> Lengths,
        int Iterations,
        List<double> LossTrace,
        double FinalLoss);

    private sealed record PendingCacheUpdate(
        bool CacheQ,
        bool IsOptimization,
        EdgeKey[]? EdgeKeys,
        HashSet<EdgeKey>? AmbiguousKeys);

    public TheseusSolveComponent()
        : base("Theseus Solve", "Theseus",
            "Solve FDM network using Theseus optimization engine.",
            "Ariadne", "Design")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Network", "Network", "FDM Network to solve", GH_ParamAccess.item);
        pManager.AddNumberParameter("Force Densities", "q", "Initial force densities", GH_ParamAccess.list, 10.0);
        pManager.AddBooleanParameter("Cache Q", "CacheQ", "Reuse optimized q values for matching edges across solves/topology rebuilds; unmatched edges use input q", GH_ParamAccess.item, false);
        pManager.AddBooleanParameter("Reset Cache", "ResetQ", "Button input that clears stored q values on the rising edge", GH_ParamAccess.item, false);
        pManager.AddVectorParameter("Loads", "Loads", "Loads on free nodes", GH_ParamAccess.list, new Vector3d(0, 0, -1));
        pManager.AddPointParameter("Load Nodes", "LN", "Nodes to apply loads to (optional; if empty, loads apply to all free nodes)", GH_ParamAccess.list);
        pManager.AddGenericParameter("Opt Config", "OptConfig", "Optimization configuration (optional, connect for optimization)", GH_ParamAccess.item);
        pManager.AddGenericParameter("Self-Weight", "SW", "Self-weight configuration (optional)", GH_ParamAccess.item);
        pManager.AddGenericParameter("Pressure", "Press", "Pressure load configuration (optional)", GH_ParamAccess.item);

        pManager[5].Optional = true;
        pManager[6].Optional = true;
        pManager[7].Optional = true;
        pManager[8].Optional = true;
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
        pManager.AddNumberParameter("Loss Trace", "Loss", "Optimization loss values recorded during the solve", GH_ParamAccess.list);
        pManager.AddNumberParameter("Final Loss", "FinalLoss", "Final optimization loss value", GH_ParamAccess.item);
        pManager.AddTextParameter("Termination Reason", "Reason", "Native optimizer termination reason", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        _uiContext ??= SynchronizationContext.Current;

        // ── Gather inputs (always, so edge detection stays in sync) ──
        FDM_Network? network = null;
        List<double> q = [];
        bool cacheQ = false;
        bool resetCache = false;
        List<Vector3d> loads = [];
        List<Point3d> loadNodes = [];
        OptimizationConfig? config = null;
        SelfWeightConfig? selfWeight = null;
        PressureConfig? pressure = null;

        if (_state != SolverState.Done)
        {
            if (!DA.GetData(0, ref network)) return;
            DA.GetDataList(1, q);
            DA.GetData(2, ref cacheQ);
            DA.GetData(3, ref resetCache);
            DA.GetDataList(4, loads);
            DA.GetDataList(5, loadNodes);
            DA.GetData(6, ref config);
            DA.GetData(7, ref selfWeight);
            DA.GetData(8, ref pressure);

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

        bool resetTriggered = resetCache && !_prevResetCacheInput;
        if (resetTriggered)
        {
            _qCache.Clear();
            _cachedResult = null;
            _hasLastInputHash = false;
            _hasLastResultInputHash = false;
            if (currentRun)
                _triggerGeneration++;
            Message = cacheQ ? "q cache reset" : "";
        }
        _prevResetCacheInput = resetCache;

        int inputHash = ComputeInputHash(network!, q, cacheQ, loads, loadNodes, config, selfWeight, pressure);

        // ── No config -> forward solve (async) ───────────────────
        if (config == null)
        {
            if (_state == SolverState.Running)
            {
                if (_lastWasOptimization)
                {
                    CancelAndDisposeCts();
                    _state = SolverState.Idle;
                    Message = "";
                }
                else
                {
                    OutputIntermediateOrCached(DA);
                    return;
                }
            }

            if (IsKnownInputHash(inputHash))
            {
                OutputIntermediateOrCached(DA);
                return;
            }

            _lastWasOptimization = false;
            var snap = CreateSnapshot(network!, q, cacheQ, loads, loadNodes, null, selfWeight, pressure);
            StartAsyncSolve(snap, streamPreview: false, inputHash);
            OutputIntermediateOrCached(DA);
            return;
        }

        // ── Optimization is in progress ──────────────────────────
        if (_state == SolverState.Running)
        {
            bool newTrigger = _triggerGeneration != _consumedGeneration;
            bool inputsChanged = !_hasLastInputHash || inputHash != _lastInputHash;

            if (newTrigger || (currentRun && inputsChanged))
            {
                _consumedGeneration = _triggerGeneration;
                StartOptimization(network!, q, cacheQ, loads, loadNodes, config, selfWeight, pressure, inputHash);
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
            if (IsKnownInputHash(inputHash))
            {
                OutputIntermediateOrCached(DA);
                return;
            }

            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                "No objectives provided. Falling back to forward solve.");
            _lastWasOptimization = false;
            var snap = CreateSnapshot(network!, q, cacheQ, loads, loadNodes, config, selfWeight, pressure);
            StartAsyncSolve(snap, streamPreview: false, inputHash);
            OutputIntermediateOrCached(DA);
            return;
        }

        if (_triggerGeneration == _consumedGeneration)
        {
            if (IsKnownInputHash(inputHash))
            {
                if (_cachedResult != null)
                    OutputResult(DA, _cachedResult);
                return;
            }
        }
        else
        {
            _consumedGeneration = _triggerGeneration;
        }

        // ── Start background optimization ────────────────────────
        StartAsyncSolve(
            CreateSnapshot(network!, q, cacheQ, loads, loadNodes, config, selfWeight, pressure),
            streamPreview: config.StreamPreview,
            inputHash);
        OutputIntermediateOrCached(DA);
    }

    // ── Async solve lifecycle ─────────────────────────────────

    private void StartAsyncSolve(SolveSnapshot snap, bool streamPreview, int inputHash)
    {
        CancelAndDisposeCts();
        _expirePending = false;
        _lastPreviewTicks = 0;
        _lastPreviewIteration = 0;
        _pendingResult = null;
        _pendingError = null;
        _pendingCacheUpdate = null;
        _lastInputHash = inputHash;
        _hasLastInputHash = true;
        _hasLastResultInputHash = false;

        lock (_intermediateLock)
        {
            _intermediateOutput = null;
            _intermediateLossTrace = [];
        }

        _cts = new CancellationTokenSource();
        var token = _cts.Token;
        var generation = ++_solveGeneration;

        bool isOptimization = snap.Config is { Objectives.Count: > 0 };
        _lastWasOptimization = isOptimization;

        Func<int, double, double[], double[], bool>? callback = null;
        (int start, int end)[]? edgeIndices = null;
        if (isOptimization)
        {
            edgeIndices = BuildEdgeIndexPairs(snap.Network);
            bool stream = streamPreview;
            callback = (iter, loss, xyz, qVals) =>
            {
                if (token.IsCancellationRequested) return false;
                if (stream)
                    OnProgress(iter, loss, xyz, qVals, edgeIndices);
                return true;
            };
        }

        _state = SolverState.Running;
        Message = "Solving...";

        Task.Run(() =>
        {
            try
            {
                var result = RunSolve(snap, callback);
                int resultInputHash = ComputeInputHash(
                    snap.Network,
                    snap.Q,
                    snap.CacheQ,
                    snap.Loads,
                    snap.LoadNodes,
                    snap.Config,
                    snap.SelfWeight,
                    snap.Pressure);

                if (token.IsCancellationRequested || generation != _solveGeneration)
                    return;

                _pendingResult = result;
                _pendingCacheUpdate = new PendingCacheUpdate(
                    snap.CacheQ,
                    isOptimization,
                    snap.EdgeKeys,
                    snap.AmbiguousKeys);
                _lastResultInputHash = resultInputHash;
                _hasLastResultInputHash = true;
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

    private void StartOptimization(
        FDM_Network network, List<double> q, bool cacheQ, List<Vector3d> loads, List<Point3d> loadNodes,
        OptimizationConfig config, SelfWeightConfig? selfWeight, PressureConfig? pressure, int inputHash)
    {
        StartAsyncSolve(
            CreateSnapshot(network, q, cacheQ, loads, loadNodes, config, selfWeight, pressure),
            streamPreview: config.StreamPreview,
            inputHash);
    }

    private void HarvestResult(IGH_DataAccess DA)
    {
        lock (_intermediateLock)
        {
            _intermediateOutput = null;
            _intermediateLossTrace = [];
        }

        _state = SolverState.Idle;
        Message = "";

        if (_pendingError != null)
        {
            var ex = _pendingError;
            _pendingError = null;
            string message = ex is TheseusException tex && tex.NativeCode != 0
                ? $"{ex.Message} (native code {tex.NativeCode})"
                : ex.Message;
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, message);
            if (_cachedResult != null)
                OutputResult(DA, _cachedResult);
            return;
        }

        if (_pendingResult != null)
        {
            _cachedResult = _pendingResult;
            _pendingResult = null;
            UpdateQCacheIfNeeded(_cachedResult, _pendingCacheUpdate);
            _pendingCacheUpdate = null;
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
                    var reason = string.IsNullOrWhiteSpace(_cachedResult.TerminationReason)
                        ? ""
                        : $": {_cachedResult.TerminationReason}";
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                        $"Did not converge after {_cachedResult.Iterations} iterations{reason}");
                }
            }
            return;
        }

        if (_cachedResult != null)
            OutputResult(DA, _cachedResult);
    }

    // ── Progress callback ────────────────────────────────────────────

    private void OnProgress(int iteration, double loss, double[] xyz, double[] qVals, (int start, int end)[] edgeIndices)
    {
        long now = Stopwatch.GetTimestamp();
        double elapsedMs = (now - _lastPreviewTicks) * 1000.0 / Stopwatch.Frequency;
        if (_lastPreviewTicks != 0
            && elapsedMs < PreviewMinIntervalMs
            && iteration - _lastPreviewIteration < 1)
            return;

        _lastPreviewTicks = now;
        _lastPreviewIteration = iteration;
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

        lock (_intermediateLock)
        {
            _intermediateLossTrace.Add(loss);
            var lossTrace = _intermediateLossTrace.ToList();
            var intermediate = new IntermediateOutput(
                nodes,
                edges,
                qList,
                forces,
                lengths,
                iteration,
                lossTrace,
                loss);
            _intermediateOutput = intermediate;
        }

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
            DA.SetData(6, intermediate.Iterations);
            DA.SetData(7, false);
            DA.SetDataList(8, intermediate.LossTrace);
            DA.SetData(9, intermediate.FinalLoss);
            DA.SetData(10, "running");
        }
        else if (_cachedResult != null)
        {
            OutputResult(DA, _cachedResult);
        }
    }

    private void CancelAndDisposeCts()
    {
        TheseusSolverService.RequestActiveCancel();
        _cts?.Cancel();
        _cts?.Dispose();
        _cts = null;
    }

    private static int ComputeInputHash(
        FDM_Network network,
        IReadOnlyList<double> q,
        bool cacheQ,
        IReadOnlyList<Vector3d> loads,
        IReadOnlyList<Point3d> loadNodes,
        OptimizationConfig? config,
        SelfWeightConfig? selfWeight,
        PressureConfig? pressure)
    {
        var hash = new HashCode();
        hash.Add(network.GetTopologyHashCode());
        hash.Add(cacheQ);
        foreach (var v in q) hash.Add(v);
        foreach (var l in loads) { hash.Add(l.X); hash.Add(l.Y); hash.Add(l.Z); }
        foreach (var ln in loadNodes) { hash.Add(ln.X); hash.Add(ln.Y); hash.Add(ln.Z); }
        hash.Add(config is not null);
        if (config is not null)
        {
            hash.Add(config.MaxIterations);
            hash.Add(config.AbsTol);
            hash.Add(config.RelTol);
            hash.Add(config.BarrierWeight);
            hash.Add(config.BarrierSharpness);
            hash.Add(config.ReportFrequency);
            hash.Add(config.AnchorSaturationLambda);
            hash.Add(config.QParameterizationMode);
            hash.Add(config.StreamPreview);
            foreach (var obj in config.Objectives)
                hash.Add(obj.GetContentHashCode());
            foreach (var support in config.VariableSupports)
                hash.Add(support.GetContentHashCode());
            foreach (var lb in config.LowerBounds) hash.Add(lb);
            foreach (var ub in config.UpperBounds) hash.Add(ub);
        }
        hash.Add(selfWeight?.GetContentHashCode() ?? 0);
        hash.Add(pressure?.GetContentHashCode() ?? 0);
        return hash.ToHashCode();
    }

    private SolveSnapshot CreateSnapshot(
        FDM_Network network,
        IReadOnlyList<double> inputQ,
        bool cacheQ,
        IReadOnlyList<Vector3d> loads,
        IReadOnlyList<Point3d> loadNodes,
        OptimizationConfig? config,
        SelfWeightConfig? selfWeight,
        PressureConfig? pressure)
    {
        var resolvedQ = ExpandInputQ(inputQ, network.Graph.Ne);
        EdgeKey[]? edgeKeys = null;
        HashSet<EdgeKey>? ambiguousKeys = null;

        if (cacheQ)
        {
            var keySet = EdgeIdentity.BuildKeys(network);
            edgeKeys = keySet.Keys;
            ambiguousKeys = keySet.AmbiguousKeys;
            for (int i = 0; i < edgeKeys.Length; i++)
            {
                var key = edgeKeys[i];
                if (!ambiguousKeys.Contains(key) && _qCache.TryGetValue(key, out double cachedQ))
                    resolvedQ[i] = cachedQ;
            }

            if (ambiguousKeys.Count > 0)
            {
                AddRuntimeMessage(
                    GH_RuntimeMessageLevel.Warning,
                    $"{ambiguousKeys.Count} ambiguous edge identity key(s); input q was used for those edges.");
            }
        }

        return new SolveSnapshot(
            network,
            resolvedQ.AsReadOnly(),
            loads.ToList().AsReadOnly(),
            loadNodes.ToList().AsReadOnly(),
            config,
            cacheQ,
            edgeKeys,
            ambiguousKeys,
            selfWeight,
            pressure);
    }

    private static List<double> ExpandInputQ(IReadOnlyList<double> inputQ, int edgeCount)
    {
        if (inputQ.Count == 0)
            throw new ArgumentException("Initial force densities cannot be empty.", nameof(inputQ));

        var result = new List<double>(edgeCount);
        for (int i = 0; i < edgeCount; i++)
            result.Add(i < inputQ.Count ? inputQ[i] : inputQ[^1]);
        return result;
    }

    private void UpdateQCacheIfNeeded(SolveResult result, PendingCacheUpdate? update)
    {
        if (update is not { CacheQ: true, EdgeKeys: { } keys })
            return;
        if (update.IsOptimization && !result.Converged)
            return;

        var ambiguous = update.AmbiguousKeys ?? [];
        var seen = new HashSet<EdgeKey>();
        int n = Math.Min(keys.Length, result.ForceDensities.Length);
        for (int i = 0; i < n; i++)
        {
            var key = keys[i];
            if (ambiguous.Contains(key))
                continue;
            _qCache[key] = result.ForceDensities[i];
            seen.Add(key);
        }

        foreach (var key in _qCache.Keys.ToList())
        {
            if (!seen.Contains(key))
                _qCache.Remove(key);
        }

        Message = $"q cache: {_qCache.Count}";
    }

    private bool IsKnownInputHash(int inputHash) =>
        (_hasLastInputHash && inputHash == _lastInputHash)
        || (_hasLastResultInputHash && inputHash == _lastResultInputHash);

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
        var loadNodeIndices = snap.LoadNodes.Count > 0
            ? TheseusSolverService.ResolveLoadNodeIndices(snap.Network, snap.LoadNodes)
            : null;

        if (snap.Config is { } cfg && cfg.Objectives.Count > 0)
        {
            var inputs = new SolverInputs
            {
                QInit = snap.Q.ToList(),
                Loads = snap.Loads.ToList(),
                LoadNodeIndices = loadNodeIndices,
                LowerBounds = cfg.LowerBounds.ToList(),
                UpperBounds = cfg.UpperBounds.ToList(),
                Objectives = cfg.Objectives.ToList(),
                QParameterizationMode = cfg.QParameterizationMode,
                SelfWeight = snap.SelfWeight,
                Pressure = snap.Pressure,
                VariableSupports = cfg.VariableSupports.ToList(),
            };
            var options = cfg.ToSolverOptions();
            return TheseusSolverService.Solve(snap.Network, inputs, options, callback);
        }

        var fwdInputs = new SolverInputs
        {
            QInit = snap.Q.ToList(),
            Loads = snap.Loads.ToList(),
            LoadNodeIndices = loadNodeIndices,
            SelfWeight = snap.SelfWeight,
            Pressure = snap.Pressure,
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
        DA.SetDataList(8, result.LossTrace);
        if (result.LossTrace.Length > 0)
            DA.SetData(9, result.FinalLoss);
        DA.SetData(10, result.TerminationReason);
    }

    // ── Cleanup ─────────────────────────────────────────────────────

    public override void RemovedFromDocument(GH_Document document)
    {
        CancelAndDisposeCts();
        lock (_intermediateLock)
        {
            _intermediateOutput = null;
            _intermediateLossTrace = [];
        }
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
    IReadOnlyList<Point3d> LoadNodes,
    OptimizationConfig? Config,
    bool CacheQ,
    EdgeKey[]? EdgeKeys = null,
    HashSet<EdgeKey>? AmbiguousKeys = null,
    SelfWeightConfig? SelfWeight = null,
    PressureConfig? Pressure = null);
