using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Display;
using Rhino.Geometry;
using Ariadne.FDM;
using Ariadne.Graphs;

namespace Ariadne.Solver;

/// <summary>
/// Grasshopper component that solves FDM networks using the Theseus optimization engine.
/// Uses GH_TaskCapableComponent to run optimization on a background thread.
/// Overrides DrawViewportWires to show intermediate geometry during optimization.
/// </summary>
public class TheseusSolveComponent : GH_TaskCapableComponent<SolveResult>
{
    private SolveResult? _cachedResult;
    private bool _lastWasOptimization;

    private CancellationTokenSource? _cts;

    private readonly object _previewLock = new();
    private Point3d[]? _previewPoints;
    private Line[]? _previewLines;
    private string? _previewStatus;

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
        if (InPreSolve)
        {
            // ── Gather inputs ──────────────────────────────────────
            FDM_Network? network = null;
            List<double> q = [];
            List<Vector3d> loads = [];
            OptimizationConfig? config = null;

            if (!DA.GetData(0, ref network)) return;
            DA.GetDataList(1, q);
            DA.GetDataList(2, loads);
            DA.GetData(3, ref config);

            if (network == null || !network.Valid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid or null network");
                return;
            }

            // ── State machine ──────────────────────────────────────
            if (config == null)
            {
                // No OptConfig: cancel any in-progress optimization, forward solve
                CancelAndDisposeCts();
                ClearPreview();
                _lastWasOptimization = false;

                var snap = new SolveSnapshot(network, q.AsReadOnly(), loads.AsReadOnly(), null);
                TaskList.Add(Task.Run(() => RunSolve(snap, null)));
                return;
            }

            if (!config.Run)
            {
                // Run = false: do nothing. Let any in-progress optimization finish.
                // Pass 2 will pick up a completed result or output cache.
                return;
            }

            // Config.Run = true: cancel old optimization, start new one
            if (config.Objectives.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                    "No objectives provided. Falling back to forward solve.");
                CancelAndDisposeCts();
                ClearPreview();
                _lastWasOptimization = false;

                var snap = new SolveSnapshot(network, q.AsReadOnly(), loads.AsReadOnly(), null);
                TaskList.Add(Task.Run(() => RunSolve(snap, null)));
                return;
            }

            // Optimization with objectives
            CancelAndDisposeCts();
            ClearPreview();
            _lastWasOptimization = true;

            _cts = new CancellationTokenSource();
            var token = _cts.Token;

            var edgeIndices = BuildEdgeIndexPairs(network);
            bool streamPreview = config.StreamPreview;

            Func<int, double, double[], bool> callback = (iter, loss, xyz) =>
            {
                if (token.IsCancellationRequested) return false;
                if (streamPreview)
                    OnProgress(iter, loss, xyz, edgeIndices);
                return true;
            };

            var optSnap = new SolveSnapshot(network, q.AsReadOnly(), loads.AsReadOnly(), config);
            TaskList.Add(Task.Run(() => RunSolve(optSnap, callback)));
            return;
        }

        // ── Pass 2 ─────────────────────────────────────────────────
        SolveResult? result;
        bool hasResult;
        try
        {
            hasResult = GetSolveResults(DA, out result);
        }
        catch
        {
            hasResult = false;
            result = null;
        }

        if (!hasResult)
        {
            if (_cachedResult != null)
                OutputResult(DA, _cachedResult);
            return;
        }

        ClearPreview();

        _cachedResult = result;
        OutputResult(DA, result!);

        if (_lastWasOptimization)
        {
            if (result!.Converged)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Remark,
                    $"Converged in {result.Iterations} iterations");
            }
            else
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                    $"Did not converge after {result.Iterations} iterations");
            }
        }
    }

    // ── Preview drawing ─────────────────────────────────────────────

    public override void DrawViewportWires(IGH_PreviewArgs args)
    {
        base.DrawViewportWires(args);

        Point3d[]? pts;
        Line[]? lines;
        string? status;

        lock (_previewLock)
        {
            pts = _previewPoints;
            lines = _previewLines;
            status = _previewStatus;
        }

        if (pts == null && lines == null) return;

        var lineColor = Color.FromArgb(200, 50, 200, 50);
        var pointColor = Color.FromArgb(200, 255, 165, 0);

        if (lines != null)
        {
            foreach (var line in lines)
                args.Display.DrawLine(line, lineColor, 2);
        }

        if (pts != null)
        {
            args.Display.DrawPoints(pts, PointStyle.RoundSimple, 4, pointColor);
        }

        if (status != null)
        {
            var bbox = BoundingBox.Empty;
            if (pts != null)
            {
                foreach (var p in pts) bbox.Union(p);
            }
            if (bbox.IsValid)
            {
                var anchor = bbox.Center + new Vector3d(0, 0, bbox.Diagonal.Z * 0.6);
                args.Display.Draw2dText(status, Color.White, anchor, false, 14);
            }
        }
    }

    public override BoundingBox ClippingBox
    {
        get
        {
            var bbox = base.ClippingBox;
            lock (_previewLock)
            {
                if (_previewPoints != null)
                {
                    foreach (var p in _previewPoints) bbox.Union(p);
                }
            }
            return bbox;
        }
    }

    // ── Callback / preview helpers ──────────────────────────────────

    private void OnProgress(int evalCount, double loss, double[] xyz, (int start, int end)[] edgeIndices)
    {
        int nn = xyz.Length / 3;
        var pts = new Point3d[nn];
        for (int i = 0; i < nn; i++)
            pts[i] = new Point3d(xyz[i * 3], xyz[i * 3 + 1], xyz[i * 3 + 2]);

        var lines = new Line[edgeIndices.Length];
        for (int i = 0; i < edgeIndices.Length; i++)
        {
            var (s, e) = edgeIndices[i];
            lines[i] = new Line(pts[s], pts[e]);
        }

        lock (_previewLock)
        {
            _previewPoints = pts;
            _previewLines = lines;
            _previewStatus = $"eval {evalCount}  loss {loss:E3}";
        }

        RhinoApp.InvokeOnUiThread(new Action<object?>(_ =>
        {
            RhinoDoc.ActiveDoc?.Views.Redraw();
        }), null);
    }

    private void ClearPreview()
    {
        lock (_previewLock)
        {
            _previewPoints = null;
            _previewLines = null;
            _previewStatus = null;
        }
    }

    private void CancelAndDisposeCts()
    {
        _cts?.Cancel();
        _cts?.Dispose();
        _cts = null;
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

    private static SolveResult RunSolve(SolveSnapshot snap, Func<int, double, double[], bool>? callback)
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
