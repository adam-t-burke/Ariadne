using System;
using System.Collections.Generic;
using System.Drawing;
using System.Threading.Tasks;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Display;
using Rhino.Geometry;
using Ariadne.FDM;
using Ariadne.Graphs;
using Theseus.Interop;

namespace Ariadne.Solver;

/// <summary>
/// Grasshopper component that solves FDM networks using the Theseus optimization engine.
/// Uses GH_TaskCapableComponent to run optimization on a background thread.
/// Overrides DrawViewportWires to show intermediate geometry during optimization.
/// </summary>
public class TheseusSolveComponent : GH_TaskCapableComponent<SolveResult>
{
    private SolveResult? _cachedResult;

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
        pManager.AddNumberParameter("Force Densities", "q", "Initial force densities (feed previous Q output for warm-start)", GH_ParamAccess.list, 10.0);
        pManager.AddVectorParameter("Loads", "Loads", "Loads on free nodes", GH_ParamAccess.list, new Vector3d(0, 0, -1));
        pManager.AddGenericParameter("Objectives", "OBJ", "Objective functions to minimize", GH_ParamAccess.list);
        pManager.AddNumberParameter("Lower Bounds", "qMin", "Lower bounds on force densities", GH_ParamAccess.list, 0.1);
        pManager.AddNumberParameter("Upper Bounds", "qMax", "Upper bounds on force densities", GH_ParamAccess.list, 100.0);
        pManager.AddIntegerParameter("Max Iterations", "MaxIter", "Maximum solver iterations", GH_ParamAccess.item, 500);
        pManager.AddBooleanParameter("Optimize", "Opt", "Run optimization (connect a button for on-demand solving)", GH_ParamAccess.item, false);
        pManager.AddIntegerParameter("Report Frequency", "ReportFreq", "Update viewport preview every N evaluations (0 = no preview)", GH_ParamAccess.item, 10);

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
            FDM_Network? network = null;
            List<double> q = [];
            List<Vector3d> loads = [];
            List<Objective> objectives = [];
            List<double> lb = [];
            List<double> ub = [];
            int maxIter = 500;
            bool optimize = false;
            int reportFreq = 10;

            if (!DA.GetData(0, ref network)) return;
            DA.GetDataList(1, q);
            DA.GetDataList(2, loads);
            DA.GetDataList(3, objectives);
            DA.GetDataList(4, lb);
            DA.GetDataList(5, ub);
            DA.GetData(6, ref maxIter);
            DA.GetData(7, ref optimize);
            DA.GetData(8, ref reportFreq);

            if (network == null || !network.Valid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid or null network");
                return;
            }

            if (!optimize && _cachedResult != null)
            {
                OutputResult(DA, _cachedResult);
                return;
            }

            var snapshot = new SolveSnapshot(network, q, loads, objectives, lb, ub, maxIter, optimize, reportFreq);

            ClearPreview();

            var edgeIndices = BuildEdgeIndexPairs(network);
            Action<int, double, double[]>? callback = null;
            if (optimize && reportFreq > 0)
            {
                callback = (iter, loss, xyz) =>
                    OnProgress(iter, loss, xyz, edgeIndices);
            }

            Task<SolveResult> task = Task.Run(() => RunSolve(snapshot, callback));
            TaskList.Add(task);
            return;
        }

        if (!GetSolveResults(DA, out var result))
        {
            if (_cachedResult != null)
                OutputResult(DA, _cachedResult);
            return;
        }

        ClearPreview();

        _cachedResult = result;
        OutputResult(DA, result);

        if (result.Converged)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Remark,
                $"Converged in {result.Iterations} iterations");
        }
        else if (result.Iterations > 1)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                $"Did not converge after {result.Iterations} iterations");
        }
    }

    // ?? Preview drawing ?????????????????????????????????????????

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

    // ?? Callback / preview helpers ??????????????????????????????

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

    // ?? Solve logic ?????????????????????????????????????????????

    private static SolveResult RunSolve(SolveSnapshot snap, Action<int, double, double[]>? progressCallback)
    {
        var solver = new TheseusSolverService();

        var inputs = new SolverInputs
        {
            QInit = snap.Q,
            Loads = snap.Loads,
            LowerBounds = snap.LowerBounds,
            UpperBounds = snap.UpperBounds,
            Objectives = snap.Objectives
        };

        if (snap.Optimize && snap.Objectives.Count > 0)
        {
            var options = new SolverOptions
            {
                MaxIterations = snap.MaxIterations,
                ReportFrequency = snap.ReportFrequency,
            };
            return solver.Solve(snap.Network, inputs, options, progressCallback);
        }

        return solver.SolveForward(snap.Network, inputs);
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
    List<double> Q,
    List<Vector3d> Loads,
    List<Objective> Objectives,
    List<double> LowerBounds,
    List<double> UpperBounds,
    int MaxIterations,
    bool Optimize,
    int ReportFrequency);
