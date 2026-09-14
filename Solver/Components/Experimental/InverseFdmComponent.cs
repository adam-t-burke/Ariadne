using System;
using System.Collections.Generic;
using System.Drawing;
using System.Globalization;
using System.Linq;
using System.Windows.Forms;
using GH_IO.Serialization;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using Ariadne.FDM;
using Ariadne.Solver;
using Theseus.Interop;

namespace Ariadne.Solver.Components.Experimental;

/// <summary>
/// Experimental inverse force-density solve at a target geometry.
/// </summary>
public class InverseFdmComponent : GH_Component
{
    private const string ParticularKey = "InverseFdmParticular";
    private const string LinearAlgebraKey = "InverseFdmLinearAlgebra";
    private const string MetricKey = "InverseFdmMetric";
    private const string Stage2Key = "InverseFdmStage2";
    private const string NondimensionalizeKey = "InverseFdmNondimensionalize";
    private ParticularMode _particular = InverseFdmUiState.DefaultParticular;
    private LinearAlgebraMode _linearAlgebra = LinearAlgebraMode.Direct;
    private MetricMode _metric = InverseFdmUiState.DefaultMetric;
    private Stage2Mode _stage2 = InverseFdmUiState.DefaultStage2;
    private bool _nondimensionalize = InverseFdmUiState.DefaultNondimensionalize;
    private double _lambda = 1e-6;
    private double _cwlsDamping = 1e-6;
    private double _seedGuardMargin = InverseFdmUiState.DefaultSeedGuardMargin;
    private int _frozenIterations = InverseFdmUiState.DefaultFrozenIterations;
    private int _gnIterations = InverseFdmUiState.DefaultGnIterations;
    private bool _solveForQ = InverseFdmUiState.DefaultSolveForQ;
    private bool _hasBox;

    public InverseFdmComponent()
        : base("Inverse FDM", "InvFDM",
            "Compliance-weighted warm start for a target geometry: Stage-1 particular, seed guard, frozen CWLS step(s), Gauss–Newton step(s), then forward-solve.",
            "Ariadne", "Experimental")
    {
        UpdateMessage();
    }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Network", "Network", "FDM Network (topology + anchors)", GH_ParamAccess.item);
        pManager.AddPointParameter("Target Points", "Target", "Desired free-node positions (one per free node, matching order)", GH_ParamAccess.list);
        pManager.AddVectorParameter("Loads", "Loads", "Loads on free nodes", GH_ParamAccess.list, new Vector3d(0, 0, -1));
        pManager.AddPointParameter("Load Nodes", "LN", "Nodes to apply loads to (optional; if empty, loads apply to all free nodes)", GH_ParamAccess.list);
        pManager.AddNumberParameter("Regularization", "λ", "Stage-1 particular regularization used by Tikhonov, Gram (sparse/dense), LSQR, Clarabel, and SPG. Ignored for Moore–Penrose and QR.", GH_ParamAccess.item, 1e-6);
        pManager.AddIntegerParameter("Frozen CWLS Iterations", "FrozenIter", "Geometric metric only: maximum frozen-target CWLS updates before Gauss–Newton. 0 skips this phase. Default 1 (the benchmark pipeline).", GH_ParamAccess.item, InverseFdmUiState.DefaultFrozenIterations);
        pManager.AddIntegerParameter("Gauss–Newton Iterations", "GNiter", "Geometric metric only: maximum CWLS-GN updates after the frozen phase. 0 skips this phase. Stops early at Tol. Default 2 (the benchmark pipeline).", GH_ParamAccess.item, InverseFdmUiState.DefaultGnIterations);
        pManager.AddBooleanParameter("Enforce Rx=0", "Rx0", "Strictly enforce zero X-reaction at supports", GH_ParamAccess.item, false);
        pManager.AddBooleanParameter("Enforce Ry=0", "Ry0", "Strictly enforce zero Y-reaction at supports", GH_ParamAccess.item, false);
        pManager.AddBooleanParameter("Enforce Rz=0", "Rz0", "Strictly enforce zero Z-reaction at supports", GH_ParamAccess.item, false);
        pManager.AddBooleanParameter("Solve Q", "SolveQ", "Stage 1 only: True solves the initial particular in force densities q. False (default) solves member forces t, then recovers q = t / target length. CWLS and L-BFGS-B operate in q.", GH_ParamAccess.item, InverseFdmUiState.DefaultSolveForQ);
        pManager.AddIntegerParameter("Signs", "Signs",
            "+1 tension (q ≥ 0), -1 compression (q ≤ 0), 0 free. Bounds always apply to q, even when Stage 1 solves member forces.",
            GH_ParamAccess.tree);
        pManager.AddNumberParameter("Lower", "Lower",
            "Lower bound on q. For a force-space Stage 1 this is internally multiplied by target edge length.",
            GH_ParamAccess.tree);
        pManager.AddNumberParameter("Upper", "Upper",
            "Upper bound on q. For a force-space Stage 1 this is internally multiplied by target edge length.",
            GH_ParamAccess.tree);
        pManager.AddIntegerParameter("Max Iterations", "MaxIter", "Iteration budget per inner solve for Clarabel, SPG, and LSQR", GH_ParamAccess.item, 500);
        pManager.AddNumberParameter("Tolerance", "Tol", "Convergence tolerance for Clarabel, SPG, and LSQR", GH_ParamAccess.item, 1e-6);
        pManager.AddNumberParameter("CWLS Damping", "λcwls", "Stage-2 Tikhonov floor λcwls‖Δq‖² on every compliance-weighted step. Separate from Stage-1 particular regularization.", GH_ParamAccess.item, 1e-6);
        pManager.AddNumberParameter("Seed Guard", "Guard",
            "Stage-1 collapse guard margin. After Stage 1 a scaled uniform sign seed is scored on the same geometric error; when Stage 1 is worse by more than this factor both seeds run Stage 2 and the better result continues. 0 disables the guard (the benchmark's pipeline_noguard / legacy rows).",
            GH_ParamAccess.item, InverseFdmUiState.DefaultSeedGuardMargin);
        pManager.AddNumberParameter("LM Damping", "λLM",
            "Levenberg–Marquardt floor for the Gauss–Newton steps, relative to the curvature diagonal. 0 (default) takes the undamped direction and halves the step on the exact merit; positive values damp the direction and grow ×10 on rejected steps.",
            GH_ParamAccess.item, InverseFdmUiState.DefaultLmDamping);
        pManager.AddNumberParameter("Reaction Weight", "wR",
            "Weight of the Rx0/Ry0/Rz0 rows relative to the equilibrium (Stage 1) and geometric (Stage 2) rows. 1 (default) counts one load unit of reaction like one target extent of geometric error.",
            GH_ParamAccess.item, InverseFdmUiState.DefaultReactionWeight);
        pManager[3].Optional = true;
        pManager[11].Optional = true;
        pManager[12].Optional = true;
        pManager[13].Optional = true;
        pManager[17].Optional = true;
        pManager[18].Optional = true;
        pManager[19].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Network", "Network", "Solved network with updated geometry", GH_ParamAccess.item);
        pManager.AddPointParameter("Nodes", "Nodes", "Solved node positions", GH_ParamAccess.list);
        pManager.AddCurveParameter("Edges", "Edges", "Solved edge curves", GH_ParamAccess.list);
        pManager.AddNumberParameter("Force Densities", "Q", "Computed force densities", GH_ParamAccess.list);
        pManager.AddNumberParameter("Member Forces", "Forces", "Target-geometry member forces", GH_ParamAccess.list);
        pManager.AddVectorParameter("Residual", "Residual", "Free-node equilibrium residual at the target", GH_ParamAccess.list);
        pManager.AddNumberParameter("Residual Ratio", "RelRes", "Residual norm divided by load norm", GH_ParamAccess.item);
        pManager.AddNumberParameter("Geometric Error", "GeomErr",
            "‖x(q) − x*‖: distance from the forward-solved geometry to the target. Unlike RelRes this is a length, not a force ratio. Divide by the target's bounding-box diagonal to compare with the benchmark's err/L.",
            GH_ParamAccess.item);
        pManager.AddBooleanParameter("Guard Used", "Guard",
            "True when the seed guard replaced the Stage-1 particular by the scaled uniform sign seed.",
            GH_ParamAccess.item);
        pManager.AddTextParameter("Diagnostics", "Diag",
            "Stage-2 diagnostics as key = value lines: Stage-1 and uniform-seed errors, accepted frozen / Gauss–Newton steps, factorisations, Clarabel fallbacks, active-set pass-limit hits, degenerate linearisations, realised reaction residual.",
            GH_ParamAccess.list);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        DA.DisableGapLogic();

        FDM_Network? network = null;
        List<Point3d> targetPoints = [];
        List<Vector3d> loads = [];
        List<Point3d> loadNodes = [];
        double regularization = 1e-6;
        const int maxL1Iter = 20;
        int frozenIterations = InverseFdmUiState.DefaultFrozenIterations;
        int gnIterations = InverseFdmUiState.DefaultGnIterations;
        bool enforceZeroRx = false;
        bool enforceZeroRy = false;
        bool enforceZeroRz = false;
        bool solveForQ = false;
        var signTree = new GH_Structure<GH_Integer>();
        var lowerTree = new GH_Structure<GH_Number>();
        var upperTree = new GH_Structure<GH_Number>();
        int maxIter = 500;
        double tol = 1e-6;
        double cwlsDamping = 1e-6;
        double seedGuardMargin = InverseFdmUiState.DefaultSeedGuardMargin;
        double lmDamping = InverseFdmUiState.DefaultLmDamping;
        double reactionWeight = InverseFdmUiState.DefaultReactionWeight;

        if (!DA.GetData(0, ref network)) return;
        if (!DA.GetDataList(1, targetPoints)) return;
        DA.GetDataList(2, loads);
        DA.GetDataList(3, loadNodes);
        DA.GetData(4, ref regularization);
        DA.GetData(5, ref frozenIterations);
        DA.GetData(6, ref gnIterations);
        DA.GetData(7, ref enforceZeroRx);
        DA.GetData(8, ref enforceZeroRy);
        DA.GetData(9, ref enforceZeroRz);
        DA.GetData(10, ref solveForQ);
        DA.GetDataTree(11, out signTree);
        DA.GetDataTree(12, out lowerTree);
        DA.GetDataTree(13, out upperTree);
        DA.GetData(14, ref maxIter);
        DA.GetData(15, ref tol);
        DA.GetData(16, ref cwlsDamping);
        DA.GetData(17, ref seedGuardMargin);
        DA.GetData(18, ref lmDamping);
        DA.GetData(19, ref reactionWeight);

        _lambda = regularization;
        _cwlsDamping = cwlsDamping;
        _seedGuardMargin = seedGuardMargin;
        _frozenIterations = Math.Max(0, frozenIterations);
        _gnIterations = Math.Max(0, gnIterations);
        _solveForQ = solveForQ;

        if (network == null || !network.Valid)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid or null network.");
            return;
        }

        if (!TryMapOptionalTree(signTree, network.Graph.EdgeInputPaths, "Signs",
                out int[] signs, out string? signError, out string? signWarning))
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, signError!);
            return;
        }
        if (!TryMapOptionalTree(lowerTree, network.Graph.EdgeInputPaths, "Lower",
                out double[] lower, out string? lowerError, out string? lowerWarning))
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, lowerError!);
            return;
        }
        if (!TryMapOptionalTree(upperTree, network.Graph.EdgeInputPaths, "Upper",
                out double[] upper, out string? upperError, out string? upperWarning))
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, upperError!);
            return;
        }

        if (signWarning is not null)
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, signWarning);
        if (lowerWarning is not null)
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, lowerWarning);
        if (upperWarning is not null)
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, upperWarning);

        _hasBox = InverseFdmUiState.HasEffectiveBounds(signs, lower, upper);
        _particular = InverseFdmUiState.UpdateParticular(_linearAlgebra, _particular, _hasBox);
        UpdateMessage();

        int numFree = network.Free.Count;
        if (targetPoints.Count != numFree)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                $"Target points count ({targetPoints.Count}) must match free node count ({numFree}).");
            return;
        }

        bool unconstrainedDirect = _linearAlgebra == LinearAlgebraMode.Direct && !_hasBox;
        if (unconstrainedDirect && _particular == ParticularMode.Tikhonov && regularization <= 0.0)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                "Tikhonov mode requires λ > 0. Switch Particular to Moore–Penrose for λ = 0.");
            return;
        }

        if (unconstrainedDirect
            && _particular == ParticularMode.GramDense
            && network.Graph.Ne > InverseFdmUiState.DenseGramWarnEdges)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                $"Gram (dense) forms EᵀE as a full {network.Graph.Ne}×{network.Graph.Ne} matrix: "
                + "O(ne²) memory and O(ne³) time. It is the reference for what sparsity buys; "
                + "use Gram (sparse) or Tikhonov for the same minimiser at sparse cost.");
        }

        double[] targetFreeXyz = new double[numFree * 3];
        for (int i = 0; i < numFree; i++)
        {
            targetFreeXyz[i * 3 + 0] = targetPoints[i].X;
            targetFreeXyz[i * 3 + 1] = targetPoints[i].Y;
            targetFreeXyz[i * 3 + 2] = targetPoints[i].Z;
        }

        var q = new List<double>(network.Graph.Ne);
        foreach (var edge in network.Graph.Edges)
            q.Add(double.IsFinite(edge.Q) ? edge.Q : 1.0);

        try
        {
            var loadNodeIndices = loadNodes.Count > 0
                ? TheseusSolverService.ResolveLoadNodeIndices(network, loadNodes)
                : null;
            var inputs = new SolverInputs
            {
                QInit = q,
                Loads = loads,
                LoadNodeIndices = loadNodeIndices,
            };

            double effectiveRegularization = unconstrainedDirect
                ? _particular switch
                {
                    ParticularMode.MoorePenrose => 0.0,
                    ParticularMode.Tikhonov => regularization,
                    ParticularMode.QrLeastSquares => 0.0,
                    _ => regularization,
                }
                : regularization;
            int particularMethod = InverseFdmUiState.NativeParticularMethod(_particular);

            int nativeMetric = InverseFdmUiState.NativeMetric(_metric);
            int frozenBudget = InverseFdmUiState.FrozenIterationBudget(_metric, _frozenIterations);
            int gnBudget = InverseFdmUiState.GaussNewtonIterationBudget(_metric, _gnIterations);
            var result = TheseusSolverService.SolveInverseFdm(
                network, inputs, targetFreeXyz, effectiveRegularization,
                true, maxL1Iter, particularMethod, (int)_linearAlgebra,
                enforceZeroRx, enforceZeroRy, enforceZeroRz, solveForQ,
                [.. signs], [.. lower], [.. upper], maxIter, tol,
                nativeMetric, gnBudget, cwlsDamping, frozenBudget,
                InverseFdmUiState.NativeStage2Method(_stage2),
                lmDamping, seedGuardMargin, _nondimensionalize, reactionWeight);

            if (_metric == MetricMode.Force
                && _hasBox
                && _linearAlgebra == LinearAlgebraMode.Iterative
                && !result.Converged)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                    $"SPG did not converge within MaxIter={maxIter}.");
            }

            var packedLoads = TheseusSolverService.PackFreeNodeLoads(
                network.FreeNodes.Count, loads, loadNodeIndices);
            var (forces, residuals, ratio) = TargetResidual(
                network, targetPoints, packedLoads, result.ForceDensities);

            var diagnostics = result.InverseDiagnostics ?? new InverseFdmDiagnostics();
            if (_metric == MetricMode.Geometric)
            {
                double loadNorm = Math.Sqrt(packedLoads.Sum(v => v.SquareLength));
                foreach (string warning in InverseFdmUiState.DiagnosticWarnings(diagnostics, loadNorm))
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, warning);
            }

            DA.SetData(0, result.Network);
            DA.SetDataList(1, result.NodePositions);
            DA.SetDataList(2, result.EdgeCurves);
            DA.SetDataList(3, result.ForceDensities);
            DA.SetDataList(4, forces);
            DA.SetDataList(5, residuals);
            DA.SetData(6, ratio);
            DA.SetData(7, result.GeometricError);
            DA.SetData(8, diagnostics.UsedUniformSeed);
            DA.SetDataList(9, InverseFdmUiState.DiagnosticLines(
                diagnostics, result.Iterations, result.Converged, result.GeometricError));

            if (ratio > 0.25 && _metric == MetricMode.Force)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                    $"Large residual ratio ({ratio:0.###}). The target/load combination may be inconsistent with equilibrium; Forces are from the particular at the target, Network is the forward solve.");
            }
        }
        catch (Exception ex)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, ex.Message);
        }
    }

    private static bool TryMapOptionalTree(
        GH_Structure<GH_Number> tree,
        IReadOnlyList<GH_Path> edgePaths,
        string label,
        out double[] values,
        out string? error,
        out string? warning)
    {
        values = [];
        error = null;
        warning = null;
        if (tree.DataCount == 0)
            return true;
        var mapping = QTreeMapper.Map(tree, edgePaths, label);
        if (!mapping.Success)
        {
            error = mapping.Error;
            return false;
        }
        values = mapping.Values!.ToArray();
        warning = mapping.Warning;
        return true;
    }

    private static bool TryMapOptionalTree(
        GH_Structure<GH_Integer> tree,
        IReadOnlyList<GH_Path> edgePaths,
        string label,
        out int[] values,
        out string? error,
        out string? warning)
    {
        values = [];
        error = null;
        warning = null;
        if (tree.DataCount == 0)
            return true;
        var mapping = QTreeMapper.Map(tree, edgePaths, label);
        if (!mapping.Success)
        {
            error = mapping.Error;
            return false;
        }
        values = mapping.Values!.Select(value => (int)Math.Round(value)).ToArray();
        warning = mapping.Warning;
        return true;
    }

    private static (double[] Forces, Vector3d[] Residuals, double Ratio) TargetResidual(
        FDM_Network network, IReadOnlyList<Point3d> target, IReadOnlyList<Vector3d> packedLoads,
        IReadOnlyList<double> q)
    {
        var positions = new Point3d[network.Graph.Nn];
        for (int i = 0; i < network.FreeNodes.Count; i++)
            positions[network.FreeNodes[i]] = target[i];
        for (int i = 0; i < network.FixedNodes.Count; i++)
            positions[network.FixedNodes[i]] = network.Fixed[i].Value;
        var residuals = new Vector3d[network.FreeNodes.Count];
        for (int i = 0; i < residuals.Length; i++)
            residuals[i] = -packedLoads[i];
        var freeLookup = new Dictionary<int, int>();
        for (int i = 0; i < network.FreeNodes.Count; i++) freeLookup[network.FreeNodes[i]] = i;
        var forces = new double[network.Graph.Ne];
        for (int e = 0; e < network.Graph.Ne; e++)
        {
            var edge = network.Graph.Edges[e];
            Vector3d delta = positions[edge.End.Index] - positions[edge.Start.Index];
            forces[e] = q[e] * delta.Length;
            if (freeLookup.TryGetValue(edge.Start.Index, out int start))
                residuals[start] -= q[e] * delta;
            if (freeLookup.TryGetValue(edge.End.Index, out int end))
                residuals[end] += q[e] * delta;
        }
        double residualNorm = 0.0, loadNorm = 0.0;
        for (int i = 0; i < residuals.Length; i++)
        {
            residualNorm += residuals[i].SquareLength;
            loadNorm += packedLoads[i].SquareLength;
        }
        return (forces, residuals, Math.Sqrt(residualNorm) / Math.Max(Math.Sqrt(loadNorm), double.Epsilon));
    }

    protected override void AppendAdditionalComponentMenuItems(ToolStripDropDown menu)
    {
        base.AppendAdditionalComponentMenuItems(menu);
        Menu_AppendSeparator(menu);
        var metricMenu = new ToolStripMenuItem("Metric");
        AppendMetricItem(metricMenu, "Force residual", MetricMode.Force);
        AppendMetricItem(metricMenu, "Geometric residual", MetricMode.Geometric);
        menu.Items.Add(metricMenu);
        Menu_AppendSeparator(menu);
        Menu_AppendItem(menu, "Linear algebra: Direct", (_, _) => SetLinearAlgebra(LinearAlgebraMode.Direct), true, _linearAlgebra == LinearAlgebraMode.Direct);
        Menu_AppendItem(menu, "Linear algebra: Iterative", (_, _) => SetLinearAlgebra(LinearAlgebraMode.Iterative), true, _linearAlgebra == LinearAlgebraMode.Iterative);
        Menu_AppendSeparator(menu);
        var directSolverMenu = new ToolStripMenuItem("Direct solver")
        {
            Enabled = _linearAlgebra == LinearAlgebraMode.Direct,
        };
        AppendDirectSolverItem(directSolverMenu, "Clarabel", ParticularMode.Clarabel);
        AppendDirectSolverItem(directSolverMenu, "Moore–Penrose", ParticularMode.MoorePenrose);
        AppendDirectSolverItem(directSolverMenu, "Tikhonov", ParticularMode.Tikhonov);
        AppendDirectSolverItem(directSolverMenu, "QR least squares", ParticularMode.QrLeastSquares);
        AppendDirectSolverItem(directSolverMenu, "Gram (sparse)", ParticularMode.Gram);
        AppendDirectSolverItem(directSolverMenu, "Gram (dense)", ParticularMode.GramDense);
        menu.Items.Add(directSolverMenu);
        Menu_AppendSeparator(menu);
        var stage2Menu = new ToolStripMenuItem("Stage 2 bounds")
        {
            Enabled = _metric == MetricMode.Geometric,
        };
        AppendStage2Item(stage2Menu, "Active set (BVLS)", Stage2Mode.ActiveSet);
        AppendStage2Item(stage2Menu, "Clarabel (interior point)", Stage2Mode.Clarabel);
        menu.Items.Add(stage2Menu);
        Menu_AppendItem(menu, "Non-dimensionalise", (_, _) => ToggleNondimensionalize(), true, _nondimensionalize);
    }

    private void AppendStage2Item(ToolStripMenuItem parent, string label, Stage2Mode mode)
    {
        var item = new ToolStripMenuItem(label)
        {
            Checked = _stage2 == mode,
        };
        item.Click += (_, _) => SetStage2(mode);
        parent.DropDownItems.Add(item);
    }

    private void SetStage2(Stage2Mode mode)
    {
        if (_stage2 == mode) return;
        RecordUndoEvent("Set Inverse FDM Stage 2 Method");
        _stage2 = mode;
        UpdateMessage();
        ExpireSolution(true);
    }

    private void ToggleNondimensionalize()
    {
        RecordUndoEvent("Toggle Inverse FDM Non-dimensionalisation");
        _nondimensionalize = !_nondimensionalize;
        UpdateMessage();
        ExpireSolution(true);
    }

    private void AppendDirectSolverItem(
        ToolStripMenuItem parent,
        string label,
        ParticularMode mode)
    {
        var item = new ToolStripMenuItem(label)
        {
            Checked = _particular == mode,
        };
        item.Click += (_, _) => SetParticular(mode);
        parent.DropDownItems.Add(item);
    }

    private void AppendMetricItem(ToolStripMenuItem parent, string label, MetricMode mode)
    {
        var item = new ToolStripMenuItem(label)
        {
            Checked = _metric == mode,
        };
        item.Click += (_, _) => SetMetric(mode);
        parent.DropDownItems.Add(item);
    }

    private void SetMetric(MetricMode mode)
    {
        if (_metric == mode) return;
        RecordUndoEvent("Set Inverse FDM Metric");
        _metric = mode;
        UpdateMessage();
        ExpireSolution(true);
    }

    private void SetParticular(ParticularMode mode)
    {
        if (_particular == mode) return;
        RecordUndoEvent("Set Inverse FDM Particular");
        _particular = mode;
        UpdateMessage();
        ExpireSolution(true);
    }

    private void SetLinearAlgebra(LinearAlgebraMode mode)
    {
        if (_linearAlgebra == mode) return;
        RecordUndoEvent("Set Inverse FDM Linear Algebra");
        _linearAlgebra = mode;
        _particular = InverseFdmUiState.UpdateParticular(
            _linearAlgebra, _particular, _hasBox);
        UpdateMessage();
        ExpireSolution(true);
    }

    private void UpdateMessage()
    {
        string unknown = _solveForQ ? "q-init" : "t-init";
        ActiveInverseEngine engine = InverseFdmUiState.ResolveEngine(
            _linearAlgebra, _particular, _hasBox);
        string engineLabel = engine switch
        {
            ActiveInverseEngine.Clarabel => $"Clarabel λ={FormatLambda(_lambda)}",
            ActiveInverseEngine.MoorePenrose => "MP",
            ActiveInverseEngine.Tikhonov => $"Tikh λ={FormatLambda(_lambda)}",
            ActiveInverseEngine.QrLeastSquares => "QR",
            ActiveInverseEngine.Gram => $"Gram λ={FormatLambda(_lambda)}",
            ActiveInverseEngine.GramDense => $"Gram-dense λ={FormatLambda(_lambda)}",
            ActiveInverseEngine.Lsqr when _lambda == 0.0 => "LSQR min-norm",
            ActiveInverseEngine.Lsqr => $"LSQR λ={FormatLambda(_lambda)}",
            _ => $"SPG λ={FormatLambda(_lambda)}",
        };
        string metricLabel = _metric == MetricMode.Geometric
            ? $" · {InverseFdmUiState.PhaseLabel(_frozenIterations, _gnIterations)} λ={FormatLambda(_cwlsDamping)}"
              + $" · {InverseFdmUiState.Stage2Label(_stage2, _seedGuardMargin, _nondimensionalize)}"
            : "";
        Message = $"{engineLabel} · L2 · {unknown}{metricLabel}";
    }

    private static string FormatLambda(double lambda)
    {
        if (lambda == 0.0) return "0";
        return lambda.ToString("0.###e0", CultureInfo.InvariantCulture);
    }

    public override bool Write(GH_IWriter writer)
    {
        writer.SetInt32(ParticularKey, (int)_particular);
        writer.SetInt32(LinearAlgebraKey, (int)_linearAlgebra);
        writer.SetInt32(MetricKey, (int)_metric);
        writer.SetInt32(Stage2Key, (int)_stage2);
        writer.SetBoolean(NondimensionalizeKey, _nondimensionalize);
        return base.Write(writer);
    }

    public override bool Read(GH_IReader reader)
    {
        if (reader.ItemExists(ParticularKey) && Enum.IsDefined(typeof(ParticularMode), reader.GetInt32(ParticularKey)))
            _particular = (ParticularMode)reader.GetInt32(ParticularKey);
        if (reader.ItemExists(LinearAlgebraKey) && Enum.IsDefined(typeof(LinearAlgebraMode), reader.GetInt32(LinearAlgebraKey)))
            _linearAlgebra = (LinearAlgebraMode)reader.GetInt32(LinearAlgebraKey);
        if (reader.ItemExists(MetricKey))
            _metric = reader.GetInt32(MetricKey) == 0 ? MetricMode.Force : MetricMode.Geometric;
        if (reader.ItemExists(Stage2Key))
            _stage2 = reader.GetInt32(Stage2Key) == 1 ? Stage2Mode.Clarabel : Stage2Mode.ActiveSet;
        if (reader.ItemExists(NondimensionalizeKey))
            _nondimensionalize = reader.GetBoolean(NondimensionalizeKey);
        UpdateMessage();
        return base.Read(reader);
    }

    protected override string HtmlHelp_Source() =>
"""
<html>
<body>
<h1>Inverse FDM</h1>
<p>
This component constructs a force-density warm start for a prescribed target
geometry <code>x*</code>, then forward-solves the network with the recovered
<code>q</code>. The pipeline is:
</p>
<p>
<b>Stage-1 particular → seed guard → Frozen CWLS → CWLS-GN → forward FDM.</b>
</p>
<p>
The two CWLS phases are optional and independently budgeted. This is an inverse
equilibrium solve at a frozen geometry, not an algebraic inverse of the square
forward matrix <code>D(q)</code>. Every option of the warm-start study is
reachable from this component; see <i>Reproducing the benchmark</i> below.
</p>

<h2>Stage 1: equilibrium particular</h2>
<p>
At <code>x*</code>, Stage 1 solves a rectangular linear least-squares problem
for one equilibrium particular. <b>SolveQ = false</b> (default) solves member
forces <code>t</code> using target unit directions, then converts
<code>q = t/L*</code>. <b>SolveQ = true</b> solves force densities directly.
The choice affects only this initializer; all CWLS updates use q.
</p>
<p>
The <b>Direct solver</b> menu selects the unboxed particular. Whenever any
Sign, Lower, or Upper is finite, the Direct path is Clarabel regardless of the
menu, since only Clarabel and SPG handle the box exactly.
</p>
<ul>
<li><b>Tikhonov (default when unboxed)</b> — one sparse LDLᵀ of the augmented
saddle <code>[I M; Mᵀ −λI]</code>, solving
<code>min ½‖Mz−p‖² + ½λ‖z‖²</code>; λ must be positive. Same minimiser as
Clarabel without bounds, at a fraction of the cost.</li>
<li><b>Clarabel</b> — convex quadratic least squares with optional q bounds and
reaction equalities. Forced whenever a box is present.</li>
<li><b>Moore–Penrose</b> — the same augmented saddle at λ = 0 for a minimum-norm
unconstrained particular. Fails on rank-deficient systems.</li>
<li><b>QR least squares</b> — sparse QR for tall, full-column-rank systems;
rank-deficient or wide systems report an error.</li>
<li><b>Gram (sparse)</b> — forms <code>MᵀM + λI</code> sparsely and factors it
with sparse LDLᵀ. Algebraically identical to Tikhonov, but the Gram product
squares the condition number and fills in; the saddle factors <code>M</code>
itself. The benchmark's <code>gram_sparse</code>.</li>
<li><b>Gram (dense)</b> — forms <code>MᵀM + λI</code> as a full
<code>ne × ne</code> matrix and factors it with dense Cholesky: O(ne²) memory
and O(ne³) time. Same minimiser as Gram (sparse); kept as the "dense trap"
reference of the talk (<code>gram_dense</code>). Refused above 6000 edges.</li>
</ul>
<p>
With <b>Iterative</b> linear algebra, unconstrained problems use LSQR and
bounded problems use SPG. <b>MaxIter</b> is the per-inner-solve iteration
budget for Clarabel, LSQR, SPG, and the Stage-2 solvers; it is not a CWLS phase
budget.
</p>

<h2>Seed guard</h2>
<p>
Stage 1 minimises the force residual at <code>x*</code>. On unbalanced or
weakly loaded nets a small force residual can hide a collapsed geometry, and
the mechanism: <code>x(q) − x* = −D(q)⁻¹r(q)</code>, a near-singular
<code>D(q)</code> amplifies a tiny <code>r</code>. After Stage 1 the guard
therefore scores a scaled uniform sign seed (one density per sign, scaled so
the reaction magnitude matches the load) on the exact geometric error. When the
clipped Stage-1 seed is worse by more than <b>Guard</b>× both seeds run Stage 2
and the better result continues; <b>Guard Used</b> reports which one. Guard = 0
disables the race (<code>pipeline_noguard</code>, <code>legacy</code>).
</p>

<h2>Metric and compliance weighting</h2>
<p>
For target equilibrium matrix <code>E(x*)</code>, target residual
<code>r(q) = E(x*)q − p</code>, and free-node FDM Laplacian
<code>D(q) = C_fᵀ diag(q) C_f</code> (applied to each coordinate), the exact
constant-load identity is
</p>
<p>
<code>r(q) = D(q)(x* − x(q))</code>, hence
<code>x(q) − x* = −D(q)⁻¹r(q)</code>.
</p>
<ul>
<li><b>Force residual</b> — returns Stage 1 directly and minimises the Euclidean
equilibrium residual. CWLS budgets are ignored.</li>
<li><b>Geometric residual (default)</b> — runs the requested compliance-weighted
phases in q. Set both phase budgets to zero to inspect Stage 1 alone.</li>
</ul>
<p>
All residual objectives in this component are L2. The former L2/L1 toggle and
its IRLS approximation were removed because IRLS was not an exact bounded L1
solve and has no consistent role in compliance-weighted CWLS.
</p>

<h2>Frozen CWLS phase</h2>
<p>
At iteration <code>q_k</code>, Frozen CWLS rebuilds the compliance
<code>D(q_k)⁻¹</code> but keeps the target Jacobian <code>E(x*)</code>. Its
step is the bounded convex least-squares model
</p>
<p>
<code>min_Δq ½‖D(q_k)⁻¹(r(q_k)+E(x*)Δq)‖²
+ ½λcwls‖Δq‖²</code>.
</p>
<p>
This is useful as compliance reweighting of the Stage-1 particular, but it is
not the exact Jacobian of the nonlinear landing map away from the target.
<b>FrozenIter</b> sets its maximum accepted-step attempts; 0 skips it.
</p>

<h2>Gauss–Newton CWLS phase</h2>
<p>
After Frozen CWLS, Gauss–Newton rebuilds both <code>D(q_k)</code> and the
Jacobian <code>E(x(q_k))</code>. For geometry-independent loads,
<code>−D(q_k)⁻¹E(x(q_k))</code> is the derivative of the forward coordinates
with respect to q, so the CWLS-GN step is a true Gauss–Newton model of
<code>½‖x(q)−x*‖²</code>. <b>GNiter</b> sets its maximum; 0 skips it.
</p>
<p>
Both phases backtrack against the exact <b>GeomErr</b>, accept only improving
trials, stop early at Tol, and retain the best point across the whole sequence.
Thus a GN phase cannot replace a better frozen result. λcwls is a fixed q-space
Tikhonov floor on every step; it does not regularize or shift <code>D(q)</code>.
<b>λLM</b> adds Levenberg–Marquardt damping relative to the curvature diagonal
of the Gauss–Newton steps; 0 (default) takes the undamped direction and halves
the step length on the exact merit instead, positive values damp the direction
and grow tenfold on rejected steps. If a Gauss–Newton linearisation point has a
collapsed edge, that step reuses the frozen target Jacobian and the
<b>Diagnostics</b> output counts it as a degenerate linearisation.
</p>

<h2>Stage 2 bounds</h2>
<p>
Each CWLS step is a bounded sparse least-squares problem. The <b>Stage 2
bounds</b> menu picks the bound handling:
</p>
<ul>
<li><b>Active set (BVLS, default)</b> — bounded-variable least squares on the
sparse weighted saddle. A pass fixes the active bounds, factors once, frees or
binds variables from the KKT sign test and repeats; each pass is one LDLᵀ. A
numerical failure falls back to Clarabel for that step (counted in
<b>Diagnostics</b>), and the pass limit returns a feasible partial descent
step (also counted).</li>
<li><b>Clarabel (interior point)</b> — the QP solver from the previous release
on the same step. Same optimum; each interior-point iteration is one
factorisation, so it costs the same as the active set at a few thousand edges
and 2.7–3.9× at 16 k–65 k. This is the <code>legacy</code> Stage 2 of the
benchmark.</li>
</ul>

<h2>Non-dimensionalisation</h2>
<p>
With <b>Non-dimensionalise</b> checked (default), positions are scaled by the
target extent and loads by their norm before anything is assembled, and q,
λ, λcwls and the guard errors are transformed exactly so the same problem is
solved in unit-free coordinates; results are mapped back. This is an exact
change of variables, not a rescaling of the answer. It makes Tol, λcwls,
λLM and the guard margin mean the same thing in metres and millimetres and
improves conditioning of the saddle factors. Uncheck it to reproduce the
benchmark's <code>legacy</code> row or to see the raw-unit behaviour.
</p>

<h2>Phase controls</h2>
<ul>
<li><b>FrozenIter = 1, GNiter = 2</b> — default; the benchmark's
<code>pipeline</code>.</li>
<li><b>FrozenIter = 1, GNiter = 0</b> — one compliance reweight; the
benchmark's <code>frozen</code>.</li>
<li><b>FrozenIter = 0, GNiter = 3</b> — Gauss–Newton only (previous default).</li>
<li><b>FrozenIter = 3, GNiter = 3</b> — frozen warm-up followed by GN.</li>
<li><b>FrozenIter = 0, GNiter = 0</b> — Stage 1 only (<code>s1</code>); the
guard is skipped.</li>
</ul>
<p>
Budgets are nonnegative and have no hard upper cap. Tol stops a phase when
GeomErr is small, or when both relative improvement and relative q-step are
small. Frozen-phase stagnation does not prevent the requested GN phase from
trying its different Jacobian.
</p>
<p>
Migration note: the former experimental <b>L2</b> input slot is now
<b>FrozenIter</b>. Remove any old Boolean wire and supply a nonnegative integer.
</p>

<h2>Reproducing the benchmark</h2>
<p>
The warm-start study (<code>crates/theseus/examples/warm_start_bench</code>)
reports every method as the clipped warm start, then after 1000 L-BFGS-B
iterations on <code>½‖x(q) − x*‖²</code> with direct box bounds. The
component's settings for each method row, with <b>Metric = Geometric</b>
unless noted, SolveQ = false, Rx0/Ry0/Rz0 off, λcwls = 1e-6, MaxIter = 4000,
Tol = 1e-8, and the same Lower / Upper box as the case:
</p>
<table border="1" cellpadding="3" cellspacing="0">
<tr><th>Bench method</th><th>Component settings</th></tr>
<tr><td><code>s1</code></td>
<td>Metric Force (or Geometric with FrozenIter = 0, GNiter = 0); Direct solver
Clarabel forced by the box; λ = 0.</td></tr>
<tr><td><code>gram_sparse</code></td>
<td>No box, SolveQ = true, Direct solver Gram (sparse), FrozenIter = GNiter = 0,
λ = 1e-8 × mean squared entry of <code>E</code>. Clip to the box afterwards.</td></tr>
<tr><td><code>gram_dense</code></td>
<td>As <code>gram_sparse</code> with Direct solver Gram (dense).</td></tr>
<tr><td><code>frozen</code></td>
<td>FrozenIter = 1, GNiter = 0, Stage 2 Active set, Guard = 3,
Non-dimensionalise on.</td></tr>
<tr><td><code>pipeline</code></td>
<td>Defaults: FrozenIter = 1, GNiter = 2, Stage 2 Active set, Guard = 3,
λLM = 0, Non-dimensionalise on, wR = 1.</td></tr>
<tr><td><code>pipeline_noguard</code></td>
<td>As <code>pipeline</code> with Guard = 0.</td></tr>
<tr><td><code>legacy</code></td>
<td>FrozenIter = 1, GNiter = 2, Stage 2 Clarabel, Guard = 0, λLM = 0,
Non-dimensionalise off.</td></tr>
<tr><td><code>uniform</code>, <code>length_ratio</code></td>
<td>Not library solves; feed the sign seed or the heuristic q directly to the
optimiser instead of this component.</td></tr>
</table>
<p>
The benchmark's <code>err/L</code> is this component's <b>GeomErr</b> divided by
the target's bounding-box diagonal. The downstream L-BFGS-B stage is the
regular Ariadne solver with the <b>Force Densities</b> output as the start,
direct box bounds, a target-geometry objective, and a 1000-iteration budget.
The bench's 64 cases are 16 nets × {jittered, bumped target} × {loose, snug
box}; the nets can be rebuilt in Grasshopper from their generator parameters in
<code>warm_start_bench/nets.rs</code>.
</p>

<h2>Bounds, signs, and invertibility</h2>
<p>
<b>Signs</b>, <b>Lower</b>, and <b>Upper</b> always constrain q in every stage.
When Stage 1 solves t, the component maps the q box through
<code>t = L* q</code> using positive target lengths. One value broadcasts;
otherwise data must match the edge tree. Empty channels are unconstrained.
</p>
<p>
Positive q gives a positive-definite D for a connected, properly anchored net.
Mixed-sign and all-compression systems use sparse LDLᵀ and are valid only when
D is nonsingular and numerically factorizable. Near-zero q and sign cancellation
can create mechanisms or failed trial factors. The solver deliberately uses the
exact compliance: no shifted inverse or pseudoinverse is substituted. A
factorisation failure on a trial is not a warning by itself: the step is
rejected and the best point kept. The component warns only on recorded events
(guard swap, degenerate linearisation, Clarabel fallback, active-set cap, or a
realised reaction well above the load along an enforced zero-reaction axis).
</p>

<h2>Other inputs</h2>
<ul>
<li><b>Loads / Load Nodes</b> — without Load Nodes, loads apply to free nodes
in order and the final load repeats. With Load Nodes, one load broadcasts or
one load per listed node is required; unlisted free nodes receive zero.</li>
<li><b>Rx0 / Ry0 / Rz0</b> — add linear zero-reaction rows to the particular
and CWLS subproblems, weighted by <b>wR</b> relative to the equilibrium and
geometric rows. They are least-squares rows, not hard equalities; the realised
reaction along those axes is reported in <b>Diagnostics</b>.</li>
<li><b>Regularization λ</b> — Stage 1 only.</li>
<li><b>CWLS Damping λcwls</b> — both geometric phases only.</li>
<li><b>Guard</b> — seed-guard margin; 0 disables.</li>
<li><b>λLM</b> — Gauss–Newton Levenberg–Marquardt floor; 0 uses step halving.</li>
<li><b>MaxIter</b> — each inner Clarabel, SPG, LSQR, or Stage-2 solve.
Independent of FrozenIter and GNiter.</li>
<li><b>Tol</b> — inner-solver tolerance and geometric phase stopping tolerance.</li>
</ul>

<h2>Outputs</h2>
<p>
<b>Network / Nodes / Edges</b> are the final forward solve.
<b>Forces / Residual / RelRes</b> are evaluated at x* using the returned q.
<b>GeomErr = ‖x(q)−x*‖</b> is a length and directly measures warm-start landing
error. A small force residual ratio does not imply a small GeomErr when the
target is not funicular. <b>Guard</b> is true when the uniform seed won the
race. <b>Diag</b> lists the Stage-2 record as <code>key = value</code> lines:
the Stage-1 and uniform-seed errors, accepted frozen and Gauss–Newton steps,
factorisations, Clarabel fallbacks, active-set cap hits, degenerate
linearisations, and the realised reaction residual.
</p>
</body>
</html>
""";

    protected override Bitmap Icon => Properties.Resources.parameters;

    public override Guid ComponentGuid => new("E1F2A3B4-C5D6-7890-E1F2-A3B4C5D60001");
}

internal enum ParticularMode
{
    MoorePenrose = 0,
    Tikhonov = 1,
    QrLeastSquares = 2,
    /// <summary>Sparse normal equations (LDLᵀ of EᵀE + λI).</summary>
    Gram = 3,
    Clarabel = 4,
    /// <summary>Dense normal equations (dense Cholesky of EᵀE + λI); the "dense trap" reference.</summary>
    GramDense = 5,
}

internal enum LinearAlgebraMode { Direct = 0, Iterative = 1 }

/// <summary>Bound handling for the Stage-2 compliance-weighted steps.</summary>
internal enum Stage2Mode
{
    /// <summary>Active-set BVLS on the sparse weighted saddle (default).</summary>
    ActiveSet = 0,
    /// <summary>Clarabel interior point (the previous default; the talk's "legacy" pipeline).</summary>
    Clarabel = 1,
}

/// <summary>Residual family exposed by the component.</summary>
internal enum MetricMode
{
    /// <summary>Minimize the force residual ‖Mx − p‖.</summary>
    Force = 0,
    /// <summary>Minimize compliance-weighted geometric error.</summary>
    Geometric = 1,
}

internal enum ActiveInverseEngine
{
    Clarabel,
    MoorePenrose,
    Tikhonov,
    QrLeastSquares,
    Gram,
    GramDense,
    Lsqr,
    Spg,
}

internal static class InverseFdmUiState
{
    /// <summary>
    /// Unboxed Stage-1 default. Any finite sign or bound routes Direct to
    /// Clarabel regardless (<see cref="UpdateParticular"/>); without bounds the
    /// augmented saddle is one sparse LDLᵀ where Clarabel is an interior-point
    /// iteration for the same minimiser.
    /// </summary>
    internal const ParticularMode DefaultParticular = ParticularMode.Tikhonov;
    internal const MetricMode DefaultMetric = MetricMode.Geometric;
    /// <summary>Talk/benchmark pipeline: one frozen step, two Gauss–Newton steps.</summary>
    internal const int DefaultFrozenIterations = 1;
    internal const int DefaultGnIterations = 2;
    internal const bool DefaultSolveForQ = false;
    internal const Stage2Mode DefaultStage2 = Stage2Mode.ActiveSet;
    internal const bool DefaultNondimensionalize = true;
    internal const double DefaultSeedGuardMargin = 3.0;
    internal const double DefaultLmDamping = 0.0;
    internal const double DefaultReactionWeight = 1.0;
    /// <summary>Edge count above which the component warns before a dense Gram solve.</summary>
    internal const int DenseGramWarnEdges = 2000;

    internal static int NativeMetric(MetricMode metric) =>
        metric == MetricMode.Force ? 0 : 2;

    internal static InverseStage2Method NativeStage2Method(Stage2Mode mode) =>
        mode == Stage2Mode.Clarabel ? InverseStage2Method.Clarabel : InverseStage2Method.ActiveSet;

    /// <summary>Short Stage-2 label for the component message.</summary>
    internal static string Stage2Label(Stage2Mode stage2, double seedGuardMargin, bool nondimensionalize)
    {
        string label = stage2 == Stage2Mode.Clarabel ? "IP" : "AS";
        if (seedGuardMargin <= 0.0) label += " · guard off";
        if (!nondimensionalize) label += " · dim";
        return label;
    }

    /// <summary>
    /// Event-driven warnings from the Stage-2 diagnostics. Replaces the former
    /// bounds-shape warning, which fired on every mixed-sign net regardless of
    /// what happened.
    /// </summary>
    internal static IEnumerable<string> DiagnosticWarnings(InverseFdmDiagnostics d, double loadNorm)
    {
        if (d.UsedUniformSeed)
            yield return "Seed guard: the Stage-1 particular was replaced by a scaled uniform sign seed "
                + $"(Stage-1 error {d.Stage1Error:0.###e0}, uniform {d.UniformSeedError:0.###e0}).";
        if (d.DegenerateLinearizations > 0)
            yield return $"{d.DegenerateLinearizations} Gauss–Newton step(s) met a collapsed edge at the "
                + "current geometry and used the frozen (target) Jacobian instead.";
        if (d.ClarabelFallbacks > 0)
            yield return $"{d.ClarabelFallbacks} Stage-2 step(s) fell back from the active set to Clarabel "
                + "after a numerical failure of the sparse saddle solve.";
        if (d.ActiveSetCapped > 0)
            yield return $"{d.ActiveSetCapped} Stage-2 step(s) hit the active-set pass limit and returned a "
                + "partial (feasible descent) step; the warm start may be further from the QP optimum.";
        if (d.ReactionResidual > 0.0 && loadNorm > 0.0 && d.ReactionResidual > 0.1 * loadNorm)
            yield return $"Realised reaction along the enforced axes is {d.ReactionResidual:0.###e0} "
                + $"({d.ReactionResidual / loadNorm:0.##}× the load norm): the zero-reaction request is "
                + "inconsistent with the q box or the topology.";
    }

    /// <summary>Human-readable diagnostics lines for the Diag output.</summary>
    internal static IReadOnlyList<string> DiagnosticLines(
        InverseFdmDiagnostics d, int iterations, bool converged, double geometricError)
    {
        return
        [
            $"geom_error = {geometricError:0.######e0}",
            $"iterations = {iterations}",
            $"converged = {converged}",
            $"stage1_error = {d.Stage1Error:0.######e0}",
            $"uniform_seed_error = {d.UniformSeedError:0.######e0}",
            $"used_uniform_seed = {d.UsedUniformSeed}",
            $"frozen_steps = {d.FrozenSteps}",
            $"newton_steps = {d.NewtonSteps}",
            $"stage2_factorizations = {d.Stage2Factorizations}",
            $"clarabel_fallbacks = {d.ClarabelFallbacks}",
            $"active_set_capped = {d.ActiveSetCapped}",
            $"degenerate_linearizations = {d.DegenerateLinearizations}",
            $"reaction_residual = {d.ReactionResidual:0.######e0}",
        ];
    }

    internal static int FrozenIterationBudget(MetricMode metric, int frozenIterations) =>
        metric == MetricMode.Force ? 0 : Math.Max(0, frozenIterations);

    internal static int GaussNewtonIterationBudget(MetricMode metric, int gnIterations) =>
        metric == MetricMode.Force ? 0 : Math.Max(0, gnIterations);

    internal static string PhaseLabel(int frozenIterations, int gnIterations)
    {
        int frozen = Math.Max(0, frozenIterations);
        int gn = Math.Max(0, gnIterations);
        if (frozen == 0 && gn == 0)
            return "Stage 1 only";
        if (frozen == 0)
            return $"GN×{gn}";
        if (gn == 0)
            return $"Frozen×{frozen}";
        return $"Frozen×{frozen} → GN×{gn}";
    }

    internal static bool HasEffectiveBounds(
        IReadOnlyList<int> signs,
        IReadOnlyList<double> lower,
        IReadOnlyList<double> upper) =>
        signs.Any(sign => sign != 0)
        || lower.Any(double.IsFinite)
        || upper.Any(double.IsFinite);

    internal static bool HasStrictSignDefiniteBounds(
        IReadOnlyList<double> lower,
        IReadOnlyList<double> upper)
    {
        bool allPositive = lower.Count > 0 && lower.All(value => value > 0.0);
        bool allNegative = upper.Count > 0 && upper.All(value => value < 0.0);
        return allPositive || allNegative;
    }

    internal static ParticularMode UpdateParticular(
        LinearAlgebraMode linearAlgebra,
        ParticularMode particular,
        bool hasEffectiveBounds) =>
        linearAlgebra == LinearAlgebraMode.Direct && hasEffectiveBounds
            ? ParticularMode.Clarabel
            : particular;

    /// <summary>
    /// Every Stage-1 backend can initialize CWLS because Stage 2 dispatches
    /// independently to a left-weightable backend.
    /// </summary>
    internal static bool SupportsGeometricMetric(
        LinearAlgebraMode linearAlgebra,
        ParticularMode particular,
        bool hasEffectiveBounds)
    {
        _ = linearAlgebra;
        _ = particular;
        _ = hasEffectiveBounds;
        return true;
    }

    internal static int NativeParticularMethod(ParticularMode particular) =>
        particular switch
        {
            ParticularMode.Gram => 0,
            ParticularMode.QrLeastSquares => 2,
            ParticularMode.Clarabel => 3,
            ParticularMode.GramDense => 4,
            _ => 1,
        };

    internal static ActiveInverseEngine ResolveEngine(
        LinearAlgebraMode linearAlgebra,
        ParticularMode particular,
        bool hasEffectiveBounds)
    {
        if (hasEffectiveBounds)
        {
            return linearAlgebra == LinearAlgebraMode.Direct
                ? ActiveInverseEngine.Clarabel
                : ActiveInverseEngine.Spg;
        }

        if (linearAlgebra == LinearAlgebraMode.Iterative)
            return ActiveInverseEngine.Lsqr;

        return particular switch
        {
            ParticularMode.Clarabel => ActiveInverseEngine.Clarabel,
            ParticularMode.MoorePenrose => ActiveInverseEngine.MoorePenrose,
            ParticularMode.Tikhonov => ActiveInverseEngine.Tikhonov,
            ParticularMode.QrLeastSquares => ActiveInverseEngine.QrLeastSquares,
            ParticularMode.GramDense => ActiveInverseEngine.GramDense,
            _ => ActiveInverseEngine.Gram,
        };
    }
}
