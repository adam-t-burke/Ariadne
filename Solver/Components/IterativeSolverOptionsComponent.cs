using System;
using System.Drawing;
using System.Windows.Forms;
using GH_IO.Serialization;
using Grasshopper.Kernel;
using Theseus.Interop;

namespace Ariadne.Solver.Components;

/// <summary>
/// Produces an <see cref="IterativeSolverOptions"/> object for the
/// "Iterative Options" input of the Optimization Config component. Numeric
/// knobs are inputs; the enumerated knobs (cycle, preconditioner precision,
/// GPU outer loop, GPU adapter preference) live in the right-click menu and
/// are persisted with the definition. Every default matches the native
/// <c>IterativeSolverOptions::default()</c>.
/// </summary>
public class IterativeSolverOptionsComponent : GH_Component
{
    private const string CycleKey = "IterativeCycle";
    private const string PrecisionKey = "IterativePrecision";
    private const string GpuOuterLoopKey = "IterativeGpuOuterLoop";
    private const string AdapterPreferenceKey = "IterativeAdapterPreference";
    private const int PrecisionDefault = -1;

    private MultigridCycle _cycle = MultigridCycle.K;
    private PreconditionerPrecision? _precision;
    private GpuOuterLoop _gpuOuterLoop = GpuOuterLoop.Auto;
    private GpuAdapterPreference _adapterPreference = GpuAdapterPreference.Discrete;

    public IterativeSolverOptionsComponent()
        : base("Iterative Solver Options", "IterOpt",
            "Options for the iterative (AMG-preconditioned) linear solvers selected in the Optimization Config context menu. " +
            "Leave inputs at their defaults unless a solve reports convergence problems.",
            "Ariadne", "Design")
    {
        UpdateMessage();
    }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddNumberParameter("Fixed Tolerance", "Tol",
            "Fixed relative residual for every linear solve. Leave empty for the adaptive schedule driven by Floor, Ceiling and Factor.",
            GH_ParamAccess.item);
        pManager.AddNumberParameter("Tolerance Floor", "Floor",
            "Adaptive schedule: tightest relative residual.",
            GH_ParamAccess.item, IterativeSolverOptions.DefaultToleranceFloor);
        pManager.AddNumberParameter("Tolerance Ceiling", "Ceil",
            "Adaptive schedule: loosest relative residual (used by the first evaluation).",
            GH_ParamAccess.item, IterativeSolverOptions.DefaultToleranceCeiling);
        pManager.AddNumberParameter("Tolerance Factor", "Factor",
            "Adaptive schedule: tol = clamp(Factor × projected-gradient ratio, Floor, Ceiling).",
            GH_ParamAccess.item, IterativeSolverOptions.DefaultToleranceFactor);
        pManager.AddIntegerParameter("Max Iterations", "MaxIter",
            "Iteration budget per linear solve.",
            GH_ParamAccess.item, (int)IterativeSolverOptions.DefaultMaxIterations);
        pManager.AddIntegerParameter("Smoother Degree", "Smooth",
            "Chebyshev smoother degree (1–255).",
            GH_ParamAccess.item, (int)IterativeSolverOptions.DefaultSmootherDegree);
        pManager.AddIntegerParameter("Aggregation Passes", "AggPass",
            "Pairwise matching passes per multigrid level (1–255).",
            GH_ParamAccess.item, (int)IterativeSolverOptions.DefaultAggregationPasses);
        pManager.AddIntegerParameter("Coarsest Size", "Coarse",
            "Stop coarsening once a level has fewer nodes than this.",
            GH_ParamAccess.item, (int)IterativeSolverOptions.DefaultCoarsestSize);
        pManager.AddIntegerParameter("Max Device MB", "DevMB",
            "GPU only: cap on device memory the solver may allocate, in megabytes. Leave empty or 0 for the adapter's own limit.",
            GH_ParamAccess.item);
        pManager[0].Optional = true;
        pManager[8].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Options", "IterOpt", "Iterative linear solver options", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        double fixedTolerance = double.NaN;
        double floor = IterativeSolverOptions.DefaultToleranceFloor;
        double ceiling = IterativeSolverOptions.DefaultToleranceCeiling;
        double factor = IterativeSolverOptions.DefaultToleranceFactor;
        int maxIter = (int)IterativeSolverOptions.DefaultMaxIterations;
        int smoother = (int)IterativeSolverOptions.DefaultSmootherDegree;
        int aggPasses = (int)IterativeSolverOptions.DefaultAggregationPasses;
        int coarsest = (int)IterativeSolverOptions.DefaultCoarsestSize;

        bool hasFixedTolerance = DA.GetData(0, ref fixedTolerance);
        DA.GetData(1, ref floor);
        DA.GetData(2, ref ceiling);
        DA.GetData(3, ref factor);
        DA.GetData(4, ref maxIter);
        DA.GetData(5, ref smoother);
        DA.GetData(6, ref aggPasses);
        DA.GetData(7, ref coarsest);
        int maxDeviceMb = 0;
        DA.GetData(8, ref maxDeviceMb);

        string? error = IterativeSolverOptionsInput.Validate(hasFixedTolerance, fixedTolerance, floor, ceiling, factor, maxIter, smoother, aggPasses, coarsest, maxDeviceMb);
        if (error is not null)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, error);
            return;
        }

        var options = new IterativeSolverOptions
        {
            ToleranceMode = hasFixedTolerance ? IterativeToleranceMode.Fixed : IterativeToleranceMode.Adaptive,
            Tolerance = hasFixedTolerance ? fixedTolerance : IterativeSolverOptions.DefaultFixedTolerance,
            ToleranceFloor = floor,
            ToleranceCeiling = ceiling,
            ToleranceFactor = factor,
            MaxIterations = (uint)maxIter,
            Cycle = _cycle,
            SmootherDegree = (uint)smoother,
            AggregationPasses = (uint)aggPasses,
            CoarsestSize = (uint)coarsest,
            PreconditionPrecision = _precision,
            GpuOuterLoop = _gpuOuterLoop,
            AdapterPreference = _adapterPreference,
            MaxDeviceBytes = IterativeSolverOptionsInput.MegabytesToDeviceBytes(maxDeviceMb),
        };

        DA.SetData(0, options);
    }

    protected override void AppendAdditionalComponentMenuItems(ToolStripDropDown menu)
    {
        base.AppendAdditionalComponentMenuItems(menu);
        Menu_AppendSeparator(menu);

        var cycle = Menu_AppendItem(menu, "Multigrid cycle");
        Menu_AppendItem(cycle.DropDown, "K-cycle (default)", (_, _) => SetCycle(MultigridCycle.K), true, _cycle == MultigridCycle.K);
        Menu_AppendItem(cycle.DropDown, "V-cycle", (_, _) => SetCycle(MultigridCycle.V), true, _cycle == MultigridCycle.V);

        var precision = Menu_AppendItem(menu, "Preconditioner precision");
        Menu_AppendItem(precision.DropDown, "Backend default (f64 CPU, f32 GPU)", (_, _) => SetPrecision(null), true, _precision is null);
        Menu_AppendItem(precision.DropDown, "f64", (_, _) => SetPrecision(PreconditionerPrecision.F64), true, _precision == PreconditionerPrecision.F64);
        Menu_AppendItem(precision.DropDown, "f32", (_, _) => SetPrecision(PreconditionerPrecision.F32), true, _precision == PreconditionerPrecision.F32);

        var outer = Menu_AppendItem(menu, "GPU outer loop");
        Menu_AppendItem(outer.DropDown, "Auto (device when f64 is available)", (_, _) => SetGpuOuterLoop(GpuOuterLoop.Auto), true, _gpuOuterLoop == GpuOuterLoop.Auto);
        Menu_AppendItem(outer.DropDown, "Device", (_, _) => SetGpuOuterLoop(GpuOuterLoop.Device), true, _gpuOuterLoop == GpuOuterLoop.Device);
        Menu_AppendItem(outer.DropDown, "Host", (_, _) => SetGpuOuterLoop(GpuOuterLoop.Host), true, _gpuOuterLoop == GpuOuterLoop.Host);

        var adapter = Menu_AppendItem(menu, "GPU adapter preference");
        Menu_AppendItem(adapter.DropDown, "Discrete first", (_, _) => SetAdapterPreference(GpuAdapterPreference.Discrete), true, _adapterPreference == GpuAdapterPreference.Discrete);
        Menu_AppendItem(adapter.DropDown, "Integrated first", (_, _) => SetAdapterPreference(GpuAdapterPreference.Integrated), true, _adapterPreference == GpuAdapterPreference.Integrated);
        Menu_AppendItem(adapter.DropDown, "Any", (_, _) => SetAdapterPreference(GpuAdapterPreference.Any), true, _adapterPreference == GpuAdapterPreference.Any);
    }

    private void SetCycle(MultigridCycle value)
    {
        if (_cycle == value) return;
        RecordUndoEvent("Set Multigrid Cycle");
        _cycle = value;
        UpdateMessage();
        ExpireSolution(true);
    }

    private void SetPrecision(PreconditionerPrecision? value)
    {
        if (_precision == value) return;
        RecordUndoEvent("Set Preconditioner Precision");
        _precision = value;
        UpdateMessage();
        ExpireSolution(true);
    }

    private void SetGpuOuterLoop(GpuOuterLoop value)
    {
        if (_gpuOuterLoop == value) return;
        RecordUndoEvent("Set GPU Outer Loop");
        _gpuOuterLoop = value;
        ExpireSolution(true);
    }

    private void SetAdapterPreference(GpuAdapterPreference value)
    {
        if (_adapterPreference == value) return;
        RecordUndoEvent("Set GPU Adapter Preference");
        _adapterPreference = value;
        ExpireSolution(true);
    }

    private void UpdateMessage()
    {
        string precision = _precision switch
        {
            PreconditionerPrecision.F64 => ", f64",
            PreconditionerPrecision.F32 => ", f32",
            _ => "",
        };
        Message = $"{_cycle}-cycle{precision}";
    }

    public override bool Write(GH_IWriter writer)
    {
        writer.SetInt32(CycleKey, (int)_cycle);
        writer.SetInt32(PrecisionKey, _precision is { } p ? (int)p : PrecisionDefault);
        writer.SetInt32(GpuOuterLoopKey, (int)_gpuOuterLoop);
        writer.SetInt32(AdapterPreferenceKey, (int)_adapterPreference);
        return base.Write(writer);
    }

    public override bool Read(GH_IReader reader)
    {
        if (reader.ItemExists(CycleKey) && Enum.IsDefined(typeof(MultigridCycle), reader.GetInt32(CycleKey)))
            _cycle = (MultigridCycle)reader.GetInt32(CycleKey);
        if (reader.ItemExists(PrecisionKey))
        {
            int value = reader.GetInt32(PrecisionKey);
            _precision = Enum.IsDefined(typeof(PreconditionerPrecision), value) ? (PreconditionerPrecision)value : null;
        }
        if (reader.ItemExists(GpuOuterLoopKey) && Enum.IsDefined(typeof(GpuOuterLoop), reader.GetInt32(GpuOuterLoopKey)))
            _gpuOuterLoop = (GpuOuterLoop)reader.GetInt32(GpuOuterLoopKey);
        if (reader.ItemExists(AdapterPreferenceKey) && Enum.IsDefined(typeof(GpuAdapterPreference), reader.GetInt32(AdapterPreferenceKey)))
            _adapterPreference = (GpuAdapterPreference)reader.GetInt32(AdapterPreferenceKey);
        UpdateMessage();
        return base.Read(reader);
    }

    protected override Bitmap Icon => Properties.Resources.parameters;

    public override Guid ComponentGuid => new("7C2D9E41-5B3A-4F68-9D1E-2A6B8C4F0A21");
}
