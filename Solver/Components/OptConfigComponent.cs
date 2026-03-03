using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Windows.Forms;
using GH_IO.Serialization;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Parameters;
using Grasshopper.Kernel.Types;
using Objective = Ariadne.Solver.Objective;

namespace Ariadne.Solver.Components;

/// <summary>
/// Grasshopper component that bundles optimization settings into an
/// <see cref="OptimizationConfig"/> object for the Theseus Solve component.
/// Objectives are always flattened so a single solve runs (avoids multiple concurrent solves from tree branches).
/// </summary>
public class OptConfigComponent : GH_Component
{
    private const string QParameterizationModeKey = "QParameterizationMode";
    private const int CoreInputCount = 10;
    private QParameterizationMode _qParameterizationMode = QParameterizationMode.DirectBoxBounds;
    private bool _hasLegacyImplicitBoundedMode;

    public OptConfigComponent()
        : base("Optimization Config", "OptConfig",
            "Bundle optimization settings for the Theseus solver.",
            "Ariadne", "Design")
    {
        UpdateMessage();
    }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Objectives", "OBJ", "Objective functions to minimize (flattened automatically)", GH_ParamAccess.tree);
        pManager.AddGenericParameter("Variable Supports", "VS", "Variable support definitions (optional, flattened)", GH_ParamAccess.tree);
        pManager.AddNumberParameter("Lower Bounds", "qMin",
            "Lower force-density bounds. Match/graft to the edge tree; one value broadcasts globally or within its branch.",
            GH_ParamAccess.tree, 0.1);
        pManager.AddNumberParameter("Upper Bounds", "qMax",
            "Upper force-density bounds. Match/graft to the edge tree; one value broadcasts globally or within its branch.",
            GH_ParamAccess.tree, 100.0);
        pManager.AddIntegerParameter("Max Iterations", "MaxIter",
            "Maximum accepted optimization iterations. Positive values are applied directly, with no hidden upper cap.",
            GH_ParamAccess.item, 500);
        pManager.AddNumberParameter("Absolute Tolerance", "AbsTol",
            "Projected-gradient infinity-norm tolerance in Direct Box Bounds. Direct Soft Bounds uses this as its L-BFGS gradient tolerance.",
            GH_ParamAccess.item, 1e-6);
        pManager.AddNumberParameter("Relative Tolerance", "RelTol",
            "Relative accepted-iterate reduction tolerance in Direct Box Bounds: stop when (f[k]-f[k+1])/max(|f[k]|,|f[k+1]|,1) <= RelTol. Direct Soft Bounds uses this as its relative cost tolerance.",
            GH_ParamAccess.item, 1e-6);
        pManager.AddIntegerParameter("Report Frequency", "ReportFreq", "Invoke progress callback every N accepted L-BFGS iterations (0 = every iteration)", GH_ParamAccess.item, 10);
        pManager.AddBooleanParameter("Run", "Run", "Toggle true for open-loop optimization; use a button for single-trigger", GH_ParamAccess.item, false);
        pManager.AddBooleanParameter("Stream Preview", "Stream", "Stream intermediate results to outputs during optimization (false = only output final result)", GH_ParamAccess.item, true);
        pManager[1].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Config", "Config", "Optimization configuration", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Objective> objectives = [];
        List<VariableSupportConfig> variableSupports = [];
        var lbTree = new GH_Structure<GH_Number>();
        var ubTree = new GH_Structure<GH_Number>();
        int maxIter = 500;
        double absTol = 1e-6;
        double relTol = 1e-6;
        double barrierWeight = 10.0;
        double barrierSharpness = 10.0;
        int reportFreq = 10;
        bool run = false;
        bool streamPreview = true;

        // Flatten objectives tree so we always get one list — avoids multiple concurrent solves when branches are not flattened
        var objTree = new GH_Structure<IGH_Goo>();
        if (!DA.GetDataTree(0, out objTree)) return;
        objTree.Flatten();
        objectives.Clear();
        foreach (var branch in objTree.Branches)
        {
            foreach (var goo in branch)
            {
                if (goo?.ScriptVariable() is Objective obj)
                    objectives.Add(obj);
            }
        }
        var vsTree = new GH_Structure<IGH_Goo>();
        if (DA.GetDataTree(1, out vsTree))
        {
            vsTree.Flatten();
            variableSupports.Clear();
            foreach (var branch in vsTree.Branches)
            {
                foreach (var goo in branch)
                {
                    if (goo?.ScriptVariable() is VariableSupportConfig vs)
                        variableSupports.Add(vs);
                }
            }
        }

        if (!DA.GetDataTree(2, out lbTree)) return;
        if (!DA.GetDataTree(3, out ubTree)) return;
        DA.GetData(4, ref maxIter);
        DA.GetData(5, ref absTol);
        DA.GetData(6, ref relTol);
        DA.GetData(7, ref reportFreq);
        DA.GetData(8, ref run);
        DA.GetData(9, ref streamPreview);
        int barrierWeightIndex = Params.IndexOfInputParam("Barrier Weight");
        int barrierSharpnessIndex = Params.IndexOfInputParam("Barrier Sharpness");
        if (barrierWeightIndex >= 0)
            DA.GetData(barrierWeightIndex, ref barrierWeight);
        if (barrierSharpnessIndex >= 0)
            DA.GetData(barrierSharpnessIndex, ref barrierSharpness);

        if (objectives.Count == 0)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                "No objectives provided. Optimization will fall back to a forward solve.");
        }

        var lowerBounds = EdgeValueTree.FromGh(lbTree);
        var upperBounds = EdgeValueTree.FromGh(ubTree);
        ResolveLegacyImplicitBoundedMode(
            lowerBounds.FlattenedValues.ToArray(),
            upperBounds.FlattenedValues.ToArray());

        var config = new OptimizationConfig
        {
            Objectives = objectives.AsReadOnly(),
            LowerBounds = lowerBounds,
            UpperBounds = upperBounds,
            MaxIterations = maxIter,
            AbsTol = absTol,
            RelTol = relTol,
            BarrierWeight = barrierWeight,
            BarrierSharpness = barrierSharpness,
            ReportFrequency = reportFreq,
            QParameterizationMode = _qParameterizationMode,
            Run = run,
            StreamPreview = streamPreview,
            VariableSupports = variableSupports.AsReadOnly(),
        };

        DA.SetData(0, config);
    }

    protected override void AppendAdditionalComponentMenuItems(ToolStripDropDown menu)
    {
        base.AppendAdditionalComponentMenuItems(menu);
        Menu_AppendSeparator(menu);
        Menu_AppendItem(
            menu,
            "q Mode: Direct Box Bounds (L-BFGS-B)",
            (_, _) => SetQParameterizationMode(QParameterizationMode.DirectBoxBounds),
            true,
            _qParameterizationMode == QParameterizationMode.DirectBoxBounds);
        Menu_AppendItem(
            menu,
            "q Mode: Direct Soft Bounds",
            (_, _) => SetQParameterizationMode(QParameterizationMode.DirectSoftBounds),
            true,
            _qParameterizationMode == QParameterizationMode.DirectSoftBounds);
    }

    private void SetQParameterizationMode(QParameterizationMode mode)
    {
        if (_qParameterizationMode == mode)
            return;

        RecordUndoEvent("Set q Parameterization Mode");
        _hasLegacyImplicitBoundedMode = false;
        _qParameterizationMode = mode;
        UpdateParameterVisibility();
        UpdateMessage();
        ExpireSolution(true);
    }

    private void UpdateMessage()
    {
        Message = _qParameterizationMode switch
        {
            QParameterizationMode.DirectSoftBounds => "q: SoftBounds",
            QParameterizationMode.DirectBoxBounds => "q: BoxBounds",
            _ => "q: SoftBounds",
        };
    }

    private void UpdateParameterVisibility()
    {
        if (Params.Input.Count < CoreInputCount)
            return;

        for (int i = Params.Input.Count - 1; i >= CoreInputCount; i--)
            Params.UnregisterInputParameter(Params.Input[i], true);

        if (_qParameterizationMode == QParameterizationMode.DirectSoftBounds)
        {
            Params.RegisterInputParam(NumberParam(
                "Barrier Weight", "BW", "Soft-bound barrier function weight"));
            Params.RegisterInputParam(NumberParam(
                "Barrier Sharpness", "BS", "Soft-bound barrier function sharpness"));
        }

        Params.OnParametersChanged();
    }

    private static Param_Number NumberParam(string name, string nickname, string description)
    {
        var param = new Param_Number
        {
            Name = name,
            NickName = nickname,
            Description = description,
            Access = GH_ParamAccess.item,
            Optional = true,
        };
        param.SetPersistentData(10.0);
        return param;
    }

    public override bool Write(GH_IWriter writer)
    {
        writer.SetInt32(QParameterizationModeKey, (int)_qParameterizationMode);
        return base.Write(writer);
    }

    public override bool Read(GH_IReader reader)
    {
        if (reader.ItemExists(QParameterizationModeKey))
        {
            int value = reader.GetInt32(QParameterizationModeKey);
            if (value == 1)
            {
                _hasLegacyImplicitBoundedMode = true;
                _qParameterizationMode = QParameterizationMode.DirectSoftBounds;
            }
            else if (Enum.IsDefined(typeof(QParameterizationMode), value))
            {
                _qParameterizationMode = (QParameterizationMode)value;
            }
        }
        UpdateParameterVisibility();
        UpdateMessage();
        return base.Read(reader);
    }

    private void ResolveLegacyImplicitBoundedMode(IReadOnlyList<double> lb, IReadOnlyList<double> ub)
    {
        if (!_hasLegacyImplicitBoundedMode)
            return;

        _qParameterizationMode = HasFiniteTwoSidedBounds(lb, ub)
            ? QParameterizationMode.DirectBoxBounds
            : QParameterizationMode.DirectSoftBounds;
        _hasLegacyImplicitBoundedMode = false;
        UpdateParameterVisibility();
        UpdateMessage();
        AddRuntimeMessage(
            GH_RuntimeMessageLevel.Remark,
            $"Migrated legacy q: ImplicitBounds mode to {_qParameterizationMode}. " +
            "Use BoxBounds for finite two-sided q bounds, or SoftBounds for one-sided/infinite bounds.");
    }

    private static bool HasFiniteTwoSidedBounds(IReadOnlyList<double> lb, IReadOnlyList<double> ub)
    {
        int n = Math.Min(lb.Count, ub.Count);
        if (n == 0)
            return false;

        for (int i = 0; i < n; i++)
        {
            if (!double.IsFinite(lb[i]) || !double.IsFinite(ub[i]) || ub[i] <= lb[i])
                return false;
        }
        return true;
    }

    protected override Bitmap Icon => Properties.Resources.parameters;

    public override Guid ComponentGuid => new("A1B2C3D4-E5F6-7890-A1B2-C3D4E5F60001");
}
