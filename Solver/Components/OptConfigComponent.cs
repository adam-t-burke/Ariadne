using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;

namespace Ariadne.Solver.Components;

/// <summary>
/// Grasshopper component that bundles optimization settings into an
/// <see cref="OptimizationConfig"/> object for the Theseus Solve component.
/// </summary>
public class OptConfigComponent : GH_Component
{
    public OptConfigComponent()
        : base("Optimization Config", "OptConfig",
            "Bundle optimization settings for the Theseus solver.",
            "Ariadne", "Design")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Objectives", "OBJ", "Objective functions to minimize", GH_ParamAccess.list);
        pManager.AddNumberParameter("Lower Bounds", "qMin", "Lower bounds on force densities", GH_ParamAccess.list, 0.1);
        pManager.AddNumberParameter("Upper Bounds", "qMax", "Upper bounds on force densities", GH_ParamAccess.list, 100.0);
        pManager.AddIntegerParameter("Max Iterations", "MaxIter", "Maximum solver iterations", GH_ParamAccess.item, 500);
        pManager.AddNumberParameter("Absolute Tolerance", "AbsTol", "Absolute convergence tolerance", GH_ParamAccess.item, 1e-6);
        pManager.AddNumberParameter("Relative Tolerance", "RelTol", "Relative convergence tolerance", GH_ParamAccess.item, 1e-6);
        pManager.AddNumberParameter("Barrier Weight", "BW", "Barrier function weight", GH_ParamAccess.item, 10.0);
        pManager.AddNumberParameter("Barrier Sharpness", "BS", "Barrier function sharpness", GH_ParamAccess.item, 10.0);
        pManager.AddIntegerParameter("Report Frequency", "ReportFreq", "Update viewport preview every N evaluations (0 = no preview)", GH_ParamAccess.item, 10);
        pManager.AddBooleanParameter("Run", "Run", "Toggle true for open-loop optimization; use a button for single-trigger", GH_ParamAccess.item, false);
        pManager.AddBooleanParameter("Stream Preview", "Preview", "Show intermediate geometry during optimization", GH_ParamAccess.item, true);
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Config", "Config", "Optimization configuration", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Objective> objectives = [];
        List<double> lb = [];
        List<double> ub = [];
        int maxIter = 500;
        double absTol = 1e-6;
        double relTol = 1e-6;
        double barrierWeight = 10.0;
        double barrierSharpness = 10.0;
        int reportFreq = 10;
        bool run = false;
        bool streamPreview = true;

        if (!DA.GetDataList(0, objectives)) return;
        DA.GetDataList(1, lb);
        DA.GetDataList(2, ub);
        DA.GetData(3, ref maxIter);
        DA.GetData(4, ref absTol);
        DA.GetData(5, ref relTol);
        DA.GetData(6, ref barrierWeight);
        DA.GetData(7, ref barrierSharpness);
        DA.GetData(8, ref reportFreq);
        DA.GetData(9, ref run);
        DA.GetData(10, ref streamPreview);

        if (objectives.Count == 0)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                "No objectives provided. Optimization will fall back to a forward solve.");
        }

        var config = new OptimizationConfig
        {
            Objectives = objectives.AsReadOnly(),
            LowerBounds = lb.AsReadOnly(),
            UpperBounds = ub.AsReadOnly(),
            MaxIterations = maxIter,
            AbsTol = absTol,
            RelTol = relTol,
            BarrierWeight = barrierWeight,
            BarrierSharpness = barrierSharpness,
            ReportFrequency = reportFreq,
            Run = run,
            StreamPreview = streamPreview,
        };

        DA.SetData(0, config);
    }

    protected override Bitmap Icon => Properties.Resources.parameters;

    public override Guid ComponentGuid => new("A1B2C3D4-E5F6-7890-A1B2-C3D4E5F60001");
}
