using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;
using Ariadne.Graphs;

namespace Ariadne.Solver.Components;

/// <summary>
/// Minimize variation in edge lengths.
/// </summary>
public class LengthVariationComponent : GH_Component
{
    public LengthVariationComponent()
        : base("Length Variation", "LenVar",
            "Minimize variation in edge lengths.",
            "Ariadne", "Objectives")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Edges", "Edges", "Edges to apply (optional, defaults to all)", GH_ParamAccess.list);
        pManager.AddNumberParameter("Weight", "Weight", "Objective weight", GH_ParamAccess.item, 1.0);
        pManager.AddNumberParameter("Sharpness", "Sharp", "Softmax sharpness parameter", GH_ParamAccess.item, 20.0);
        pManager[0].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Objective", "OBJ", "Length Variation Objective", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Edge> edges = [];
        double weight = 1.0;
        double sharpness = 20.0;

        DA.GetDataList(0, edges);
        DA.GetData(1, ref weight);
        DA.GetData(2, ref sharpness);

        var objective = new LengthVariationObjective(weight, edges.Count > 0 ? edges : null, sharpness);
        DA.SetData(0, objective);
    }

    protected override Bitmap Icon => Properties.Resources.lengthvar;
    public override Guid ComponentGuid => new("B3C4D5E6-F7A8-9012-C3D4-E5F607128901");
}

/// <summary>
/// Minimize variation in member forces.
/// </summary>
public class ForceVariationComponent : GH_Component
{
    public ForceVariationComponent()
        : base("Force Variation", "ForceVar",
            "Minimize variation in member forces.",
            "Ariadne", "Objectives")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Edges", "Edges", "Edges to apply (optional, defaults to all)", GH_ParamAccess.list);
        pManager.AddNumberParameter("Weight", "Weight", "Objective weight", GH_ParamAccess.item, 1.0);
        pManager.AddNumberParameter("Sharpness", "Sharp", "Softmax sharpness parameter", GH_ParamAccess.item, 20.0);
        pManager[0].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Objective", "OBJ", "Force Variation Objective", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Edge> edges = [];
        double weight = 1.0;
        double sharpness = 20.0;

        DA.GetDataList(0, edges);
        DA.GetData(1, ref weight);
        DA.GetData(2, ref sharpness);

        var objective = new ForceVariationObjective(weight, edges.Count > 0 ? edges : null, sharpness);
        DA.SetData(0, objective);
    }

    protected override Bitmap Icon => Properties.Resources.Forcevar;
    public override Guid ComponentGuid => new("C4D5E6F7-A8B9-0123-D4E5-F60712890123");
}

/// <summary>
/// Minimize sum of force x length products (structural efficiency).
/// </summary>
public class PerformanceComponent : GH_Component
{
    public PerformanceComponent()
        : base("Performance", "Perf",
            "Minimize sum of force x length products (structural efficiency).",
            "Ariadne", "Objectives")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Edges", "Edges", "Edges to apply (optional, defaults to all)", GH_ParamAccess.list);
        pManager.AddNumberParameter("Weight", "Weight", "Objective weight", GH_ParamAccess.item, 1.0);
        pManager[0].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Objective", "OBJ", "Performance Objective", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Edge> edges = [];
        double weight = 1.0;

        DA.GetDataList(0, edges);
        DA.GetData(1, ref weight);

        var objective = new PerformanceObjective(weight, edges.Count > 0 ? edges : null);
        DA.SetData(0, objective);
    }

    protected override Bitmap Icon => Properties.Resources.performance;
    public override Guid ComponentGuid => new("D5E6F7A8-B9C0-1234-E5F6-071289012345");
}

/// <summary>
/// Target specific edge lengths.
/// </summary>
public class TargetLengthComponent : GH_Component
{
    public TargetLengthComponent()
        : base("Target Length", "TargetLen",
            "Target specific edge lengths.",
            "Ariadne", "Objectives")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Edges", "Edges", "Edges to apply (optional, defaults to all)", GH_ParamAccess.list);
        pManager.AddNumberParameter("Targets", "Targets", "Target lengths", GH_ParamAccess.list);
        pManager.AddNumberParameter("Weight", "Weight", "Objective weight", GH_ParamAccess.item, 1.0);
        pManager[0].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Objective", "OBJ", "Target Length Objective", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Edge> edges = [];
        List<double> targets = [];
        double weight = 1.0;

        DA.GetDataList(0, edges);
        if (!DA.GetDataList(1, targets)) return;
        DA.GetData(2, ref weight);

        var objective = new TargetLengthObjective(weight, targets, edges.Count > 0 ? edges : null);
        DA.SetData(0, objective);
    }

    protected override Bitmap Icon => Properties.Resources.Target_Length;
    public override Guid ComponentGuid => new("E6F7A8B9-C0D1-2345-F607-128901234567");
}

/// <summary>
/// Minimum length constraint with barrier penalty.
/// </summary>
public class MinLengthComponent : GH_Component
{
    public MinLengthComponent()
        : base("Min Length", "MinLen",
            "Barrier penalty for edges below minimum length.",
            "Ariadne", "Objectives")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Edges", "Edges", "Edges to apply (optional, defaults to all)", GH_ParamAccess.list);
        pManager.AddNumberParameter("Thresholds", "Min", "Minimum length thresholds", GH_ParamAccess.list);
        pManager.AddNumberParameter("Weight", "Weight", "Objective weight", GH_ParamAccess.item, 1.0);
        pManager.AddNumberParameter("Sharpness", "Sharp", "Barrier sharpness", GH_ParamAccess.item, 10.0);
        pManager[0].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Objective", "OBJ", "Min Length Objective", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Edge> edges = [];
        List<double> thresholds = [];
        double weight = 1.0;
        double sharpness = 10.0;

        DA.GetDataList(0, edges);
        if (!DA.GetDataList(1, thresholds)) return;
        DA.GetData(2, ref weight);
        DA.GetData(3, ref sharpness);

        var objective = new MinLengthObjective(weight, thresholds, edges.Count > 0 ? edges : null, sharpness);
        DA.SetData(0, objective);
    }

    protected override Bitmap Icon => Properties.Resources.minlength;
    public override Guid ComponentGuid => new("F7A8B9C0-D1E2-3456-0712-890123456789");
}

/// <summary>
/// Maximum length constraint with barrier penalty.
/// </summary>
public class MaxLengthComponent : GH_Component
{
    public MaxLengthComponent()
        : base("Max Length", "MaxLen",
            "Barrier penalty for edges above maximum length.",
            "Ariadne", "Objectives")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Edges", "Edges", "Edges to apply (optional, defaults to all)", GH_ParamAccess.list);
        pManager.AddNumberParameter("Thresholds", "Max", "Maximum length thresholds", GH_ParamAccess.list);
        pManager.AddNumberParameter("Weight", "Weight", "Objective weight", GH_ParamAccess.item, 1.0);
        pManager.AddNumberParameter("Sharpness", "Sharp", "Barrier sharpness", GH_ParamAccess.item, 10.0);
        pManager[0].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Objective", "OBJ", "Max Length Objective", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Edge> edges = [];
        List<double> thresholds = [];
        double weight = 1.0;
        double sharpness = 10.0;

        DA.GetDataList(0, edges);
        if (!DA.GetDataList(1, thresholds)) return;
        DA.GetData(2, ref weight);
        DA.GetData(3, ref sharpness);

        var objective = new MaxLengthObjective(weight, thresholds, edges.Count > 0 ? edges : null, sharpness);
        DA.SetData(0, objective);
    }

    protected override Bitmap Icon => Properties.Resources.maxlength;
    public override Guid ComponentGuid => new("A8B9C0D1-E2F3-4567-1289-012345678901");
}

/// <summary>
/// Minimum force constraint with barrier penalty.
/// </summary>
public class MinForceComponent : GH_Component
{
    public MinForceComponent()
        : base("Min Force", "MinForce",
            "Barrier penalty for forces below minimum.",
            "Ariadne", "Objectives")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Edges", "Edges", "Edges to apply (optional, defaults to all)", GH_ParamAccess.list);
        pManager.AddNumberParameter("Thresholds", "Min", "Minimum force thresholds", GH_ParamAccess.list);
        pManager.AddNumberParameter("Weight", "Weight", "Objective weight", GH_ParamAccess.item, 1.0);
        pManager.AddNumberParameter("Sharpness", "Sharp", "Barrier sharpness", GH_ParamAccess.item, 10.0);
        pManager[0].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Objective", "OBJ", "Min Force Objective", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Edge> edges = [];
        List<double> thresholds = [];
        double weight = 1.0;
        double sharpness = 10.0;

        DA.GetDataList(0, edges);
        if (!DA.GetDataList(1, thresholds)) return;
        DA.GetData(2, ref weight);
        DA.GetData(3, ref sharpness);

        var objective = new MinForceObjective(weight, thresholds, edges.Count > 0 ? edges : null, sharpness);
        DA.SetData(0, objective);
    }

    protected override Bitmap Icon => Properties.Resources.minforce;
    public override Guid ComponentGuid => new("B9C0D1E2-F3A4-5678-2890-123456789012");
}

/// <summary>
/// Maximum force constraint with barrier penalty.
/// </summary>
public class MaxForceComponent : GH_Component
{
    public MaxForceComponent()
        : base("Max Force", "MaxForce",
            "Barrier penalty for forces above maximum.",
            "Ariadne", "Objectives")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Edges", "Edges", "Edges to apply (optional, defaults to all)", GH_ParamAccess.list);
        pManager.AddNumberParameter("Thresholds", "Max", "Maximum force thresholds", GH_ParamAccess.list);
        pManager.AddNumberParameter("Weight", "Weight", "Objective weight", GH_ParamAccess.item, 1.0);
        pManager.AddNumberParameter("Sharpness", "Sharp", "Barrier sharpness", GH_ParamAccess.item, 10.0);
        pManager[0].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Objective", "OBJ", "Max Force Objective", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Edge> edges = [];
        List<double> thresholds = [];
        double weight = 1.0;
        double sharpness = 10.0;

        DA.GetDataList(0, edges);
        if (!DA.GetDataList(1, thresholds)) return;
        DA.GetData(2, ref weight);
        DA.GetData(3, ref sharpness);

        var objective = new MaxForceObjective(weight, thresholds, edges.Count > 0 ? edges : null, sharpness);
        DA.SetData(0, objective);
    }

    protected override Bitmap Icon => Properties.Resources.maxforce;
    public override Guid ComponentGuid => new("C0D1E2F3-A4B5-6789-3901-234567890123");
}
