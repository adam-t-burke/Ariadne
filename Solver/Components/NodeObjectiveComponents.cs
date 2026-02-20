using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;
using Rhino.Geometry;
using Ariadne.Graphs;

namespace Ariadne.Solver.Components;

/// <summary>
/// Minimize XY deviation from starting nodal positions (TNA-like behavior).
/// </summary>
public class TargetXYComponent : GH_Component
{
    public TargetXYComponent()
        : base("Target Geometry XY", "TargetXY",
            "Minimize deviation from starting nodal positions in the XY plane (TNA-like).",
            "Ariadne", "Objectives")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Nodes", "Nodes", "Nodes to target (optional)", GH_ParamAccess.list);
        pManager.AddNumberParameter("Weight", "Weight", "Objective weight", GH_ParamAccess.item, 1.0);
        pManager[0].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Objective", "OBJ", "Target XY Objective", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Node> nodes = [];
        double weight = 1.0;

        DA.GetDataList(0, nodes);
        DA.GetData(1, ref weight);

        var objective = new TargetXYObjective(weight, nodes.Count > 0 ? nodes : null);
        DA.SetData(0, objective);
    }

    protected override Bitmap Icon => Properties.Resources.Target_XY;
    public override Guid ComponentGuid => new("B8F7A2C1-3D4E-5F60-9182-A3B4C5D6E7F8");
}

/// <summary>
/// Minimize 3D deviation from starting nodal positions.
/// </summary>
public class TargetXYZComponent : GH_Component
{
    public TargetXYZComponent()
        : base("Target Geometry XYZ", "TargetXYZ",
            "Minimize deviation from starting nodal positions in 3D.",
            "Ariadne", "Objectives")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Nodes", "Nodes", "Nodes to target (optional)", GH_ParamAccess.list);
        pManager.AddNumberParameter("Weight", "Weight", "Objective weight", GH_ParamAccess.item, 1.0);
        pManager[0].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Objective", "OBJ", "Target XYZ Objective", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Node> nodes = [];
        double weight = 1.0;

        DA.GetDataList(0, nodes);
        DA.GetData(1, ref weight);

        var objective = new TargetXYZObjective(weight, nodes.Count > 0 ? nodes : null);
        DA.SetData(0, objective);
    }

    protected override Bitmap Icon => Properties.Resources.Target;
    public override Guid ComponentGuid => new("A2B3C4D5-E6F7-8901-B2C3-D4E5F6071234");
}

/// <summary>
/// Maintain relative distances within a point set (rigid body).
/// </summary>
public class RigidPointSetComponent : GH_Component
{
    public RigidPointSetComponent()
        : base("Rigid Point Set", "RigidSet",
            "Maintain relative distances within a point set.",
            "Ariadne", "Objectives")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Nodes", "Nodes", "Nodes to keep as rigid set", GH_ParamAccess.list);
        pManager.AddNumberParameter("Weight", "Weight", "Objective weight", GH_ParamAccess.item, 1.0);
        pManager[0].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Objective", "OBJ", "Rigid Point Set Objective", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Node> nodes = [];
        double weight = 1.0;

        DA.GetDataList(0, nodes);
        DA.GetData(1, ref weight);

        var objective = new RigidPointSetObjective(weight, nodes.Count > 0 ? nodes : null);
        DA.SetData(0, objective);
    }

    protected override Bitmap Icon => Properties.Resources.PointSet;
    public override Guid ComponentGuid => new("C9E8B7A6-4D5C-6E7F-8091-F2E3D4C5B6A7");
}

/// <summary>
/// Align anchor reaction force directions with target directions.
/// </summary>
public class ReactionDirectionComponent : GH_Component
{
    public ReactionDirectionComponent()
        : base("Reaction Direction", "ReactDir",
            "Align anchor reaction force directions with target directions.",
            "Ariadne", "Objectives")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Anchor Nodes", "Anchors", "Anchor nodes to constrain", GH_ParamAccess.list);
        pManager.AddVectorParameter("Target Directions", "Dirs", "Target reaction directions", GH_ParamAccess.list);
        pManager.AddNumberParameter("Weight", "Weight", "Objective weight", GH_ParamAccess.item, 1.0);
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Objective", "OBJ", "Reaction Direction Objective", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Node> anchors = [];
        List<Vector3d> dirs = [];
        double weight = 1.0;

        if (!DA.GetDataList(0, anchors)) return;
        if (!DA.GetDataList(1, dirs)) return;
        DA.GetData(2, ref weight);

        var objective = new ReactionDirectionObjective(weight, anchors, dirs);
        DA.SetData(0, objective);
    }

    protected override Bitmap? Icon => null;
    public override Guid ComponentGuid => new("D1E2F3A4-B5C6-7890-A1B2-C3D4E5F60789");
}

/// <summary>
/// Align anchor reaction force directions and match target magnitudes.
/// </summary>
public class ReactionDirectionMagnitudeComponent : GH_Component
{
    public ReactionDirectionMagnitudeComponent()
        : base("Reaction Direction + Magnitude", "ReactDirMag",
            "Align reaction directions and match target magnitudes at anchors.",
            "Ariadne", "Objectives")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Anchor Nodes", "Anchors", "Anchor nodes to constrain", GH_ParamAccess.list);
        pManager.AddVectorParameter("Target Directions", "Dirs", "Target reaction directions", GH_ParamAccess.list);
        pManager.AddNumberParameter("Target Magnitudes", "Mags", "Target reaction magnitudes", GH_ParamAccess.list);
        pManager.AddNumberParameter("Weight", "Weight", "Objective weight", GH_ParamAccess.item, 1.0);
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Objective", "OBJ", "Reaction Direction + Magnitude Objective", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Node> anchors = [];
        List<Vector3d> dirs = [];
        List<double> mags = [];
        double weight = 1.0;

        if (!DA.GetDataList(0, anchors)) return;
        if (!DA.GetDataList(1, dirs)) return;
        if (!DA.GetDataList(2, mags)) return;
        DA.GetData(3, ref weight);

        var objective = new ReactionDirectionMagnitudeObjective(weight, anchors, dirs, mags);
        DA.SetData(0, objective);
    }

    protected override Bitmap? Icon => null;
    public override Guid ComponentGuid => new("E2F3A4B5-C6D7-8901-B2C3-D4E5F6078901");
}
