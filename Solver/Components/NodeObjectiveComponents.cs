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
/// Minimize deviation from starting nodal positions projected onto an arbitrary plane.
/// Defaults to the world XY plane when no plane is connected.
/// </summary>
public class TargetPlaneComponent : GH_Component
{
    public TargetPlaneComponent()
        : base("Target Plane", "TargetPlane",
            "Minimize deviation from starting positions projected onto a plane (default: World XY).",
            "Ariadne", "Objectives")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Nodes", "Nodes", "Nodes to target (optional)", GH_ParamAccess.list);
        pManager.AddPlaneParameter("Plane", "P", "Target plane (default: World XY)", GH_ParamAccess.item);
        pManager.AddNumberParameter("Weight", "Weight", "Objective weight", GH_ParamAccess.item, 1.0);
        pManager[0].Optional = true;
        pManager[1].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Objective", "OBJ", "Target Plane Objective", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Node> nodes = [];
        Plane plane = Plane.WorldXY;
        double weight = 1.0;

        DA.GetDataList(0, nodes);
        bool hasPlane = DA.GetData(1, ref plane);
        DA.GetData(2, ref weight);

        Plane? planeOpt = hasPlane ? plane : null;
        var objective = new TargetPlaneObjective(weight, nodes.Count > 0 ? nodes : null, planeOpt);
        DA.SetData(0, objective);
    }

    protected override Bitmap Icon => Properties.Resources.Target_UV;
    public override Guid ComponentGuid => new("E7F8A9B0-C1D2-4E5F-6A7B-8C9D0E1F2A3B");
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
/// Align anchor reaction force directions with target directions,
/// optionally also matching target magnitudes.
/// </summary>
public class ReactionComponent : GH_Component
{
    public ReactionComponent()
        : base("Reaction", "Reaction",
            "Align anchor reaction directions, optionally matching magnitudes.",
            "Ariadne", "Objectives")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Anchor Nodes", "Anchors", "Anchor nodes to constrain", GH_ParamAccess.list);
        pManager.AddVectorParameter("Target Directions", "Dirs", "Target reaction directions", GH_ParamAccess.list);
        pManager.AddBooleanParameter("Include Magnitude", "Mag?", "Also match reaction magnitudes", GH_ParamAccess.item, false);
        pManager.AddNumberParameter("Target Magnitudes", "Mags", "Target reaction magnitudes (used when Mag? is true)", GH_ParamAccess.list);
        pManager.AddNumberParameter("Weight", "Weight", "Objective weight", GH_ParamAccess.item, 1.0);
        pManager[3].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Objective", "OBJ", "Reaction Objective", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Node> anchors = [];
        List<Vector3d> dirs = [];
        bool includeMag = false;
        List<double> mags = [];
        double weight = 1.0;

        if (!DA.GetDataList(0, anchors)) return;
        if (!DA.GetDataList(1, dirs)) return;
        DA.GetData(2, ref includeMag);
        DA.GetDataList(3, mags);
        DA.GetData(4, ref weight);

        if (includeMag && mags.Count == 0)
        {
            AddRuntimeMessage(Grasshopper.Kernel.GH_RuntimeMessageLevel.Warning,
                "Include Magnitude is true but no magnitudes provided; direction only will be used.");
            includeMag = false;
        }

        var objective = new ReactionObjective(weight, anchors, dirs, includeMag, mags.Count > 0 ? mags : null);
        DA.SetData(0, objective);
    }

    protected override Bitmap Icon => Properties.Resources.Reactions;
    public override Guid ComponentGuid => new("D1E2F3A4-B5C6-7890-A1B2-C3D4E5F60789");
}
