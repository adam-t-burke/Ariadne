using System;
using System.Collections.Generic;
using System.Drawing;
using System.Windows.Forms;
using Grasshopper.Kernel;
using GH_IO.Serialization;
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
/// Planar constraint: pull nodes onto a plane along a direction. No target positions —
/// minimizes distance (along the direction) from each node to the plane. Default direction is the plane normal.
/// </summary>
public class ProjectToPlaneComponent : GH_Component
{
    public ProjectToPlaneComponent()
        : base("Project to Plane", "ProjPlane",
            "Pull nodes onto a plane along a direction (default: plane normal). No target positions.",
            "Ariadne", "Objectives")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Nodes", "Nodes", "Nodes to constrain (optional)", GH_ParamAccess.list);
        pManager.AddPlaneParameter("Plane", "P", "Plane to project onto (default: World XY)", GH_ParamAccess.item);
        pManager.AddVectorParameter("Direction", "D", "Projection direction (default: plane normal)", GH_ParamAccess.item);
        pManager.AddNumberParameter("Weight", "Weight", "Objective weight", GH_ParamAccess.item, 1.0);
        pManager[0].Optional = true;
        pManager[1].Optional = true;
        pManager[2].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Objective", "OBJ", "Project to Plane Objective", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Node> nodes = [];
        Plane plane = Plane.WorldXY;
        Vector3d dir = default;
        double weight = 1.0;

        DA.GetDataList(0, nodes);
        bool hasPlane = DA.GetData(1, ref plane);
        bool hasDir = DA.GetData(2, ref dir) && dir.Length > 1e-10;
        DA.GetData(3, ref weight);

        Vector3d? directionOpt = hasDir ? dir : null;
        Plane? planeOpt = hasPlane ? plane : null;
        var objective = new PlanarConstraintAlongDirectionObjective(weight, nodes.Count > 0 ? nodes : null, planeOpt, directionOpt);
        if (!objective.IsValid)
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Direction is parallel to plane; objective will be skipped.");
        DA.SetData(0, objective);
    }

    protected override Bitmap Icon => Properties.Resources.Project_Nodes_To_Plane;
    public override Guid ComponentGuid => new("F8E9B0A1-D2C3-4E5F-6A7B-8C9D0E1F2B4C");
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
    private const string MagnitudeBehaviorKey = "ReactionMagnitudeBehavior";
    private const string MagnitudeSignKey = "ReactionMagnitudeSign";
    private ReactionMagnitudeBehavior _magnitudeBehavior = ReactionMagnitudeBehavior.Max;
    private ReactionMagnitudeSign _magnitudeSign = ReactionMagnitudeSign.Unsigned;

    public ReactionComponent()
        : base("Reaction", "Reaction",
            "Align anchor reaction directions, optionally matching magnitudes.",
            "Ariadne", "Objectives")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Anchor Nodes", "Anchors", "Anchor nodes to constrain", GH_ParamAccess.list);
        pManager.AddVectorParameter("Target Directions", "Dirs", "Target reaction directions and signed magnitude axes", GH_ParamAccess.list);
        pManager.AddBooleanParameter("Include Direction", "Dirs?", "Align reaction directions", GH_ParamAccess.item, true);
        pManager.AddBooleanParameter("Include Magnitude", "Mag?", "Constrain reaction magnitudes", GH_ParamAccess.item, false);
        pManager.AddNumberParameter("Target Magnitudes", "Mags", "Reaction magnitude values (used when Mag? is true)", GH_ParamAccess.list);
        pManager.AddNumberParameter("Weight", "Weight", "Objective weight", GH_ParamAccess.item, 1.0);
        pManager[1].Optional = true;
        pManager[4].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Objective", "OBJ", "Reaction Objective", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Node> anchors = [];
        List<Vector3d> dirs = [];
        bool includeDir = true;
        bool includeMag = false;
        List<double> mags = [];
        double weight = 1.0;

        if (!DA.GetDataList(0, anchors)) return;
        DA.GetDataList(1, dirs);
        DA.GetData(2, ref includeDir);
        DA.GetData(3, ref includeMag);
        DA.GetDataList(4, mags);
        DA.GetData(5, ref weight);

        if (includeDir && dirs.Count == 0)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                "Dirs? is true but no directions were provided; direction mode will be disabled.");
            includeDir = false;
        }

        if (includeMag && mags.Count == 0)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                "Mag? is true but no magnitudes were provided; magnitude mode will be disabled.");
            includeMag = false;
        }

        if (includeMag &&
            _magnitudeSign == ReactionMagnitudeSign.SignedProjected &&
            !HasUsableDirections(dirs, anchors.Count))
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                "Signed projected reaction magnitude requires non-zero directions; magnitude mode will be disabled.");
            includeMag = false;
        }

        if (!includeDir && !includeMag)
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                "Dirs? and Mag? are both inactive; reaction objective will have no effect.");

        var objective = new ReactionObjective(
            weight,
            anchors,
            dirs,
            includeDir,
            includeMag,
            mags.Count > 0 ? mags : null,
            _magnitudeBehavior,
            _magnitudeSign);
        if (!objective.IsValid)
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Reaction objective data is invalid; objective will be skipped.");
        DA.SetData(0, objective);
    }

    protected override void AppendAdditionalComponentMenuItems(ToolStripDropDown menu)
    {
        base.AppendAdditionalComponentMenuItems(menu);
        Menu_AppendSeparator(menu);
        Menu_AppendItem(
            menu,
            "Target Reaction Force",
            (_, _) => SetMagnitudeBehavior(ReactionMagnitudeBehavior.Target),
            true,
            _magnitudeBehavior == ReactionMagnitudeBehavior.Target);
        Menu_AppendItem(
            menu,
            "Max Reaction Force",
            (_, _) => SetMagnitudeBehavior(ReactionMagnitudeBehavior.Max),
            true,
            _magnitudeBehavior == ReactionMagnitudeBehavior.Max);
        Menu_AppendItem(
            menu,
            "Min Reaction Force",
            (_, _) => SetMagnitudeBehavior(ReactionMagnitudeBehavior.Min),
            true,
            _magnitudeBehavior == ReactionMagnitudeBehavior.Min);
        Menu_AppendSeparator(menu);
        Menu_AppendItem(
            menu,
            "Unsigned Magnitude (|R|)",
            (_, _) => SetMagnitudeSign(ReactionMagnitudeSign.Unsigned),
            true,
            _magnitudeSign == ReactionMagnitudeSign.Unsigned);
        Menu_AppendItem(
            menu,
            "Signed Projected Magnitude (dot(R, dir))",
            (_, _) => SetMagnitudeSign(ReactionMagnitudeSign.SignedProjected),
            true,
            _magnitudeSign == ReactionMagnitudeSign.SignedProjected);
    }

    private void SetMagnitudeBehavior(ReactionMagnitudeBehavior behavior)
    {
        if (_magnitudeBehavior == behavior)
            return;

        RecordUndoEvent("Set Reaction Magnitude Behavior");
        _magnitudeBehavior = behavior;
        ExpireSolution(true);
    }

    private void SetMagnitudeSign(ReactionMagnitudeSign sign)
    {
        if (_magnitudeSign == sign)
            return;

        RecordUndoEvent("Set Reaction Magnitude Sign");
        _magnitudeSign = sign;
        ExpireSolution(true);
    }

    public override bool Write(GH_IWriter writer)
    {
        writer.SetInt32(MagnitudeBehaviorKey, (int)_magnitudeBehavior);
        writer.SetInt32(MagnitudeSignKey, (int)_magnitudeSign);
        return base.Write(writer);
    }

    public override bool Read(GH_IReader reader)
    {
        if (reader.ItemExists(MagnitudeBehaviorKey))
        {
            int value = reader.GetInt32(MagnitudeBehaviorKey);
            if (Enum.IsDefined(typeof(ReactionMagnitudeBehavior), value))
                _magnitudeBehavior = (ReactionMagnitudeBehavior)value;
        }
        if (reader.ItemExists(MagnitudeSignKey))
        {
            int value = reader.GetInt32(MagnitudeSignKey);
            if (Enum.IsDefined(typeof(ReactionMagnitudeSign), value))
                _magnitudeSign = (ReactionMagnitudeSign)value;
        }

        return base.Read(reader);
    }

    private static bool HasUsableDirections(List<Vector3d> dirs, int count)
    {
        if (dirs.Count == 0)
            return false;

        for (int i = 0; i < count; i++)
        {
            var dir = i < dirs.Count ? dirs[i] : dirs[^1];
            if (dir.Length < 1e-10)
                return false;
        }
        return true;
    }

    protected override Bitmap Icon => Properties.Resources.Reactions;
    public override Guid ComponentGuid => new("D1E2F3A4-B5C6-7890-A1B2-C3D4E5F60789");
}
