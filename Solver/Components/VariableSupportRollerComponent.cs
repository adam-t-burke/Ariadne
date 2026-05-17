using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;
using Rhino.Geometry;
using Ariadne.Graphs;

namespace Ariadne.Solver.Components;

/// <summary>
/// Defines roller-style variable supports for fixed nodes with per-axis relative domains.
/// </summary>
public class VariableSupportRollerComponent : GH_Component
{
    public VariableSupportRollerComponent()
        : base("Variable Support Roller", "VarRoller",
            "Create roller-style variable supports with per-axis relative bounds from initial position.",
            "Ariadne", "Variable Support")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Nodes", "Nodes", "Fixed nodes to make variable", GH_ParamAccess.list);
        pManager.AddBooleanParameter("Free X", "Fx", "Allow movement in X", GH_ParamAccess.item, true);
        pManager.AddBooleanParameter("Free Y", "Fy", "Allow movement in Y", GH_ParamAccess.item, true);
        pManager.AddBooleanParameter("Free Z", "Fz", "Allow movement in Z", GH_ParamAccess.item, true);
        pManager.AddIntervalParameter("X Domain", "Dx", "Relative X bounds from initial position", GH_ParamAccess.item, new Interval(-1.0, 1.0));
        pManager.AddIntervalParameter("Y Domain", "Dy", "Relative Y bounds from initial position", GH_ParamAccess.item, new Interval(-1.0, 1.0));
        pManager.AddIntervalParameter("Z Domain", "Dz", "Relative Z bounds from initial position", GH_ParamAccess.item, new Interval(-1.0, 1.0));
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Variable Support", "VS", "Roller variable support definition", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Node> nodes = [];
        bool freeX = true, freeY = true, freeZ = true;
        Interval dx = new(-1.0, 1.0), dy = new(-1.0, 1.0), dz = new(-1.0, 1.0);

        if (!DA.GetDataList(0, nodes)) return;
        DA.GetData(1, ref freeX);
        DA.GetData(2, ref freeY);
        DA.GetData(3, ref freeZ);
        DA.GetData(4, ref dx);
        DA.GetData(5, ref dy);
        DA.GetData(6, ref dz);

        if (nodes.Count == 0)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No nodes provided.");
            return;
        }
        if (!(freeX || freeY || freeZ))
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "At least one axis must be free.");
            return;
        }

        if ((freeX && !dx.IsIncreasing) || (freeY && !dy.IsIncreasing) || (freeZ && !dz.IsIncreasing))
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Enabled axis domains must be increasing.");
            return;
        }

        var cfg = new RollerVariableSupport
        {
            Nodes = nodes,
            FreeX = freeX,
            FreeY = freeY,
            FreeZ = freeZ,
            DomainX = dx,
            DomainY = dy,
            DomainZ = dz
        };
        DA.SetData(0, cfg);
    }

    protected override Bitmap Icon => Properties.Resources.Var_Anchor;
    public override Guid ComponentGuid => new("27AE0E9F-F57C-47B8-A312-3B665ECF4E02");
}

