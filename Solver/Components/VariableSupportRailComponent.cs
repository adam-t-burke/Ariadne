using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;
using Rhino.Geometry;
using Ariadne.Graphs;

namespace Ariadne.Solver.Components;

/// <summary>
/// Defines line-segment constrained variable supports for fixed nodes.
/// </summary>
public class VariableSupportRailComponent : GH_Component
{
    public VariableSupportRailComponent()
        : base("Variable Support Rail", "VarRail",
            "Create rail-constrained variable supports along an input line segment.",
            "Ariadne", "Variable Support")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Nodes", "Nodes", "Fixed nodes to make variable", GH_ParamAccess.list);
        pManager.AddLineParameter("Rail", "Rail", "Line segment rail constraint", GH_ParamAccess.item);
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Variable Support", "VS", "Rail variable support definition", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Node> nodes = [];
        Line rail = default;
        if (!DA.GetDataList(0, nodes)) return;
        if (!DA.GetData(1, ref rail)) return;

        if (nodes.Count == 0)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No nodes provided.");
            return;
        }
        if (!rail.IsValid || rail.Length <= 1e-12)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Rail must be a valid non-degenerate line segment.");
            return;
        }

        var cfg = new RailVariableSupport
        {
            Nodes = nodes,
            Rail = rail
        };
        DA.SetData(0, cfg);
    }

    protected override Bitmap Icon => Properties.Resources.Var_Anchor;
    public override Guid ComponentGuid => new("E754A445-83F3-430D-A5DD-0AC8234A7503");
}

