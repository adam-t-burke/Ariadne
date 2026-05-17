using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;
using Ariadne.Graphs;

namespace Ariadne.Solver.Components;

/// <summary>
/// Defines spherical variable supports around each selected fixed node.
/// Radius is interpreted as a relative bound from each node's initial position.
/// </summary>
public class VariableSupportSphereComponent : GH_Component
{
    public VariableSupportSphereComponent()
        : base("Variable Support Sphere", "VarSphere",
            "Create spherical variable support constraints for fixed nodes.",
            "Ariadne", "Variable Support")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Nodes", "Nodes", "Fixed nodes to make variable", GH_ParamAccess.list);
        pManager.AddNumberParameter("Radius", "R", "Relative radius from initial node position", GH_ParamAccess.item, 1.0);
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Variable Support", "VS", "Sphere variable support definition", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Node> nodes = [];
        double radius = 1.0;
        if (!DA.GetDataList(0, nodes)) return;
        if (!DA.GetData(1, ref radius)) return;

        if (nodes.Count == 0)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No nodes provided.");
            return;
        }
        if (!double.IsFinite(radius) || radius <= 0.0)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Radius must be a positive finite number.");
            return;
        }

        var cfg = new SphereVariableSupport
        {
            Nodes = nodes,
            Radius = radius
        };
        DA.SetData(0, cfg);
    }

    protected override Bitmap Icon => Properties.Resources.Var_Anchor;
    public override Guid ComponentGuid => new("0E6D7A11-6829-42F6-A367-64D5B09E1E01");
}

