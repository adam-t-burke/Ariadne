using System;
using System.Collections.Generic;
using System.Drawing;
using System.Windows.Forms;
using GH_IO.Serialization;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Parameters;
using Rhino.Geometry;
using Ariadne.Graphs;

namespace Ariadne.Solver.Components;

public enum VariableSupportMode
{
    Sphere = 0,
    Roller = 1,
    Rail = 2,
    NurbsCurve = 3,
    NurbsSurface = 4,
}

/// <summary>
/// Defines variable supports for fixed nodes using a right-click selectable mode.
/// </summary>
public class VariableSupportComponent : GH_Component
{
    private const string ModeKey = "VariableSupportMode";

    private const int IdxNodes = 0;
    private VariableSupportMode _mode = VariableSupportMode.Sphere;

    public VariableSupportComponent()
        : base("Variable Support", "VarSupport",
            "Create variable support constraints for fixed nodes.",
            "Ariadne", "Variable Support")
    {
        UpdateMessage();
    }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Nodes", "Nodes", "Fixed nodes to make variable", GH_ParamAccess.list);
        pManager.AddNumberParameter("Radius", "R", "Relative radius from initial node position", GH_ParamAccess.item, 1.0);
        pManager[1].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Variable Support", "VS", "Variable support definition", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        List<Node> nodes = [];
        if (!DA.GetDataList(IdxNodes, nodes)) return;
        if (nodes.Count == 0)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No nodes provided.");
            return;
        }

        switch (_mode)
        {
            case VariableSupportMode.Sphere:
                SolveSphere(DA, nodes);
                break;
            case VariableSupportMode.Roller:
                SolveRoller(DA, nodes);
                break;
            case VariableSupportMode.Rail:
                SolveRail(DA, nodes);
                break;
            case VariableSupportMode.NurbsCurve:
                SolveCurve(DA, nodes);
                break;
            case VariableSupportMode.NurbsSurface:
                SolveSurface(DA, nodes);
                break;
        }
    }

    private void SolveSphere(IGH_DataAccess DA, List<Node> nodes)
    {
        double radius = 1.0;
        DA.GetData(1, ref radius);
        if (!double.IsFinite(radius) || radius <= 0.0)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Radius must be a positive finite number.");
            return;
        }
        DA.SetData(0, new SphereVariableSupport { Nodes = nodes, Radius = radius });
    }

    private void SolveRoller(IGH_DataAccess DA, List<Node> nodes)
    {
        bool freeX = true, freeY = true, freeZ = true;
        Interval dx = new(-1.0, 1.0), dy = new(-1.0, 1.0), dz = new(-1.0, 1.0);
        DA.GetData(1, ref freeX);
        DA.GetData(2, ref freeY);
        DA.GetData(3, ref freeZ);
        DA.GetData(4, ref dx);
        DA.GetData(5, ref dy);
        DA.GetData(6, ref dz);

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

        DA.SetData(0, new RollerVariableSupport
        {
            Nodes = nodes,
            FreeX = freeX,
            FreeY = freeY,
            FreeZ = freeZ,
            DomainX = dx,
            DomainY = dy,
            DomainZ = dz
        });
    }

    private void SolveRail(IGH_DataAccess DA, List<Node> nodes)
    {
        Line rail = default;
        if (!DA.GetData(1, ref rail)) return;
        if (!rail.IsValid || rail.Length <= 1e-12)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Rail must be a valid non-degenerate line segment.");
            return;
        }
        DA.SetData(0, new RailVariableSupport { Nodes = nodes, Rail = rail });
    }

    private void SolveCurve(IGH_DataAccess DA, List<Node> nodes)
    {
        Curve? curve = null;
        Interval domain = Interval.Unset;
        if (!DA.GetData(1, ref curve) || curve == null) return;
        bool hasDomain = DA.GetData(2, ref domain);
        if (!curve.IsValid)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Curve must be valid.");
            return;
        }
        if (hasDomain && !domain.IsIncreasing)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Curve domain must be increasing.");
            return;
        }
        DA.SetData(0, new CurveVariableSupport
        {
            Nodes = nodes,
            Curve = curve.DuplicateCurve(),
            Domain = hasDomain ? domain : null
        });
    }

    private void SolveSurface(IGH_DataAccess DA, List<Node> nodes)
    {
        Surface? surface = null;
        Interval domainU = Interval.Unset;
        Interval domainV = Interval.Unset;
        if (!DA.GetData(1, ref surface) || surface == null) return;
        bool hasDomainU = DA.GetData(2, ref domainU);
        bool hasDomainV = DA.GetData(3, ref domainV);
        if (!surface.IsValid)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Surface must be valid.");
            return;
        }
        if ((hasDomainU && !domainU.IsIncreasing) || (hasDomainV && !domainV.IsIncreasing))
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Surface domains must be increasing.");
            return;
        }
        DA.SetData(0, new SurfaceVariableSupport
        {
            Nodes = nodes,
            Surface = surface.ToNurbsSurface(),
            DomainU = hasDomainU ? domainU : null,
            DomainV = hasDomainV ? domainV : null
        });
    }

    protected override void AppendAdditionalComponentMenuItems(ToolStripDropDown menu)
    {
        base.AppendAdditionalComponentMenuItems(menu);
        Menu_AppendSeparator(menu);
        AppendModeItem(menu, "Sphere", VariableSupportMode.Sphere);
        AppendModeItem(menu, "Roller", VariableSupportMode.Roller);
        AppendModeItem(menu, "Rail", VariableSupportMode.Rail);
        AppendModeItem(menu, "NURBS Curve", VariableSupportMode.NurbsCurve);
        AppendModeItem(menu, "NURBS Surface", VariableSupportMode.NurbsSurface);
    }

    private void AppendModeItem(ToolStripDropDown menu, string label, VariableSupportMode mode)
    {
        Menu_AppendItem(menu, label, (_, _) => SetMode(mode), true, _mode == mode);
    }

    private void SetMode(VariableSupportMode mode)
    {
        if (_mode == mode)
            return;

        RecordUndoEvent("Set Variable Support Mode");
        _mode = mode;
        UpdateParameterVisibility();
        UpdateMessage();
        ExpireSolution(true);
    }

    private void UpdateMessage()
    {
        Message = _mode switch
        {
            VariableSupportMode.Sphere => "VS: Sphere",
            VariableSupportMode.Roller => "VS: Roller",
            VariableSupportMode.Rail => "VS: Rail",
            VariableSupportMode.NurbsCurve => "VS: Curve",
            VariableSupportMode.NurbsSurface => "VS: Surface",
            _ => "VS"
        };
    }

    private void UpdateParameterVisibility()
    {
        if (Params.Input.Count == 0)
            return;

        for (int i = Params.Input.Count - 1; i >= 1; i--)
            Params.UnregisterInputParameter(Params.Input[i], true);

        foreach (var param in CreateModeParams())
            Params.RegisterInputParam(param);

        Params.OnParametersChanged();
    }

    private IEnumerable<IGH_Param> CreateModeParams()
    {
        switch (_mode)
        {
            case VariableSupportMode.Sphere:
                yield return NumberParam("Radius", "R", "Relative radius from initial node position", true);
                break;
            case VariableSupportMode.Roller:
                yield return BooleanParam("Free X", "Fx", "Allow movement in X", true);
                yield return BooleanParam("Free Y", "Fy", "Allow movement in Y", true);
                yield return BooleanParam("Free Z", "Fz", "Allow movement in Z", true);
                yield return IntervalParam("X Domain", "Dx", "Relative X bounds from initial position", true);
                yield return IntervalParam("Y Domain", "Dy", "Relative Y bounds from initial position", true);
                yield return IntervalParam("Z Domain", "Dz", "Relative Z bounds from initial position", true);
                break;
            case VariableSupportMode.Rail:
                yield return LineParam("Rail", "Rail", "Line segment rail constraint", false);
                break;
            case VariableSupportMode.NurbsCurve:
                yield return CurveParam("Curve", "C", "NURBS curve constraint", false);
                yield return IntervalParam("Curve Domain", "Dc", "Optional curve parameter domain", true);
                break;
            case VariableSupportMode.NurbsSurface:
                yield return SurfaceParam("Surface", "S", "NURBS surface constraint", false);
                yield return IntervalParam("U Domain", "Du", "Optional surface U parameter domain", true);
                yield return IntervalParam("V Domain", "Dv", "Optional surface V parameter domain", true);
                break;
        }
    }

    private static Param_Number NumberParam(string name, string nickname, string description, bool optional) => new()
    {
        Name = name,
        NickName = nickname,
        Description = description,
        Access = GH_ParamAccess.item,
        Optional = optional,
    };

    private static Param_Boolean BooleanParam(string name, string nickname, string description, bool optional) => new()
    {
        Name = name,
        NickName = nickname,
        Description = description,
        Access = GH_ParamAccess.item,
        Optional = optional,
    };

    private static Param_Interval IntervalParam(string name, string nickname, string description, bool optional) => new()
    {
        Name = name,
        NickName = nickname,
        Description = description,
        Access = GH_ParamAccess.item,
        Optional = optional,
    };

    private static Param_Line LineParam(string name, string nickname, string description, bool optional) => new()
    {
        Name = name,
        NickName = nickname,
        Description = description,
        Access = GH_ParamAccess.item,
        Optional = optional,
    };

    private static Param_Curve CurveParam(string name, string nickname, string description, bool optional) => new()
    {
        Name = name,
        NickName = nickname,
        Description = description,
        Access = GH_ParamAccess.item,
        Optional = optional,
    };

    private static Param_Surface SurfaceParam(string name, string nickname, string description, bool optional) => new()
    {
        Name = name,
        NickName = nickname,
        Description = description,
        Access = GH_ParamAccess.item,
        Optional = optional,
    };

    public override bool Write(GH_IWriter writer)
    {
        writer.SetInt32(ModeKey, (int)_mode);
        return base.Write(writer);
    }

    public override bool Read(GH_IReader reader)
    {
        if (reader.ItemExists(ModeKey))
        {
            int value = reader.GetInt32(ModeKey);
            if (Enum.IsDefined(typeof(VariableSupportMode), value))
                _mode = (VariableSupportMode)value;
        }
        UpdateParameterVisibility();
        UpdateMessage();
        return base.Read(reader);
    }

    protected override Bitmap Icon => Properties.Resources.Var_Anchor;
    public override Guid ComponentGuid => new("6D83C7B3-D66B-43CF-AF5E-D94873094B04");
}
