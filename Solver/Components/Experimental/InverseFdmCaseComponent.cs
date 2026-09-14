using System;
using System.Collections.Generic;
using System.Drawing;
using GH_IO.Serialization;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Special;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

namespace Ariadne.Solver.Components.Experimental;

/// <summary>
/// Loads one external inverse-FDM benchmark case (jax_fdm / compas_cem export) as an FDM Network
/// plus the free-node Target, Loads, Load Nodes, and per-edge Signs / bounds that Inverse FDM and
/// Theseus Solve take directly.
/// </summary>
public class InverseFdmCaseComponent : GH_Component
{
    private bool _dropdownAdded;

    public InverseFdmCaseComponent()
        : base("Load Inverse FDM Case", "InvFdmCase",
            "Load a jax_fdm / compas_cem benchmark case as an FDM Network with free-node Target, Loads, Load Nodes, and per-edge Signs / Lower / Upper for Inverse FDM and Theseus Solve.",
            "Ariadne", "Experimental")
    {
    }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddTextParameter("Case", "Case",
            "Case name (filename stem). A dropdown of the bundled cases is added when the component is placed.",
            GH_ParamAccess.item);
        pManager.AddBooleanParameter("Designer Target", "Designer",
            "Use the example's own design target (target_original) when every free node has one. Falls back to the exact target with a warning otherwise.",
            GH_ParamAccess.item, false);
        pManager.AddTextParameter("Folder", "Folder",
            "Optional folder of case JSON files (e.g. bench/external/cases after re-export). Empty uses the cases bundled in the plugin.",
            GH_ParamAccess.item);
        pManager[2].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Network", "Network", "FDM Network at the case geometry (topology + anchors). Wire into Deconstruct Network, Inverse FDM, or Theseus Solve.", GH_ParamAccess.item);
        pManager.AddPointParameter("Target Points", "Target", "Target position of every free node, in free-node order (Inverse FDM Target).", GH_ParamAccess.list);
        pManager.AddVectorParameter("Loads", "Loads", "Load on every free node, in free-node order.", GH_ParamAccess.list);
        pManager.AddPointParameter("Load Nodes", "LN", "Free-node positions parallel to Loads (Inverse FDM / Theseus Load Nodes).", GH_ParamAccess.list);
        pManager.AddIntegerParameter("Signs", "Signs", "+1 tension, -1 compression, 0 free; one per edge on path {0}.", GH_ParamAccess.tree);
        pManager.AddNumberParameter("Lower", "Lower", "Lower bound on q per edge on path {0}; empty when the case has no lower bounds.", GH_ParamAccess.tree);
        pManager.AddNumberParameter("Upper", "Upper", "Upper bound on q per edge on path {0}; empty when the case has no upper bounds.", GH_ParamAccess.tree);
        pManager.AddNumberParameter("Reference Q", "Qref", "Force densities that reproduce the case geometry. Not written to the Network, so Inverse FDM does not start at the solution.", GH_ParamAccess.list);
        pManager.AddTextParameter("Name", "Name", "Loaded case id.", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        string caseName = "";
        bool designer = false;
        string folder = "";
        if (!DA.GetData(0, ref caseName)) return;
        DA.GetData(1, ref designer);
        DA.GetData(2, ref folder);

        if (string.IsNullOrWhiteSpace(caseName))
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Pick a case name.");
            return;
        }

        InverseFdmCaseNetwork built;
        try
        {
            var c = InverseFdmCaseLibrary.Load(caseName, folder);
            built = InverseFdmCaseNetworkBuilder.Build(c, designer);
        }
        catch (Exception ex)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, ex.Message);
            Message = null;
            return;
        }

        foreach (string warning in built.Warnings)
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, warning);
        if (!built.Network.Valid)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Loaded network is not valid.");
            return;
        }

        var signs = new GH_Structure<GH_Integer>();
        var path = new GH_Path(0);
        foreach (int s in built.Signs)
            signs.Append(new GH_Integer(s), path);

        DA.SetData(0, built.Network);
        DA.SetDataList(1, built.Target);
        DA.SetDataList(2, built.Loads);
        DA.SetDataList(3, built.LoadNodes);
        DA.SetDataTree(4, signs);
        DA.SetDataTree(5, BoundsTree(built.Lower, double.NegativeInfinity));
        DA.SetDataTree(6, BoundsTree(built.Upper, double.PositiveInfinity));
        DA.SetDataList(7, built.QRef);
        DA.SetData(8, built.Name);

        Message = $"{built.Network.Free.Count} free · {built.Network.Graph.Ne} edges"
            + (built.UsedDesignerTarget ? " · designer" : "");
    }

    private static GH_Structure<GH_Number> BoundsTree(double?[] side, double fill)
    {
        var tree = new GH_Structure<GH_Number>();
        double[]? dense = InverseFdmCaseNetwork.DenseBounds(side, fill);
        if (dense is null) return tree;
        var path = new GH_Path(0);
        foreach (double v in dense)
            tree.Append(new GH_Number(v), path);
        return tree;
    }

    /// <summary>
    /// Drops a value list of the bundled case names onto the Case input the first time the component
    /// is placed. Components read back from a file (or pasted) skip this so saved wiring is untouched.
    /// </summary>
    public override void AddedToDocument(GH_Document document)
    {
        base.AddedToDocument(document);
        if (_dropdownAdded) return;
        _dropdownAdded = true;
        if (Params.Input[0].SourceCount > 0) return;

        var names = InverseFdmCaseLibrary.BundledNames;
        if (names.Count == 0) return;

        var list = new GH_ValueList();
        list.CreateAttributes();
        list.NickName = "Case";
        list.ListMode = GH_ValueListMode.DropDown;
        list.ListItems.Clear();
        foreach (string name in names)
            list.ListItems.Add(new GH_ValueListItem(name, $"\"{name}\""));
        list.SelectItem(0);

        var pivot = Attributes.Pivot;
        list.Attributes.Pivot = new PointF(pivot.X - 260, pivot.Y - 10);
        document.AddObject(list, false);
        Params.Input[0].AddSource(list);
    }

    public override bool Read(GH_IReader reader)
    {
        // Anything deserialised was placed before; never spawn a second dropdown.
        _dropdownAdded = true;
        return base.Read(reader);
    }

    protected override string HtmlHelp_Source() =>
        "<h1>Load Inverse FDM Case</h1>" +
        "<p>Loads one of the external benchmark structures exported from <b>jax_fdm</b> or <b>compas_cem</b> " +
        "(<code>bench/external/cases/*.json</code>) as an FDM Network at its reference equilibrium.</p>" +
        "<ul>" +
        "<li><b>Network</b> sits at the case geometry, so Theseus <b>Target XYZ</b> fits the same shape.</li>" +
        "<li><b>Target</b>, <b>Loads</b>, and <b>Load Nodes</b> are per free node in free-node order and plug straight into " +
        "Inverse FDM / Theseus Solve.</li>" +
        "<li><b>Signs</b>, <b>Lower</b>, <b>Upper</b> are per edge on path {0}. Unbounded sides are left empty; a side with only " +
        "some unbounded entries uses ±∞ for those entries.</li>" +
        "<li><b>Qref</b> is the known solution. It is not stored on the Network edges, so Inverse FDM is not seeded at the answer. " +
        "Wire it into Theseus <b>q</b> for a verification forward solve.</li>" +
        "</ul>" +
        "<p><b>Designer Target</b> swaps in the example's own goal geometry (<code>target_original</code>) where every free node has one " +
        "(e.g. <code>jaxfdm_creased_shell</code>); otherwise the exact target is used with a warning.</p>" +
        "<p>Sign convention: q &gt; 0 tension. Node indices are 0-based and preserved from the JSON, with free nodes first.</p>";

    protected override Bitmap Icon => Properties.Resources.parameters;
    public override Guid ComponentGuid => new("7C1E9B2A-4F3D-4E58-9A6B-2D8F0C5E1A47");
}
