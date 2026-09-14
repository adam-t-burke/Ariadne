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

    /// <summary>Target kinds offered in the auto-dropped Target value list.</summary>
    private static readonly string[] TargetChoices =
    [
        InverseFdmCaseNetworkBuilder.ExactTarget,
        "jit2pctd",
        "bump10pctd",
        InverseFdmCaseNetworkBuilder.DesignerTarget,
    ];

    public InverseFdmCaseComponent()
        : base("Load Inverse FDM Case", "InvFdmCase",
            "Load a warm-start benchmark case (synthetic suite, jax_fdm, compas_cem) as an FDM Network with free-node Target, Loads, Load Nodes, and per-edge Signs / Lower / Upper for Inverse FDM and Theseus Solve.",
            "Ariadne", "Experimental")
    {
    }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddTextParameter("Case", "Case",
            "Case name (filename stem). A dropdown of the bundled cases is added when the component is placed.",
            GH_ParamAccess.item);
        pManager.AddTextParameter("Target", "Target",
            "Which target to load: 'exact' (the case's reference equilibrium), 'jit2pctd' / 'bump10pctd' (synthetic suite perturbations), or 'designer' (the upstream example's own goal). Unavailable kinds fall back to exact with a warning; see the Targets output.",
            GH_ParamAccess.item, InverseFdmCaseNetworkBuilder.ExactTarget);
        pManager.AddBooleanParameter("Snug Box", "Snug",
            "Synthetic suite: use the snug ×/÷1 box around q_ref instead of the loose ×/÷100 box.",
            GH_ParamAccess.item, false);
        pManager.AddTextParameter("Folder", "Folder",
            "Optional folder of case JSON files (e.g. bench/synthetic/cases after re-export). Empty uses the cases bundled in the plugin.",
            GH_ParamAccess.item);
        pManager[3].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Network", "Network", "FDM Network at the chosen target geometry (topology + anchors). Wire into Deconstruct Network, Inverse FDM, or Theseus Solve.", GH_ParamAccess.item);
        pManager.AddCurveParameter("Edges", "Edges", "Network edges as line curves (preview / inspection).", GH_ParamAccess.list);
        pManager.AddPointParameter("Anchors", "Anchors", "Fixed (support) node positions.", GH_ParamAccess.list);
        pManager.AddPointParameter("Target Points", "Target", "Target position of every free node, in free-node order (Inverse FDM Target).", GH_ParamAccess.list);
        pManager.AddVectorParameter("Loads", "Loads", "Load on every free node, in free-node order.", GH_ParamAccess.list);
        pManager.AddPointParameter("Load Nodes", "LN", "Free-node positions parallel to Loads (Inverse FDM / Theseus Load Nodes).", GH_ParamAccess.list);
        pManager.AddIntegerParameter("Signs", "Signs", "+1 tension, -1 compression, 0 free; one per edge on path {0}.", GH_ParamAccess.tree);
        pManager.AddNumberParameter("Lower", "Lower", "Lower bound on q per edge on path {0}; empty when the case has no lower bounds.", GH_ParamAccess.tree);
        pManager.AddNumberParameter("Upper", "Upper", "Upper bound on q per edge on path {0}; empty when the case has no upper bounds.", GH_ParamAccess.tree);
        pManager.AddNumberParameter("Reference Q", "Qref", "Force densities that reproduce the exact geometry. Not written to the Network, so Inverse FDM does not start at the solution.", GH_ParamAccess.list);
        pManager.AddTextParameter("Name", "Name", "Loaded case id and target kind.", GH_ParamAccess.item);
        pManager.AddTextParameter("Targets", "Targets", "Target kinds this case can serve.", GH_ParamAccess.list);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        string caseName = "";
        string targetKind = InverseFdmCaseNetworkBuilder.ExactTarget;
        bool snug = false;
        string folder = "";
        if (!DA.GetData(0, ref caseName)) return;
        DA.GetData(1, ref targetKind);
        DA.GetData(2, ref snug);
        DA.GetData(3, ref folder);

        if (string.IsNullOrWhiteSpace(caseName))
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Pick a case name.");
            return;
        }

        InverseFdmCaseNetwork built;
        try
        {
            var c = InverseFdmCaseLibrary.Load(caseName, folder);
            built = InverseFdmCaseNetworkBuilder.Build(c, targetKind, snug);
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

        var edges = new List<GH_Curve>(built.Network.Graph.Ne);
        foreach (var edge in built.Network.Graph.Edges)
            edges.Add(new GH_Curve(edge.Value));

        DA.SetData(0, built.Network);
        DA.SetDataList(1, edges);
        DA.SetDataList(2, built.Network.Anchors);
        DA.SetDataList(3, built.Target);
        DA.SetDataList(4, built.Loads);
        DA.SetDataList(5, built.LoadNodes);
        DA.SetDataTree(6, signs);
        DA.SetDataTree(7, BoundsTree(built.Lower, double.NegativeInfinity));
        DA.SetDataTree(8, BoundsTree(built.Upper, double.PositiveInfinity));
        DA.SetDataList(9, built.QRef);
        DA.SetData(10, $"{built.Name} · {built.TargetKind}");
        DA.SetDataList(11, built.AvailableTargets);

        Message = $"{built.Network.Free.Count} free · {built.Network.Graph.Ne} edges · {built.TargetKind}"
            + (built.UsedSnugBox ? " · snug" : "");
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
    /// Drops value lists for the bundled case names (Case) and the target kinds (Target) the first
    /// time the component is placed. Components read back from a file (or pasted) skip this so saved
    /// wiring is untouched.
    /// </summary>
    public override void AddedToDocument(GH_Document document)
    {
        base.AddedToDocument(document);
        if (_dropdownAdded) return;
        _dropdownAdded = true;

        var pivot = Attributes.Pivot;
        if (Params.Input[0].SourceCount == 0 && InverseFdmCaseLibrary.BundledNames.Count > 0)
            AddDropdown(document, 0, "Case", InverseFdmCaseLibrary.BundledNames, new PointF(pivot.X - 280, pivot.Y - 12));
        if (Params.Input[1].SourceCount == 0)
            AddDropdown(document, 1, "Target", TargetChoices, new PointF(pivot.X - 280, pivot.Y + 20));
    }

    private void AddDropdown(GH_Document document, int input, string nickname, IReadOnlyList<string> items, PointF pivot)
    {
        var list = new GH_ValueList();
        list.CreateAttributes();
        list.NickName = nickname;
        list.ListMode = GH_ValueListMode.DropDown;
        list.ListItems.Clear();
        foreach (string item in items)
            list.ListItems.Add(new GH_ValueListItem(item, $"\"{item}\""));
        list.SelectItem(0);
        list.Attributes.Pivot = pivot;
        document.AddObject(list, false);
        Params.Input[input].AddSource(list);
    }

    public override bool Read(GH_IReader reader)
    {
        // Anything deserialised was placed before; never spawn a second dropdown.
        _dropdownAdded = true;
        return base.Read(reader);
    }

    protected override string HtmlHelp_Source() =>
        "<h1>Load Inverse FDM Case</h1>" +
        "<p>Loads one warm-start benchmark structure as an FDM Network: the 16 synthetic suite nets from the write-up " +
        "(<code>bench/synthetic/cases</code>, exported by <code>warm_start_bench export_cases</code>) and the " +
        "<b>jax_fdm</b> / <b>compas_cem</b> examples (<code>bench/external/cases</code>).</p>" +
        "<ul>" +
        "<li><b>Network</b> sits at the chosen target geometry, so Theseus <b>Target XYZ</b> fits the same shape Inverse FDM receives. " +
        "<b>Edges</b> and <b>Anchors</b> are the same geometry for preview.</li>" +
        "<li><b>Target</b>, <b>Loads</b>, and <b>Load Nodes</b> are per free node in free-node order and plug straight into " +
        "Inverse FDM / Theseus Solve.</li>" +
        "<li><b>Signs</b>, <b>Lower</b>, <b>Upper</b> are per edge on path {0}. Unbounded sides are left empty; a side with only " +
        "some unbounded entries uses ±∞ for those entries. <b>Snug Box</b> switches the suite nets to the ×/÷1 box.</li>" +
        "<li><b>Qref</b> is the force density that reproduces the <i>exact</i> geometry. It is not stored on the Network edges, " +
        "so Inverse FDM is not seeded at the answer. Wire it into Theseus <b>q</b> for a verification forward solve.</li>" +
        "</ul>" +
        "<p><b>Target</b>: <code>exact</code> is the reference equilibrium; <code>jit2pctd</code> (2 % white noise on the depth) and " +
        "<code>bump10pctd</code> (smooth 10 % push) are the suite's perturbed targets; <code>designer</code> is the upstream example's own goal " +
        "(<code>target_original</code>, e.g. <code>jaxfdm_creased_shell</code>). Kinds a case lacks fall back to exact with a warning; " +
        "the <b>Targets</b> output lists what is available.</p>" +
        "<p>Sign convention: q &gt; 0 tension. Node indices are 0-based and preserved from the JSON, with free nodes first.</p>";

    protected override Bitmap Icon => Properties.Resources.parameters;
    public override Guid ComponentGuid => new("7C1E9B2A-4F3D-4E58-9A6B-2D8F0C5E1A47");
}
