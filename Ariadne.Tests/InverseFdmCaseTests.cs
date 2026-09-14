namespace Ariadne.Tests;

using System;
using System.IO;
using System.Linq;
using Ariadne.Solver;
using Ariadne.Solver.Components.Experimental;
using Grasshopper.Kernel.Data;
using Xunit;

public sealed class InverseFdmCaseTests
{
    [Fact]
    public void BundlesAllExternalCases()
    {
        var names = InverseFdmCaseLibrary.BundledNames;

        Assert.Equal(18, names.Count);
        Assert.Contains("jaxfdm_arch_loadpath", names);
        Assert.Contains("cem_bridge_2d", names);
        Assert.Equal(names.OrderBy(n => n, StringComparer.Ordinal), names);
    }

    [Fact]
    public void ArchBuildsValidFreeFirstNetwork()
    {
        var built = InverseFdmCaseNetworkBuilder.Build(
            InverseFdmCaseLibrary.LoadBundled("jaxfdm_arch_loadpath"), useDesignerTarget: false);
        var network = built.Network;

        Assert.True(network.Valid);
        Assert.Equal(11, network.Graph.Nn);
        Assert.Equal(10, network.Graph.Ne);
        Assert.Equal(9, network.Free.Count);
        Assert.Equal([9, 10], network.FixedNodes);
        Assert.Equal(Enumerable.Range(0, 9), network.FreeNodes);
        Assert.All(network.Fixed, n => Assert.True(n.Anchor));
        Assert.All(network.Free, n => Assert.False(n.Anchor));

        // JSON edge (9,10): node 9 is the last free node (index 8), node 10 the second anchor (index 10).
        var last = network.Graph.Edges[9];
        Assert.Equal(8, last.Start.Index);
        Assert.Equal(10, last.End.Index);
        Assert.Equal(0, network.Graph.Edges[0].Start.Index);
        Assert.Equal(1, network.Graph.Edges[0].End.Index);

        Assert.Equal(9, built.Target.Count);
        Assert.Equal(9, built.Loads.Count);
        Assert.Equal(9, built.LoadNodes.Count);
        Assert.Equal(built.Target, built.LoadNodes);
        Assert.All(built.Loads, v => Assert.Equal(-0.3, v.Z, 12));

        Assert.Equal(10, built.QRef.Length);
        Assert.Equal(10, built.Signs.Length);
        Assert.All(built.Signs, s => Assert.Equal(-1, s));
        Assert.All(built.Lower, v => Assert.True(v.HasValue && v < 0));
        Assert.All(built.Upper, v => Assert.True(v.HasValue && v < 0));
        Assert.False(built.UsedDesignerTarget);
        Assert.Empty(built.Warnings);
    }

    [Fact]
    public void ReferenceQIsNotSeededOnEdges()
    {
        var built = InverseFdmCaseNetworkBuilder.Build(
            InverseFdmCaseLibrary.LoadBundled("jaxfdm_arch_loadpath"), useDesignerTarget: false);

        Assert.All(built.Network.Graph.Edges, e => Assert.Equal(0.0, e.Q));
        Assert.Contains(built.QRef, q => q != 0.0);
    }

    [Fact]
    public void EdgePathsMapFlatTreesOntoEveryEdge()
    {
        var built = InverseFdmCaseNetworkBuilder.Build(
            InverseFdmCaseLibrary.LoadBundled("jaxfdm_truss_equal_force"), useDesignerTarget: false);
        var paths = built.Network.Graph.EdgeInputPaths;

        Assert.Equal(built.Network.Graph.Ne, paths.Count);
        Assert.All(paths, p => Assert.Equal(new GH_Path(0), p));

        var signs = new EdgeValueTree([new EdgeValueBranch("{0}", built.Signs.Select(s => (double)s).ToArray())]);
        var mapped = QTreeMapper.Map(signs, paths, "Signs");
        Assert.True(mapped.Success);
        Assert.Null(mapped.Warning);
        Assert.Contains(mapped.Values!, s => s > 0);
        Assert.Contains(mapped.Values!, s => s < 0);
    }

    [Fact]
    public void LoadNodesResolveOneToOneOntoFreeNodes()
    {
        var built = InverseFdmCaseNetworkBuilder.Build(
            InverseFdmCaseLibrary.LoadBundled("cem_tree_canopy_3d"), useDesignerTarget: false);

        var indices = TheseusSolverService.ResolveLoadNodeIndices(built.Network, built.LoadNodes);
        Assert.Equal(Enumerable.Range(0, built.Network.Free.Count), indices);

        var packed = TheseusSolverService.PackFreeNodeLoads(built.Network.Free.Count, built.Loads, indices);
        Assert.Equal(built.Loads, packed);
    }

    [Fact]
    public void NullBoundSidesStayEmpty()
    {
        var cablenet = InverseFdmCaseNetworkBuilder.Build(
            InverseFdmCaseLibrary.LoadBundled("jaxfdm_cablenet"), useDesignerTarget: false);
        Assert.Null(InverseFdmCaseNetwork.DenseBounds(cablenet.Upper, double.PositiveInfinity));
        var lower = InverseFdmCaseNetwork.DenseBounds(cablenet.Lower, double.NegativeInfinity);
        Assert.NotNull(lower);
        Assert.All(lower!, v => Assert.True(double.IsFinite(v) && v > 0));

        var cem = InverseFdmCaseNetworkBuilder.Build(
            InverseFdmCaseLibrary.LoadBundled("cem_bridge_2d"), useDesignerTarget: false);
        Assert.Null(InverseFdmCaseNetwork.DenseBounds(cem.Lower, double.NegativeInfinity));
        Assert.Null(InverseFdmCaseNetwork.DenseBounds(cem.Upper, double.PositiveInfinity));

        double?[] mixed = [1.0, null, 2.0];
        var dense = InverseFdmCaseNetwork.DenseBounds(mixed, double.NegativeInfinity);
        Assert.NotNull(dense);
        Assert.Equal([1.0, double.NegativeInfinity, 2.0], dense!);
    }

    [Fact]
    public void PartialDesignerTargetFallsBackWithWarning()
    {
        var built = InverseFdmCaseNetworkBuilder.Build(
            InverseFdmCaseLibrary.LoadBundled("cem_bridge_2d"), useDesignerTarget: true);

        Assert.False(built.UsedDesignerTarget);
        Assert.Single(built.Warnings);
        Assert.Contains("designer target", built.Warnings[0]);
    }

    [Fact]
    public void CompleteDesignerTargetMovesFreeNodesAndTarget()
    {
        var exact = InverseFdmCaseNetworkBuilder.Build(
            InverseFdmCaseLibrary.LoadBundled("jaxfdm_creased_shell"), useDesignerTarget: false);
        var designer = InverseFdmCaseNetworkBuilder.Build(
            InverseFdmCaseLibrary.LoadBundled("jaxfdm_creased_shell"), useDesignerTarget: true);

        Assert.True(designer.UsedDesignerTarget);
        Assert.Empty(designer.Warnings);
        Assert.NotEqual(exact.Target, designer.Target);
        Assert.Equal(designer.Target, designer.Network.Free.Select(n => n.Value).ToList());
        Assert.Equal(exact.Network.Anchors, designer.Network.Anchors);
    }

    [Fact]
    public void ParseRejectsInconsistentLengths()
    {
        const string json = """
            {"name":"bad","nodes":[[0,0,0],[1,0,0],[2,0,0]],"edges":[[0,1],[1,2]],"fixed":[0,2],
             "loads":[[0,0,0],[0,0,-1],[0,0,0]],"q_ref":[1.0],"target":[[0,0,0],[1,0,0],[2,0,0]]}
            """;

        var ex = Assert.Throws<InvalidOperationException>(() => InverseFdmCaseLibrary.Parse(json));
        Assert.Contains("q_ref", ex.Message);
    }

    [Fact]
    public void LoadsFromFolderOverride()
    {
        string folder = Path.Combine(Path.GetTempPath(), "ariadne-cases-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(folder);
        try
        {
            const string json = """
                {"name":"tiny","nodes":[[0,0,0],[1,0,1],[2,0,0]],"edges":[[0,1],[1,2]],"fixed":[0,2],
                 "loads":[[0,0,0],[0,0,-1],[0,0,0]],"q_ref":[-1.0,-1.0],"target":[[0,0,0],[1,0,1],[2,0,0]]}
                """;
            File.WriteAllText(Path.Combine(folder, "tiny.json"), json);

            var c = InverseFdmCaseLibrary.Load("tiny", folder);
            var built = InverseFdmCaseNetworkBuilder.Build(c, useDesignerTarget: false);

            Assert.Equal("tiny", built.Name);
            Assert.True(built.Network.Valid);
            Assert.Single(built.Target);
            Assert.Equal([-1, -1], built.Signs);
            Assert.All(built.Lower, v => Assert.Null(v));
            Assert.Throws<FileNotFoundException>(() => InverseFdmCaseLibrary.Load("missing", folder));
        }
        finally
        {
            Directory.Delete(folder, recursive: true);
        }
    }
}
