using System;
using System.Collections.Generic;
using Ariadne.FDM;
using Ariadne.Graphs;
using Ariadne.Solver;
using Rhino.Geometry;
using Xunit;

namespace Ariadne.Tests;

public class EdgeIdentityTests
{
    [Fact]
    public void BuildKeys_PrefersReferenceId()
    {
        var referenceId = Guid.NewGuid();
        var network = NetworkWithSingleEdge(referenceId);

        var keys = EdgeIdentity.BuildKeys(network);

        Assert.Equal(EdgeKeyKind.ReferenceGuid, keys.Keys[0].Kind);
        Assert.Equal(referenceId, keys.Keys[0].ReferenceId);
        Assert.Empty(keys.AmbiguousKeys);
    }

    [Fact]
    public void FdmNetworkCopy_PreservesReferenceIds()
    {
        var referenceId = Guid.NewGuid();
        var network = NetworkWithSingleEdge(referenceId);

        var copy = new FDM_Network(network);

        Assert.NotSame(network.Graph.Edges[0], copy.Graph.Edges[0]);
        Assert.Equal(referenceId, copy.Graph.Edges[0].ReferenceID);
    }

    [Fact]
    public void BuildKeys_UsesGeometryForProceduralEdges()
    {
        var network = NetworkWithSingleEdge(Guid.Empty);

        var keys = EdgeIdentity.BuildKeys(network);

        Assert.Equal(EdgeKeyKind.Geometry, keys.Keys[0].Kind);
        Assert.Empty(keys.AmbiguousKeys);
    }

    private static FDM_Network NetworkWithSingleEdge(Guid referenceId)
    {
        var start = new Node { Value = new Point3d(0, 0, 0), Index = 0 };
        var end = new Node { Value = new Point3d(1, 0, 0), Index = 1 };
        start.Neighbors.Add(end);
        end.Neighbors.Add(start);

        var edge = new Edge
        {
            Start = start,
            End = end,
            Value = new LineCurve(start.Value, end.Value),
            ReferenceID = referenceId,
        };

        return new FDM_Network
        {
            Graph = new Graph
            {
                Tolerance = 0.001,
                Nodes = new List<Node> { start, end },
                Edges = new List<Edge> { edge },
                EdgeInputMap = new List<(int, int)> { (0, 0) },
            },
            ETol = 0.001,
            ATol = 0.001,
            Anchors = [],
            Free = [start, end],
            Fixed = [],
            FreeNodes = [0, 1],
            FixedNodes = [],
            Valid = true,
        };
    }
}
