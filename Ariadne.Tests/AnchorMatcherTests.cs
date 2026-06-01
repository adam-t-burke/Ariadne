using System.Collections.Generic;
using Ariadne.Graphs;
using Ariadne.Utilities;
using Rhino.Geometry;
using Xunit;

namespace Ariadne.Tests;

public class AnchorMatcherTests
{
    [Fact]
    public void FindNodeIndex_MatchesWithinTolerance()
    {
        var nodes = new List<Node>
        {
            new() { Value = new Point3d(0, 0, 0) },
            new() { Value = new Point3d(1, 0, 0) },
            new() { Value = new Point3d(0, 1, 0) },
        };

        int idx = AnchorMatcher.FindNodeIndex(nodes, new Point3d(1.0005, 0, 0), 0.01);
        Assert.Equal(1, idx);
    }

    [Fact]
    public void FindNodeIndex_ReturnsNegativeWhenNoMatch()
    {
        var nodes = new List<Node> { new() { Value = new Point3d(0, 0, 0) } };
        int idx = AnchorMatcher.FindNodeIndex(nodes, new Point3d(5, 0, 0), 0.001);
        Assert.Equal(-1, idx);
    }
}
