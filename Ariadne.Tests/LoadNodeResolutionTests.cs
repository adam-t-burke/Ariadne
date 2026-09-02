namespace Ariadne.Tests;

using System;
using System.Collections.Generic;
using Ariadne.FDM;
using Ariadne.Graphs;
using Ariadne.Solver;
using Rhino.Geometry;
using Xunit;

public sealed class LoadNodeResolutionTests
{
    [Fact]
    public void ResolvesPointsWithinNetworkTolerance()
    {
        var network = Network(
            tolerance: 0.01,
            new Point3d(0, 0, 0),
            new Point3d(1, 0, 0));

        var indices = TheseusSolverService.ResolveLoadNodeIndices(
            network,
            [new Point3d(1.005, 0, 0)]);

        Assert.Equal([1], indices);
    }

    [Fact]
    public void ResolvesEachPointToNearestFreeNode()
    {
        var network = Network(
            tolerance: 0.25,
            new Point3d(0, 0, 0),
            new Point3d(1, 0, 0),
            new Point3d(2, 0, 0));

        var indices = TheseusSolverService.ResolveLoadNodeIndices(
            network,
            [new Point3d(1.1, 0, 0), new Point3d(0.05, 0, 0)]);

        Assert.Equal([1, 0], indices);
    }

    [Fact]
    public void RejectsPointOutsideNetworkTolerance()
    {
        var network = Network(
            tolerance: 0.01,
            new Point3d(0, 0, 0),
            new Point3d(1, 0, 0));

        var exception = Assert.Throws<ArgumentException>(
            () => TheseusSolverService.ResolveLoadNodeIndices(
                network,
                [new Point3d(1.02, 0, 0)]));

        Assert.Contains("does not match any free node", exception.Message);
        Assert.Contains("0.01", exception.Message);
    }

    private static FDM_Network Network(double tolerance, params Point3d[] points)
    {
        var nodes = new List<Node>(points.Length);
        for (int i = 0; i < points.Length; i++)
            nodes.Add(new Node { Index = i, Value = points[i] });

        return new FDM_Network
        {
            Graph = new Graph
            {
                Nodes = nodes,
                Edges = [],
                EdgeInputMap = [],
                Tolerance = tolerance,
            },
            Anchors = [],
            ATol = tolerance,
            ETol = tolerance,
            Free = nodes,
            Fixed = [],
            FreeNodes = [.. nodes.ConvertAll(node => node.Index)],
            FixedNodes = [],
            Valid = true,
        };
    }
}
