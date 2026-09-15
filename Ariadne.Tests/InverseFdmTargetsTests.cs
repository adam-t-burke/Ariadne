using System;
using System.Collections.Generic;
using System.Linq;
using Ariadne.Solver.Components.Experimental;
using Rhino.Geometry;
using Xunit;

namespace Ariadne.Tests;

public sealed class InverseFdmTargetsTests
{
    private static readonly Point3d[] Free =
    [
        new(0, 0, 0),
        new(1, 0, 0),
        new(2, 0, 0),
        new(3, 0, 0),
    ];

    [Fact]
    public void FullListWithoutNodesIsPackedInOrderWithoutWeights()
    {
        var targets = Free.Select(p => new Point3d(p.X, p.Y, 1.0)).ToList();

        var set = InverseFdmTargets.Complete(Free, targets, null, 0.5);

        Assert.False(set.IsPartial);
        Assert.Null(set.Weights);
        Assert.Equal([0, 0, 1, 1, 0, 1, 2, 0, 1, 3, 0, 1], set.FreeXyz);
        Assert.Equal(Enumerable.Range(0, 4), set.TargetedFree);
        Assert.Equal(targets, set.CompletedPoints);
    }

    [Fact]
    public void FullListWithWrongCountIsRejected()
    {
        var ex = Assert.Throws<ArgumentException>(
            () => InverseFdmTargets.Complete(Free, [new Point3d(0, 0, 1)], null, 0.01));
        Assert.Contains("must match free node count (4)", ex.Message);
    }

    [Fact]
    public void SubsetCompletesTheOthersFromCurrentPositionsWithHoldWeight()
    {
        var targets = new List<Point3d> { new(3, 0, 2), new(1, 0, 1) };

        var set = InverseFdmTargets.Complete(Free, targets, [3, 1], 0.01);

        Assert.True(set.IsPartial);
        Assert.Equal([0.01, 1.0, 0.01, 1.0], set.Weights);
        Assert.Equal([1, 3], set.TargetedFree);
        Assert.Equal(Free[0], set.CompletedPoints[0]);
        Assert.Equal(new Point3d(1, 0, 1), set.CompletedPoints[1]);
        Assert.Equal(Free[2], set.CompletedPoints[2]);
        Assert.Equal(new Point3d(3, 0, 2), set.CompletedPoints[3]);
        Assert.Equal([0, 0, 0, 1, 0, 1, 2, 0, 0, 3, 0, 2], set.FreeXyz);
    }

    [Fact]
    public void SubsetCountMismatchIsRejected()
    {
        var ex = Assert.Throws<ArgumentException>(
            () => InverseFdmTargets.Complete(Free, [new Point3d(1, 0, 1)], [1, 2], 0.01));
        Assert.Contains("must match Target Nodes count (2)", ex.Message);
    }

    [Fact]
    public void DuplicateAndOutOfRangeNodesAreRejected()
    {
        var two = new List<Point3d> { new(1, 0, 1), new(1, 0, 2) };
        var duplicate = Assert.Throws<ArgumentException>(
            () => InverseFdmTargets.Complete(Free, two, [1, 1], 0.01));
        Assert.Contains("more than once", duplicate.Message);

        var outside = Assert.Throws<ArgumentException>(
            () => InverseFdmTargets.Complete(Free, two, [1, 4], 0.01));
        Assert.Contains("outside the free-node range", outside.Message);
    }

    [Theory]
    [InlineData(0.0)]
    [InlineData(-0.1)]
    [InlineData(1.5)]
    [InlineData(double.NaN)]
    [InlineData(double.PositiveInfinity)]
    public void HoldWeightOutsideUnitIntervalIsRejected(double holdWeight)
    {
        var ex = Assert.Throws<ArgumentException>(
            () => InverseFdmTargets.Complete(Free, [new Point3d(1, 0, 1)], [1], holdWeight));
        Assert.Contains("Hold Weight must be in (0, 1]", ex.Message);
    }

    [Fact]
    public void HoldWeightIsIgnoredForAFullTarget()
    {
        var targets = Free.ToList();
        var set = InverseFdmTargets.Complete(Free, targets, [], 0.0);
        Assert.False(set.IsPartial);
    }

    [Fact]
    public void TargetErrorCoversOnlyTheTargetedNodes()
    {
        var set = InverseFdmTargets.Complete(Free, [new Point3d(1, 0, 1)], [1], 0.01);
        // Graph nodes: free node i sits at graph index i + 1 (index 0 is an anchor).
        int[] freeNodes = [1, 2, 3, 4];
        var solved = new List<Point3d>
        {
            new(-1, 0, 0),   // anchor
            new(0, 0, 5),    // free 0, untargeted: large miss must not count
            new(1, 0, 1.5),  // free 1, targeted: 0.5 off
            new(2, 0, 5),
            new(3, 0, 5),
        };

        Assert.Equal(0.5, InverseFdmTargets.TargetError(solved, freeNodes, set), 12);
    }
}
