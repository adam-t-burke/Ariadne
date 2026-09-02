namespace Ariadne.Tests;

using System;
using System.Collections.Generic;
using Ariadne.Solver;
using Xunit;

public sealed class QCachePolicyTests
{
    [Fact]
    public void FinitePartialOptimizationResultUpdatesCache()
    {
        var first = Key("first");
        var second = Key("second");
        var cache = new Dictionary<EdgeKey, double>();

        QCacheUpdater.Apply(cache, [first, second], [0.5, 1.25], new HashSet<EdgeKey>());

        Assert.Equal(0.5, cache[first]);
        Assert.Equal(1.25, cache[second]);
    }

    [Fact]
    public void NonFiniteResultPreservesPriorValueAndPrunesAbsentEdges()
    {
        var current = Key("current");
        var stale = Key("stale");
        var cache = new Dictionary<EdgeKey, double>
        {
            [current] = 2.0,
            [stale] = 3.0,
        };

        QCacheUpdater.Apply(cache, [current], [double.NaN], new HashSet<EdgeKey>());

        Assert.Equal(2.0, cache[current]);
        Assert.DoesNotContain(stale, cache);
    }

    [Fact]
    public void AmbiguousEdgeIsNeverCached()
    {
        var ambiguous = Key("ambiguous");
        var cache = new Dictionary<EdgeKey, double>();

        QCacheUpdater.Apply(cache, [ambiguous], [4.0], new HashSet<EdgeKey> { ambiguous });

        Assert.Empty(cache);
    }

    private static EdgeKey Key(string id) =>
        new(EdgeKeyKind.UserId, Guid.Empty, id, -1, -1, 0, 0, 0, 0, 0, 0, 0);
}
