using System;
using System.Collections.Generic;
using Ariadne.Graphs;
using Rhino.Geometry;

namespace Ariadne.Utilities;

/// <summary>
/// Spatial-hash anchor-to-node matching (same cell strategy as graph construction).
/// </summary>
internal static class AnchorMatcher
{
    public static int FindNodeIndex(IReadOnlyList<Node> nodes, Point3d anchor, double tolerance)
    {
        if (nodes.Count == 0)
            return -1;

        if (tolerance <= 0)
        {
            for (int i = 0; i < nodes.Count; i++)
            {
                if (nodes[i].Value.EpsilonEquals(anchor, 0))
                    return i;
            }
            return -1;
        }

        double tol2 = tolerance * tolerance;
        double cell = tolerance;
        var cells = BuildCellIndex(nodes, cell);

        int bestIdx = -1;
        double bestDist = double.MaxValue;
        var anchorKey = Key(anchor, cell);

        foreach (var nk in NeighborKeys(anchorKey))
        {
            if (!cells.TryGetValue(nk, out var list))
                continue;

            for (int j = 0; j < list.Count; j++)
            {
                int nodeIdx = list[j];
                double dist2 = nodes[nodeIdx].Value.DistanceToSquared(anchor);
                if (dist2 <= tol2 && dist2 < bestDist)
                {
                    bestDist = dist2;
                    bestIdx = nodeIdx;
                }
            }
        }

        return bestIdx;
    }

    private static Dictionary<(long, long, long), List<int>> BuildCellIndex(
        IReadOnlyList<Node> nodes, double cell)
    {
        var cells = new Dictionary<(long, long, long), List<int>>();
        for (int i = 0; i < nodes.Count; i++)
        {
            var key = Key(nodes[i].Value, cell);
            if (!cells.TryGetValue(key, out var bucket))
            {
                bucket = new List<int>(4);
                cells[key] = bucket;
            }
            bucket.Add(i);
        }
        return cells;
    }

    private static (long, long, long) Key(Point3d p, double c) =>
        ((long)Math.Floor(p.X / c),
         (long)Math.Floor(p.Y / c),
         (long)Math.Floor(p.Z / c));

    private static IEnumerable<(long, long, long)> NeighborKeys((long, long, long) k)
    {
        for (long dx = -1; dx <= 1; dx++)
            for (long dy = -1; dy <= 1; dy++)
                for (long dz = -1; dz <= 1; dz++)
                    yield return (k.Item1 + dx, k.Item2 + dy, k.Item3 + dz);
    }
}
