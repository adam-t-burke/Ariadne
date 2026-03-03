namespace Ariadne.Solver;

using System;
using System.Collections.Generic;
using Ariadne.FDM;
using Ariadne.Graphs;
using Rhino.Geometry;

internal enum EdgeKeyKind
{
    ReferenceGuid,
    UserId,
    Geometry,
    TreeLocation,
}

internal readonly record struct EdgeKey(
    EdgeKeyKind Kind,
    Guid ReferenceId,
    string UserId,
    int Branch,
    int Item,
    long Ax,
    long Ay,
    long Az,
    long Bx,
    long By,
    long Bz,
    long Length);

internal sealed record EdgeKeySet(EdgeKey[] Keys, HashSet<EdgeKey> AmbiguousKeys);

internal static class EdgeIdentity
{
    public static EdgeKeySet BuildKeys(FDM_Network network)
    {
        var graph = network.Graph;
        var keys = new EdgeKey[graph.Ne];
        var counts = new Dictionary<EdgeKey, int>();

        for (int i = 0; i < graph.Ne; i++)
        {
            keys[i] = BuildKey(network, i);
            counts.TryGetValue(keys[i], out int count);
            counts[keys[i]] = count + 1;
        }

        var ambiguous = new HashSet<EdgeKey>();
        foreach (var (key, count) in counts)
        {
            if (count > 1)
                ambiguous.Add(key);
        }

        return new EdgeKeySet(keys, ambiguous);
    }

    private static EdgeKey BuildKey(FDM_Network network, int edgeIndex)
    {
        var graph = network.Graph;
        var edge = graph.Edges[edgeIndex];

        if (edge.ReferenceID != Guid.Empty)
        {
            return new EdgeKey(
                EdgeKeyKind.ReferenceGuid,
                edge.ReferenceID,
                "",
                -1,
                -1,
                0,
                0,
                0,
                0,
                0,
                0,
                0);
        }

        if (!string.IsNullOrWhiteSpace(edge.UserKey))
        {
            return new EdgeKey(
                EdgeKeyKind.UserId,
                Guid.Empty,
                edge.UserKey,
                -1,
                -1,
                0,
                0,
                0,
                0,
                0,
                0,
                0);
        }

        if (TryBuildGeometryKey(network, edge, out var geometryKey))
            return geometryKey;

        if (graph.EdgeInputMap is { Count: var count } && count == graph.Ne)
        {
            var (branch, item) = graph.EdgeInputMap[edgeIndex];
            return new EdgeKey(
                EdgeKeyKind.TreeLocation,
                Guid.Empty,
                "",
                branch,
                item,
                0,
                0,
                0,
                0,
                0,
                0,
                0);
        }

        return new EdgeKey(
            EdgeKeyKind.TreeLocation,
            Guid.Empty,
            "",
            0,
            edgeIndex,
            0,
            0,
            0,
            0,
            0,
            0,
            0);
    }

    private static bool TryBuildGeometryKey(FDM_Network network, Edge edge, out EdgeKey key)
    {
        double tol = network.ETol > 0 ? network.ETol : network.Graph.Tolerance;
        if (!double.IsFinite(tol) || tol <= 0)
            tol = 1e-6;

        var a = Quantize(edge.Start.Value, tol);
        var b = Quantize(edge.End.Value, tol);
        if (Compare(a, b) > 0)
            (a, b) = (b, a);

        long length = (long)Math.Round(edge.Value.GetLength() / tol);
        key = new EdgeKey(
            EdgeKeyKind.Geometry,
            Guid.Empty,
            "",
            -1,
            -1,
            a.X,
            a.Y,
            a.Z,
            b.X,
            b.Y,
            b.Z,
            length);
        return true;
    }

    private static QuantizedPoint Quantize(Point3d point, double tol) =>
        new(
            (long)Math.Round(point.X / tol),
            (long)Math.Round(point.Y / tol),
            (long)Math.Round(point.Z / tol));

    private static int Compare(QuantizedPoint left, QuantizedPoint right)
    {
        int x = left.X.CompareTo(right.X);
        if (x != 0) return x;
        int y = left.Y.CompareTo(right.Y);
        if (y != 0) return y;
        return left.Z.CompareTo(right.Z);
    }

    private readonly record struct QuantizedPoint(long X, long Y, long Z);
}
