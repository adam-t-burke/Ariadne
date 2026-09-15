using System;
using System.Collections.Generic;
using System.Linq;
using Rhino.Geometry;

namespace Ariadne.Solver.Components.Experimental;

/// <summary>
/// A complete free-node target for the inverse pipeline together with the
/// per-node weights that pose it on a subset of the free nodes.
/// </summary>
/// <param name="FreeXyz">Row-major XYZ for every free node, in free-node order.</param>
/// <param name="Weights">
/// Per-free-node geometric weights, or <c>null</c> when every free node carries
/// a designer target (all weights one).
/// </param>
/// <param name="CompletedPoints">The target per free node after completion.</param>
/// <param name="TargetedFree">Free-node indices that carry a designer target.</param>
internal sealed record InverseFdmTargetSet(
    double[] FreeXyz,
    double[]? Weights,
    IReadOnlyList<Point3d> CompletedPoints,
    IReadOnlyList<int> TargetedFree)
{
    /// <summary>True when only a subset of the free nodes carries a target.</summary>
    public bool IsPartial => Weights is not null;
}

/// <summary>
/// Builds the dense target the native pipeline expects from either a full
/// free-node list or a subset of targeted nodes.
/// </summary>
internal static class InverseFdmTargets
{
    /// <summary>
    /// Default weight of free nodes without a target toward their current
    /// position (dimensionless; one unit of displacement of a held node costs
    /// <c>w²</c> of a targeted node's).
    /// </summary>
    internal const double DefaultHoldWeight = 0.01;

    /// <summary>
    /// Complete a target over all free nodes.
    /// </summary>
    /// <param name="freePositions">Current free-node positions, in free-node order.</param>
    /// <param name="targetPoints">
    /// Desired positions: one per free node when <paramref name="targetNodeIndices"/>
    /// is null or empty, otherwise one per listed node.
    /// </param>
    /// <param name="targetNodeIndices">
    /// Free-node indices the target points refer to, or null/empty for a full target.
    /// </param>
    /// <param name="holdWeight">
    /// Weight of the untargeted rows, in <c>(0, 1]</c>. Ignored for a full target.
    /// </param>
    /// <exception cref="ArgumentException">Counts mismatch, duplicate or out-of-range node, invalid weight.</exception>
    internal static InverseFdmTargetSet Complete(
        IReadOnlyList<Point3d> freePositions,
        IReadOnlyList<Point3d> targetPoints,
        IReadOnlyList<int>? targetNodeIndices,
        double holdWeight)
    {
        int numFree = freePositions.Count;
        if (targetNodeIndices is null || targetNodeIndices.Count == 0)
        {
            if (targetPoints.Count != numFree)
                throw new ArgumentException(
                    $"Target points count ({targetPoints.Count}) must match free node count ({numFree}) "
                    + "when no Target Nodes are given.");
            return new InverseFdmTargetSet(
                Pack(targetPoints), null, targetPoints, Enumerable.Range(0, numFree).ToArray());
        }

        if (targetPoints.Count != targetNodeIndices.Count)
            throw new ArgumentException(
                $"Target points count ({targetPoints.Count}) must match Target Nodes count ({targetNodeIndices.Count}).");
        if (!double.IsFinite(holdWeight) || holdWeight <= 0.0 || holdWeight > 1.0)
            throw new ArgumentException(
                $"Hold Weight must be in (0, 1], got {holdWeight}. Use a small positive value, not zero, "
                + "for free nodes without a target.");

        var completed = freePositions.ToArray();
        var weights = new double[numFree];
        Array.Fill(weights, holdWeight);
        var seen = new bool[numFree];
        var targeted = new List<int>(targetNodeIndices.Count);
        for (int k = 0; k < targetNodeIndices.Count; k++)
        {
            int idx = targetNodeIndices[k];
            if (idx < 0 || idx >= numFree)
                throw new ArgumentException($"Target node index {idx} is outside the free-node range [0, {numFree}).");
            if (seen[idx])
                throw new ArgumentException($"Target node {idx} is listed more than once.");
            seen[idx] = true;
            completed[idx] = targetPoints[k];
            weights[idx] = 1.0;
            targeted.Add(idx);
        }
        targeted.Sort();
        return new InverseFdmTargetSet(Pack(completed), weights, completed, targeted);
    }

    /// <summary>
    /// <c>‖x_S(q) − x*_S‖</c> over the targeted free nodes: the plain landing
    /// distance, without the hold weights.
    /// </summary>
    /// <param name="solvedNodes">Solved positions indexed by graph node.</param>
    /// <param name="freeNodes">Graph index of each free node.</param>
    /// <param name="targets">The completed target and its targeted subset.</param>
    internal static double TargetError(
        IReadOnlyList<Point3d> solvedNodes,
        IReadOnlyList<int> freeNodes,
        InverseFdmTargetSet targets)
    {
        double sum = 0.0;
        foreach (int i in targets.TargetedFree)
        {
            double d = solvedNodes[freeNodes[i]].DistanceTo(targets.CompletedPoints[i]);
            sum += d * d;
        }
        return Math.Sqrt(sum);
    }

    private static double[] Pack(IReadOnlyList<Point3d> points)
    {
        var xyz = new double[points.Count * 3];
        for (int i = 0; i < points.Count; i++)
        {
            xyz[i * 3 + 0] = points[i].X;
            xyz[i * 3 + 1] = points[i].Y;
            xyz[i * 3 + 2] = points[i].Z;
        }
        return xyz;
    }
}
