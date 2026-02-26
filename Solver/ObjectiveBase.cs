namespace Ariadne.Solver;

using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using Ariadne.FDM;
using Ariadne.Graphs;
using Theseus.Interop;

/// <summary>
/// Base class for all objectives. Each objective knows how to apply itself to the Theseus solver.
/// </summary>
public abstract class Objective
{
    /// <summary>Weight of this objective in the combined loss (default 1.0).</summary>
    public double Weight { get; init; }
    /// <summary>False if the objective could not be built (e.g. invalid inputs); such objectives are skipped.</summary>
    public bool IsValid { get; protected set; } = true;

    /// <summary>Registers this objective with the native solver (e.g. target length, force variation).</summary>
    /// <param name="solver">The Theseus solver instance to configure.</param>
    /// <param name="context">Resolved node/edge indices and network data.</param>
    public abstract void ApplyTo(TheseusSolver solver, SolverContext context);

    /// <summary>Content hash for cache invalidation; must reflect node/edge subset and parameters.</summary>
    public virtual int GetContentHashCode()
    {
        var h = new HashCode();
        h.Add(Weight);
        h.Add(GetType());
        return h.ToHashCode();
    }
}

/// <summary>
/// Context passed to objectives containing resolved indices and network data.
/// </summary>
public sealed class SolverContext
{
    /// <summary>The FDM network being solved.</summary>
    public required FDM_Network Network { get; init; }
    /// <summary>Maps each graph node to its 0-based index in the solver.</summary>
    public required Dictionary<Node, int> NodeIndexMap { get; init; }
    /// <summary>Maps each graph edge to its 0-based index in the solver.</summary>
    public required Dictionary<Edge, int> EdgeIndexMap { get; init; }

    /// <summary>
    /// Resolves node references to their indices. Returns all free nodes if list is null/empty.
    /// </summary>
    /// <param name="nodes">Nodes to resolve, or null/empty for all free nodes.</param>
    /// <returns>Array of 0-based node indices.</returns>
    public int[] ResolveNodeIndices(List<Node>? nodes)
    {
        if (nodes == null || nodes.Count == 0)
            return [.. Network.FreeNodes];

        return nodes.Select(n => NodeIndexMap[n]).ToArray();
    }

    /// <summary>
    /// Resolves edge references to their indices. Returns all edges if list is null/empty.
    /// </summary>
    /// <param name="edges">Edges to resolve, or null/empty for all edges.</param>
    /// <returns>Array of 0-based edge indices.</returns>
    public int[] ResolveEdgeIndices(List<Edge>? edges)
    {
        if (edges == null || edges.Count == 0)
            return Enumerable.Range(0, Network.Graph.Ne).ToArray();

        return edges.Select(e => EdgeIndexMap[e]).ToArray();
    }

    /// <summary>Gets the 3D position of a node by its index.</summary>
    /// <param name="nodeIndex">0-based node index in Graph.Nodes.</param>
    /// <returns>The node position.</returns>
    public Rhino.Geometry.Point3d GetNodePosition(int nodeIndex) =>
        Network.Graph.Nodes[nodeIndex].Value;
}

/// <summary>
/// Base class for node-based objectives.
/// </summary>
public abstract class NodeObjective : Objective
{
    /// <summary>Nodes to apply the objective to; null or empty means all free nodes.</summary>
    public List<Node>? TargetNodes { get; init; }

    /// <inheritdoc />
    public override int GetContentHashCode()
    {
        var h = new HashCode();
        h.Add(Weight);
        h.Add(GetType());
        if (TargetNodes == null || TargetNodes.Count == 0)
            h.Add(0); // sentinel for "all nodes"
        else
            foreach (var n in TargetNodes)
                h.Add(RuntimeHelpers.GetHashCode(n));
        return h.ToHashCode();
    }
}

/// <summary>
/// Base class for edge-based objectives.
/// </summary>
public abstract class EdgeObjective : Objective
{
    /// <summary>Edges to apply the objective to; null or empty means all edges.</summary>
    public List<Edge>? TargetEdges { get; init; }
    /// <summary>Sharpness parameter for barrier/smoothing (e.g. in min/max constraints).</summary>
    public double Sharpness { get; init; } = 20.0;

    /// <inheritdoc />
    public override int GetContentHashCode()
    {
        var h = new HashCode();
        h.Add(Weight);
        h.Add(Sharpness);
        h.Add(GetType());
        if (TargetEdges == null || TargetEdges.Count == 0)
            h.Add(0); // sentinel for "all edges"
        else
            foreach (var e in TargetEdges)
                h.Add(RuntimeHelpers.GetHashCode(e));
        return h.ToHashCode();
    }
}

/// <summary>
/// Base class for edge objectives that have threshold values.
/// </summary>
public abstract class ThresholdEdgeObjective : EdgeObjective
{
    /// <summary>Threshold values per edge (or one value applied to all).</summary>
    public List<double> Thresholds { get; init; } = [];

    /// <inheritdoc />
    public override int GetContentHashCode()
    {
        var h = new HashCode();
        h.Add(base.GetContentHashCode());
        foreach (var t in Thresholds)
            h.Add(t);
        return h.ToHashCode();
    }

    /// <summary>Expands Thresholds to an array of length count (repeating last value if needed).</summary>
    /// <param name="count">Required length.</param>
    /// <returns>Array of threshold values.</returns>
    protected double[] GetExpandedThresholds(int count)
    {
        if (Thresholds.Count == 0)
            throw new InvalidOperationException("Thresholds cannot be empty");

        var result = new double[count];
        for (int i = 0; i < count; i++)
            result[i] = i < Thresholds.Count ? Thresholds[i] : Thresholds[^1];
        return result;
    }
}
