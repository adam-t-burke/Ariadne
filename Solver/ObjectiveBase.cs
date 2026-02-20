namespace Ariadne.Solver;

using System;
using System.Collections.Generic;
using System.Linq;
using Ariadne.FDM;
using Ariadne.Graphs;
using Theseus.Interop;

/// <summary>
/// Base class for all objectives. Each objective knows how to apply itself to the Theseus solver.
/// </summary>
public abstract class Objective
{
    public double Weight { get; init; }
    public bool IsValid { get; protected set; } = true;

    public abstract void ApplyTo(TheseusSolver solver, SolverContext context);
}

/// <summary>
/// Context passed to objectives containing resolved indices and network data.
/// </summary>
public sealed class SolverContext
{
    public required FDM_Network Network { get; init; }
    public required Dictionary<Node, int> NodeIndexMap { get; init; }
    public required Dictionary<Edge, int> EdgeIndexMap { get; init; }

    /// <summary>
    /// Resolves node references to their indices. Returns all free nodes if list is null/empty.
    /// </summary>
    public int[] ResolveNodeIndices(List<Node>? nodes)
    {
        if (nodes == null || nodes.Count == 0)
            return [.. Network.FreeNodes];

        return nodes.Select(n => NodeIndexMap[n]).ToArray();
    }

    /// <summary>
    /// Resolves edge references to their indices. Returns all edges if list is null/empty.
    /// </summary>
    public int[] ResolveEdgeIndices(List<Edge>? edges)
    {
        if (edges == null || edges.Count == 0)
            return Enumerable.Range(0, Network.Graph.Ne).ToArray();

        return edges.Select(e => EdgeIndexMap[e]).ToArray();
    }

    public Rhino.Geometry.Point3d GetNodePosition(int nodeIndex) =>
        Network.Graph.Nodes[nodeIndex].Value;
}

/// <summary>
/// Base class for node-based objectives.
/// </summary>
public abstract class NodeObjective : Objective
{
    public List<Node>? TargetNodes { get; init; }
}

/// <summary>
/// Base class for edge-based objectives.
/// </summary>
public abstract class EdgeObjective : Objective
{
    public List<Edge>? TargetEdges { get; init; }
    public double Sharpness { get; init; } = 20.0;
}

/// <summary>
/// Base class for edge objectives that have threshold values.
/// </summary>
public abstract class ThresholdEdgeObjective : EdgeObjective
{
    public List<double> Thresholds { get; init; } = [];

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
