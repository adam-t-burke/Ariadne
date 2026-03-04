namespace Ariadne.FEA;

using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.CompilerServices;
using Ariadne.Graphs;
using Theseus.Interop;

/// <summary>
/// Base class for FEA objectives. Each objective registers itself with the native FEA solver.
/// </summary>
public abstract class FeaObjective
{
    public double Weight { get; init; } = 1.0;
    public bool IsValid { get; protected set; } = true;

    public abstract void ApplyTo(FeaSolver solver, FeaSolverContext context);

    public virtual int GetContentHashCode()
    {
        var h = new HashCode();
        h.Add(Weight);
        h.Add(GetType());
        return h.ToHashCode();
    }
}

/// <summary>
/// Context for resolving node/edge references to solver indices.
/// </summary>
public sealed class FeaSolverContext
{
    public required FEA_Network Network { get; init; }
    public required Dictionary<Node, int> NodeIndexMap { get; init; }
    public required Dictionary<Edge, int> EdgeIndexMap { get; init; }

    public int[] ResolveNodeIndices(List<Node>? nodes)
    {
        if (nodes == null || nodes.Count == 0)
            return Enumerable.Range(0, Network.Graph.Nn).ToArray();
        return nodes.Select(n => NodeIndexMap[n]).ToArray();
    }

    public int[] ResolveEdgeIndices(List<Edge>? edges)
    {
        if (edges == null || edges.Count == 0)
            return Enumerable.Range(0, Network.Graph.Ne).ToArray();
        return edges.Select(e => EdgeIndexMap[e]).ToArray();
    }
}

/// <summary>
/// Node-based FEA objective.
/// </summary>
public abstract class FeaNodeObjective : FeaObjective
{
    public List<Node>? TargetNodes { get; init; }

    public override int GetContentHashCode()
    {
        var h = new HashCode();
        h.Add(Weight);
        h.Add(GetType());
        if (TargetNodes == null || TargetNodes.Count == 0)
            h.Add(0);
        else
            foreach (var n in TargetNodes)
                h.Add(RuntimeHelpers.GetHashCode(n));
        return h.ToHashCode();
    }
}

/// <summary>
/// Edge-based FEA objective.
/// </summary>
public abstract class FeaEdgeObjective : FeaObjective
{
    public List<Edge>? TargetEdges { get; init; }
    public double Sharpness { get; init; } = 20.0;

    public override int GetContentHashCode()
    {
        var h = new HashCode();
        h.Add(Weight);
        h.Add(Sharpness);
        h.Add(GetType());
        if (TargetEdges == null || TargetEdges.Count == 0)
            h.Add(0);
        else
            foreach (var e in TargetEdges)
                h.Add(RuntimeHelpers.GetHashCode(e));
        return h.ToHashCode();
    }
}
