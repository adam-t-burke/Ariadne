namespace Ariadne.Solver;

using System.Collections.Generic;
using Ariadne.Graphs;
using Theseus.Interop;

/// <summary>
/// Minimize variation in edge lengths.
/// </summary>
public sealed class LengthVariationObjective : EdgeObjective
{
    public LengthVariationObjective(double weight, List<Edge>? edges = null, double sharpness = 20.0)
    {
        Weight = weight;
        TargetEdges = edges;
        Sharpness = sharpness;
    }

    public override void ApplyTo(TheseusSolver solver, SolverContext context)
    {
        int[] indices = context.ResolveEdgeIndices(TargetEdges);
        solver.AddLengthVariation(Weight, indices, Sharpness);
    }
}

/// <summary>
/// Minimize variation in member forces.
/// </summary>
public sealed class ForceVariationObjective : EdgeObjective
{
    public ForceVariationObjective(double weight, List<Edge>? edges = null, double sharpness = 20.0)
    {
        Weight = weight;
        TargetEdges = edges;
        Sharpness = sharpness;
    }

    public override void ApplyTo(TheseusSolver solver, SolverContext context)
    {
        int[] indices = context.ResolveEdgeIndices(TargetEdges);
        solver.AddForceVariation(Weight, indices, Sharpness);
    }
}

/// <summary>
/// Minimize sum of force x length products (structural efficiency).
/// </summary>
public sealed class PerformanceObjective : EdgeObjective
{
    public PerformanceObjective(double weight, List<Edge>? edges = null)
    {
        Weight = weight;
        TargetEdges = edges;
    }

    public override void ApplyTo(TheseusSolver solver, SolverContext context)
    {
        int[] indices = context.ResolveEdgeIndices(TargetEdges);
        solver.AddSumForceLength(Weight, indices);
    }
}

/// <summary>
/// Target specific edge lengths.
/// </summary>
public sealed class TargetLengthObjective : ThresholdEdgeObjective
{
    public TargetLengthObjective(double weight, List<double> targets, List<Edge>? edges = null)
    {
        Weight = weight;
        Thresholds = targets;
        TargetEdges = edges;
    }

    public override void ApplyTo(TheseusSolver solver, SolverContext context)
    {
        int[] indices = context.ResolveEdgeIndices(TargetEdges);
        double[] targets = GetExpandedThresholds(indices.Length);
        solver.AddTargetLength(Weight, indices, targets);
    }
}

/// <summary>
/// Barrier penalty for edges below minimum length threshold.
/// </summary>
public sealed class MinLengthObjective : ThresholdEdgeObjective
{
    public MinLengthObjective(double weight, List<double> thresholds, List<Edge>? edges = null, double sharpness = 10.0)
    {
        Weight = weight;
        Thresholds = thresholds;
        TargetEdges = edges;
        Sharpness = sharpness;
    }

    public override void ApplyTo(TheseusSolver solver, SolverContext context)
    {
        int[] indices = context.ResolveEdgeIndices(TargetEdges);
        double[] thresholds = GetExpandedThresholds(indices.Length);
        solver.AddMinLength(Weight, indices, thresholds, Sharpness);
    }
}

/// <summary>
/// Barrier penalty for edges above maximum length threshold.
/// </summary>
public sealed class MaxLengthObjective : ThresholdEdgeObjective
{
    public MaxLengthObjective(double weight, List<double> thresholds, List<Edge>? edges = null, double sharpness = 10.0)
    {
        Weight = weight;
        Thresholds = thresholds;
        TargetEdges = edges;
        Sharpness = sharpness;
    }

    public override void ApplyTo(TheseusSolver solver, SolverContext context)
    {
        int[] indices = context.ResolveEdgeIndices(TargetEdges);
        double[] thresholds = GetExpandedThresholds(indices.Length);
        solver.AddMaxLength(Weight, indices, thresholds, Sharpness);
    }
}

/// <summary>
/// Barrier penalty for forces below minimum threshold.
/// </summary>
public sealed class MinForceObjective : ThresholdEdgeObjective
{
    public MinForceObjective(double weight, List<double> thresholds, List<Edge>? edges = null, double sharpness = 10.0)
    {
        Weight = weight;
        Thresholds = thresholds;
        TargetEdges = edges;
        Sharpness = sharpness;
    }

    public override void ApplyTo(TheseusSolver solver, SolverContext context)
    {
        int[] indices = context.ResolveEdgeIndices(TargetEdges);
        double[] thresholds = GetExpandedThresholds(indices.Length);
        solver.AddMinForce(Weight, indices, thresholds, Sharpness);
    }
}

/// <summary>
/// Barrier penalty for forces above maximum threshold.
/// </summary>
public sealed class MaxForceObjective : ThresholdEdgeObjective
{
    public MaxForceObjective(double weight, List<double> thresholds, List<Edge>? edges = null, double sharpness = 10.0)
    {
        Weight = weight;
        Thresholds = thresholds;
        TargetEdges = edges;
        Sharpness = sharpness;
    }

    public override void ApplyTo(TheseusSolver solver, SolverContext context)
    {
        int[] indices = context.ResolveEdgeIndices(TargetEdges);
        double[] thresholds = GetExpandedThresholds(indices.Length);
        solver.AddMaxForce(Weight, indices, thresholds, Sharpness);
    }
}
