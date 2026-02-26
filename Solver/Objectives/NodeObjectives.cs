namespace Ariadne.Solver;

using System.Collections.Generic;
using System.Linq;
using Ariadne.Graphs;
using Rhino.Geometry;
using Theseus.Interop;

/// <summary>
/// Minimize XY deviation from initial node positions (TNA-like behavior).
/// </summary>
public sealed class TargetXYObjective : NodeObjective
{
    public TargetXYObjective(double weight, List<Node>? nodes = null)
    {
        Weight = weight;
        TargetNodes = nodes;
    }

    public override void ApplyTo(TheseusSolver solver, SolverContext context)
    {
        int[] indices = context.ResolveNodeIndices(TargetNodes);

        // Rust FFI expects num_nodes * 3 (row-major X,Y,Z); Z is ignored by the loss but required for layout.
        double[] targetXy = new double[indices.Length * 3];
        for (int i = 0; i < indices.Length; i++)
        {
            var pos = context.GetNodePosition(indices[i]);
            targetXy[i * 3 + 0] = pos.X;
            targetXy[i * 3 + 1] = pos.Y;
            targetXy[i * 3 + 2] = pos.Z;
        }

        solver.AddTargetXy(Weight, indices, targetXy);
    }
}

/// <summary>
/// Minimize 3D deviation from initial node positions.
/// </summary>
public sealed class TargetXYZObjective : NodeObjective
{
    public TargetXYZObjective(double weight, List<Node>? nodes = null)
    {
        Weight = weight;
        TargetNodes = nodes;
    }

    public override void ApplyTo(TheseusSolver solver, SolverContext context)
    {
        int[] indices = context.ResolveNodeIndices(TargetNodes);

        double[] targetXyz = new double[indices.Length * 3];
        for (int i = 0; i < indices.Length; i++)
        {
            var pos = context.GetNodePosition(indices[i]);
            targetXyz[i * 3 + 0] = pos.X;
            targetXyz[i * 3 + 1] = pos.Y;
            targetXyz[i * 3 + 2] = pos.Z;
        }

        solver.AddTargetXyz(Weight, indices, targetXyz);
    }
}

/// <summary>
/// Minimize deviation from initial node positions projected onto an arbitrary plane.
/// Defaults to the world XY plane when no plane is supplied.
/// </summary>
public sealed class TargetPlaneObjective : NodeObjective
{
    public Plane? Plane { get; }

    public TargetPlaneObjective(double weight, List<Node>? nodes = null, Plane? plane = null)
    {
        Weight = weight;
        TargetNodes = nodes;
        Plane = plane;
    }

    public override void ApplyTo(TheseusSolver solver, SolverContext context)
    {
        int[] indices = context.ResolveNodeIndices(TargetNodes);

        double[] targetXyz = new double[indices.Length * 3];
        for (int i = 0; i < indices.Length; i++)
        {
            var pos = context.GetNodePosition(indices[i]);
            targetXyz[i * 3 + 0] = pos.X;
            targetXyz[i * 3 + 1] = pos.Y;
            targetXyz[i * 3 + 2] = pos.Z;
        }

        var plane = Plane ?? Rhino.Geometry.Plane.WorldXY;
        double[] origin = { plane.Origin.X, plane.Origin.Y, plane.Origin.Z };
        double[] xAxis = { plane.XAxis.X, plane.XAxis.Y, plane.XAxis.Z };
        double[] yAxis = { plane.YAxis.X, plane.YAxis.Y, plane.YAxis.Z };

        solver.AddTargetPlane(Weight, indices, targetXyz, origin, xAxis, yAxis);
    }
}

/// <summary>
/// Maintain relative distances within a point set (rigid body constraint).
/// </summary>
public sealed class RigidPointSetObjective : NodeObjective
{
    public RigidPointSetObjective(double weight, List<Node>? nodes = null)
    {
        Weight = weight;
        TargetNodes = nodes;
    }

    public override void ApplyTo(TheseusSolver solver, SolverContext context)
    {
        int[] indices = context.ResolveNodeIndices(TargetNodes);

        double[] targetXyz = new double[indices.Length * 3];
        for (int i = 0; i < indices.Length; i++)
        {
            var pos = context.GetNodePosition(indices[i]);
            targetXyz[i * 3 + 0] = pos.X;
            targetXyz[i * 3 + 1] = pos.Y;
            targetXyz[i * 3 + 2] = pos.Z;
        }

        solver.AddRigidSetCompare(Weight, indices, targetXyz);
    }
}

/// <summary>
/// Align anchor reaction force directions with target directions,
/// optionally also matching target magnitudes.
/// </summary>
public sealed class ReactionObjective : NodeObjective
{
    public List<Vector3d> TargetDirections { get; init; }
    public List<double>? TargetMagnitudes { get; init; }
    public bool IncludeMagnitude { get; init; }

    public ReactionObjective(
        double weight,
        List<Node> anchorNodes,
        List<Vector3d> targetDirections,
        bool includeMagnitude = false,
        List<double>? targetMagnitudes = null)
    {
        Weight = weight;
        TargetNodes = anchorNodes;
        TargetDirections = targetDirections;
        IncludeMagnitude = includeMagnitude;
        TargetMagnitudes = targetMagnitudes;
    }

    public override void ApplyTo(TheseusSolver solver, SolverContext context)
    {
        if (TargetNodes == null || TargetNodes.Count == 0) return;

        int[] indices = TargetNodes.Select(n => context.NodeIndexMap[n]).ToArray();

        double[] targetDirs = new double[indices.Length * 3];
        for (int i = 0; i < indices.Length; i++)
        {
            var dir = i < TargetDirections.Count ? TargetDirections[i] : TargetDirections[^1];
            targetDirs[i * 3 + 0] = dir.X;
            targetDirs[i * 3 + 1] = dir.Y;
            targetDirs[i * 3 + 2] = dir.Z;
        }

        if (IncludeMagnitude && TargetMagnitudes is { Count: > 0 })
        {
            double[] targetMags = new double[indices.Length];
            for (int i = 0; i < indices.Length; i++)
                targetMags[i] = i < TargetMagnitudes.Count ? TargetMagnitudes[i] : TargetMagnitudes[^1];

            solver.AddReactionDirectionMagnitude(Weight, indices, targetDirs, targetMags);
        }
        else
        {
            solver.AddReactionDirection(Weight, indices, targetDirs);
        }
    }
}
