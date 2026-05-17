namespace Ariadne.Solver;

using System;
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
    /// <summary>Plane to project target positions onto; null uses world XY.</summary>
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

    public override int GetContentHashCode()
    {
        var h = new HashCode();
        h.Add(base.GetContentHashCode());
        if (Plane is { } p)
        {
            h.Add(p.Origin.X); h.Add(p.Origin.Y); h.Add(p.Origin.Z);
            h.Add(p.XAxis.X); h.Add(p.XAxis.Y); h.Add(p.XAxis.Z);
            h.Add(p.YAxis.X); h.Add(p.YAxis.Y); h.Add(p.YAxis.Z);
        }
        else
            h.Add(0);
        return h.ToHashCode();
    }
}

/// <summary>
/// Planar constraint: pull nodes onto a plane along a given direction. No target positions —
/// at each step minimizes squared distance (along the direction) from each node to the plane.
/// When Direction is null, uses the plane normal (orthographic pull onto the plane).
/// </summary>
public sealed class PlanarConstraintAlongDirectionObjective : NodeObjective
{
    /// <summary>Plane to constrain nodes to; null uses world XY.</summary>
    public Plane? Plane { get; }
    /// <summary>Direction along which distance to the plane is measured; null uses plane normal.</summary>
    public Vector3d? Direction { get; }

    public PlanarConstraintAlongDirectionObjective(double weight, List<Node>? nodes = null, Plane? plane = null, Vector3d? direction = null)
    {
        Weight = weight;
        TargetNodes = nodes;
        Plane = plane;
        Direction = direction;
    }

    public override void ApplyTo(TheseusSolver solver, SolverContext context)
    {
        var plane = Plane ?? Rhino.Geometry.Plane.WorldXY;
        Vector3d dir = Direction ?? plane.ZAxis;
        double nDotD = Math.Abs(plane.ZAxis.X * dir.X + plane.ZAxis.Y * dir.Y + plane.ZAxis.Z * dir.Z);
        if (nDotD < 1e-6)
        {
            IsValid = false;
            return;
        }

        int[] indices = context.ResolveNodeIndices(TargetNodes);
        double[] origin = { plane.Origin.X, plane.Origin.Y, plane.Origin.Z };
        double[] xAxis = { plane.XAxis.X, plane.XAxis.Y, plane.XAxis.Z };
        double[] yAxis = { plane.YAxis.X, plane.YAxis.Y, plane.YAxis.Z };
        double[] directionArr = { dir.X, dir.Y, dir.Z };

        solver.AddPlanarConstraintAlongDirection(Weight, indices, origin, xAxis, yAxis, directionArr);
    }

    public override int GetContentHashCode()
    {
        var h = new HashCode();
        h.Add(base.GetContentHashCode());
        if (Plane is { } p)
        {
            h.Add(p.Origin.X); h.Add(p.Origin.Y); h.Add(p.Origin.Z);
            h.Add(p.XAxis.X); h.Add(p.XAxis.Y); h.Add(p.XAxis.Z);
            h.Add(p.YAxis.X); h.Add(p.YAxis.Y); h.Add(p.YAxis.Z);
        }
        else
            h.Add(0);
        if (Direction is { } d)
        {
            h.Add(d.X); h.Add(d.Y); h.Add(d.Z);
        }
        else
            h.Add(0);
        return h.ToHashCode();
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
public enum ReactionMagnitudeBehavior
{
    Target = 0,
    Max = 1,
    Min = 2,
}

public enum ReactionMagnitudeSign
{
    Unsigned = 0,
    SignedProjected = 1,
}

public sealed class ReactionObjective : NodeObjective
{
    public List<Vector3d> TargetDirections { get; init; }
    public List<double>? TargetMagnitudes { get; init; }
    public bool IncludeDirection { get; init; }
    public bool IncludeMagnitude { get; init; }
    public ReactionMagnitudeBehavior MagnitudeBehavior { get; init; }
    public ReactionMagnitudeSign MagnitudeSign { get; init; }

    public ReactionObjective(
        double weight,
        List<Node> anchorNodes,
        List<Vector3d> targetDirections,
        bool includeDirection = true,
        bool includeMagnitude = false,
        List<double>? targetMagnitudes = null,
        ReactionMagnitudeBehavior magnitudeBehavior = ReactionMagnitudeBehavior.Max,
        ReactionMagnitudeSign magnitudeSign = ReactionMagnitudeSign.Unsigned)
    {
        Weight = weight;
        TargetNodes = anchorNodes;
        TargetDirections = targetDirections;
        IncludeDirection = includeDirection;
        IncludeMagnitude = includeMagnitude;
        TargetMagnitudes = targetMagnitudes;
        MagnitudeBehavior = magnitudeBehavior;
        MagnitudeSign = magnitudeSign;
    }

    public override void ApplyTo(TheseusSolver solver, SolverContext context)
    {
        if (TargetNodes == null || TargetNodes.Count == 0) return;
        if (!IncludeDirection && !IncludeMagnitude) return;
        if (IncludeDirection && TargetDirections.Count == 0)
        {
            IsValid = false;
            return;
        }
        if (IncludeMagnitude && TargetMagnitudes is not { Count: > 0 })
        {
            IsValid = false;
            return;
        }

        int[] indices = TargetNodes.Select(n => context.NodeIndexMap[n]).ToArray();

        double[] targetDirs = new double[indices.Length * 3];
        for (int i = 0; i < indices.Length; i++)
        {
            var dir = TargetDirections.Count > 0
                ? (i < TargetDirections.Count ? TargetDirections[i] : TargetDirections[^1])
                : Vector3d.XAxis;
            targetDirs[i * 3 + 0] = dir.X;
            targetDirs[i * 3 + 1] = dir.Y;
            targetDirs[i * 3 + 2] = dir.Z;
        }

        if (MagnitudeSign == ReactionMagnitudeSign.SignedProjected && !HasUsableDirections(targetDirs))
        {
            IsValid = false;
            return;
        }

        if (IncludeMagnitude)
        {
            var magnitudes = TargetMagnitudes!;
            double[] targetMags = new double[indices.Length];
            for (int i = 0; i < indices.Length; i++)
                targetMags[i] = i < magnitudes.Count ? magnitudes[i] : magnitudes[^1];

            if (IncludeDirection)
            {
                solver.AddReactionDirectionMagnitude(
                    Weight,
                    indices,
                    targetDirs,
                    targetMags,
                    (int)MagnitudeBehavior,
                    (int)MagnitudeSign);
            }
            else
            {
                solver.AddReactionMagnitude(
                    Weight,
                    indices,
                    targetDirs,
                    targetMags,
                    (int)MagnitudeBehavior,
                    (int)MagnitudeSign);
            }
        }
        else
        {
            solver.AddReactionDirection(Weight, indices, targetDirs);
        }
    }

    public override int GetContentHashCode()
    {
        var h = new HashCode();
        h.Add(base.GetContentHashCode());
        h.Add(IncludeDirection);
        h.Add(IncludeMagnitude);
        h.Add(MagnitudeBehavior);
        h.Add(MagnitudeSign);
        if (TargetDirections != null)
            foreach (var d in TargetDirections)
            { h.Add(d.X); h.Add(d.Y); h.Add(d.Z); }
        if (TargetMagnitudes != null)
            foreach (var m in TargetMagnitudes)
                h.Add(m);
        return h.ToHashCode();
    }

    private static bool HasUsableDirections(double[] targetDirs)
    {
        for (int i = 0; i < targetDirs.Length; i += 3)
        {
            double norm = Math.Sqrt(
                targetDirs[i + 0] * targetDirs[i + 0] +
                targetDirs[i + 1] * targetDirs[i + 1] +
                targetDirs[i + 2] * targetDirs[i + 2]);
            if (norm < 1e-10)
                return false;
        }
        return true;
    }
}
