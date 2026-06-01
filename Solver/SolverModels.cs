namespace Ariadne.Solver;

using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using Ariadne.FDM;
using Ariadne.Graphs;
using Rhino.Geometry;

/// <summary>
/// Solver configuration options.
/// </summary>
public sealed record SolverOptions
{
    /// <summary>Maximum optimization iterations.</summary>
    public int MaxIterations { get; init; } = 500;
    /// <summary>Absolute convergence tolerance.</summary>
    public double AbsTol { get; init; } = 1e-6;
    /// <summary>Relative convergence tolerance.</summary>
    public double RelTol { get; init; } = 1e-6;
    /// <summary>Barrier function weight for bound constraints.</summary>
    public double BarrierWeight { get; init; } = 10.0;
    /// <summary>Barrier function sharpness.</summary>
    public double BarrierSharpness { get; init; } = 10.0;
    /// <summary>Invoke progress callback every N accepted L-BFGS iterations (0 = every iteration).</summary>
    public int ReportFrequency { get; init; } = 10;
}

public enum QParameterizationMode
{
    DirectSoftBounds = 0,
    ImplicitBounded = 1,
}

/// <summary>
/// Bundled optimization configuration passed from the OptConfig component
/// to the solve component. When absent, the solver runs forward-only.
/// </summary>
public sealed record OptimizationConfig
{
    /// <summary>Objective functions to minimize (e.g. target length, force variation).</summary>
    public required IReadOnlyList<Objective> Objectives { get; init; }
    /// <summary>Lower bounds on force densities per edge.</summary>
    public IReadOnlyList<double> LowerBounds { get; init; } = [0.1];
    /// <summary>Upper bounds on force densities per edge.</summary>
    public IReadOnlyList<double> UpperBounds { get; init; } = [100.0];
    /// <summary>Maximum optimization iterations.</summary>
    public int MaxIterations { get; init; } = 500;
    /// <summary>Absolute convergence tolerance.</summary>
    public double AbsTol { get; init; } = 1e-6;
    /// <summary>Relative convergence tolerance.</summary>
    public double RelTol { get; init; } = 1e-6;
    /// <summary>Barrier function weight.</summary>
    public double BarrierWeight { get; init; } = 10.0;
    /// <summary>Barrier function sharpness.</summary>
    public double BarrierSharpness { get; init; } = 10.0;
    /// <summary>Progress callback frequency (accepted L-BFGS iterations between callbacks).</summary>
    public int ReportFrequency { get; init; } = 10;
    /// <summary>q optimization mode: direct soft bounds (default) or implicit hard-bounded mapping.</summary>
    public QParameterizationMode QParameterizationMode { get; init; } = QParameterizationMode.DirectSoftBounds;
    /// <summary>When true, optimization runs (e.g. from a button or toggle).</summary>
    public bool Run { get; init; } = false;
    /// <summary>When true, stream intermediate results to outputs during optimization.</summary>
    public bool StreamPreview { get; init; } = true;
    /// <summary>Variable support definitions (optional).</summary>
    public IReadOnlyList<VariableSupportConfig> VariableSupports { get; init; } = [];
}

/// <summary>
/// Base definition for variable supports attached to fixed nodes.
/// </summary>
public abstract record VariableSupportConfig
{
    /// <summary>Fixed nodes this support definition applies to.</summary>
    public required List<Node> Nodes { get; init; }

    public virtual int GetContentHashCode()
    {
        var h = new HashCode();
        h.Add(GetType());
        foreach (var node in Nodes)
            h.Add(RuntimeHelpers.GetHashCode(node));
        return h.ToHashCode();
    }
}

/// <summary>
/// Spherical support around each node's initial position.
/// </summary>
public sealed record SphereVariableSupport : VariableSupportConfig
{
    /// <summary>Maximum relative radius from initial node position.</summary>
    public required double Radius { get; init; }

    public override int GetContentHashCode()
    {
        var h = new HashCode();
        h.Add(base.GetContentHashCode());
        h.Add(Radius);
        return h.ToHashCode();
    }
}

/// <summary>
/// Roller-style support with per-axis relative bounds from initial node position.
/// Disabled axes remain fixed.
/// </summary>
public sealed record RollerVariableSupport : VariableSupportConfig
{
    public bool FreeX { get; init; } = true;
    public bool FreeY { get; init; } = true;
    public bool FreeZ { get; init; } = true;
    public Interval DomainX { get; init; } = new(-1.0, 1.0);
    public Interval DomainY { get; init; } = new(-1.0, 1.0);
    public Interval DomainZ { get; init; } = new(-1.0, 1.0);

    public override int GetContentHashCode()
    {
        var h = new HashCode();
        h.Add(base.GetContentHashCode());
        h.Add(FreeX);
        h.Add(FreeY);
        h.Add(FreeZ);
        h.Add(DomainX.T0);
        h.Add(DomainX.T1);
        h.Add(DomainY.T0);
        h.Add(DomainY.T1);
        h.Add(DomainZ.T0);
        h.Add(DomainZ.T1);
        return h.ToHashCode();
    }
}

/// <summary>
/// Rail support constrained to a provided line segment.
/// </summary>
public sealed record RailVariableSupport : VariableSupportConfig
{
    /// <summary>Line segment that constrains node movement.</summary>
    public required Line Rail { get; init; }

    public override int GetContentHashCode()
    {
        var h = new HashCode();
        h.Add(base.GetContentHashCode());
        h.Add(Rail.FromX);
        h.Add(Rail.FromY);
        h.Add(Rail.FromZ);
        h.Add(Rail.ToX);
        h.Add(Rail.ToY);
        h.Add(Rail.ToZ);
        return h.ToHashCode();
    }
}

/// <summary>
/// Self-weight configuration.
/// </summary>
public abstract record SelfWeightConfig
{
    /// <summary>Gravity vector (default [0, 0, -9.81]).</summary>
    public Vector3d Gravity { get; init; } = new(0, 0, -9.81);
    /// <summary>Maximum self-weight iterations.</summary>
    public int MaxIters { get; init; } = 50;
    /// <summary>Convergence tolerance for load iteration.</summary>
    public double Tolerance { get; init; } = 1e-6;
    /// <summary>Relaxation factor (1.0 = no relaxation).</summary>
    public double Relaxation { get; init; } = 1.0;

    public virtual int GetContentHashCode()
    {
        var h = new HashCode();
        h.Add(GetType());
        h.Add(Gravity.X);
        h.Add(Gravity.Y);
        h.Add(Gravity.Z);
        h.Add(MaxIters);
        h.Add(Tolerance);
        h.Add(Relaxation);
        return h.ToHashCode();
    }

    /// <summary>Prescribed linear density (mass/length) per edge.</summary>
    public sealed record Prescribed : SelfWeightConfig
    {
        /// <summary>Linear density per edge (kg/m). One value = uniform.</summary>
        public required List<double> LinearDensities { get; init; }

        public override int GetContentHashCode()
        {
            var h = new HashCode();
            h.Add(base.GetContentHashCode());
            foreach (var density in LinearDensities)
                h.Add(density);
            return h.ToHashCode();
        }
    }

    /// <summary>Force-based sizing: A_k = |F_k| / sigma.</summary>
    public sealed record Sizing : SelfWeightConfig
    {
        /// <summary>Material density (kg/m^3).</summary>
        public required double Rho { get; init; }
        /// <summary>Allowable stress (Pa or consistent units).</summary>
        public required double Sigma { get; init; }

        public override int GetContentHashCode()
        {
            var h = new HashCode();
            h.Add(base.GetContentHashCode());
            h.Add(Rho);
            h.Add(Sigma);
            return h.ToHashCode();
        }
    }
}

/// <summary>
/// Pressure load configuration.
/// </summary>
public abstract record PressureConfig
{
    /// <summary>Face topology: each inner list is ordered vertex indices for one face.</summary>
    public required List<List<int>> Faces { get; init; }
    /// <summary>Maximum pressure iteration count.</summary>
    public int MaxIters { get; init; } = 50;
    /// <summary>Convergence tolerance.</summary>
    public double Tolerance { get; init; } = 1e-6;
    /// <summary>Relaxation factor.</summary>
    public double Relaxation { get; init; } = 1.0;

    public virtual int GetContentHashCode()
    {
        var h = new HashCode();
        h.Add(GetType());
        h.Add(MaxIters);
        h.Add(Tolerance);
        h.Add(Relaxation);
        foreach (var face in Faces)
        {
            h.Add(face.Count);
            foreach (var vertex in face)
                h.Add(vertex);
        }
        return h.ToHashCode();
    }

    /// <summary>Constant pressure along each face's outward normal.</summary>
    public sealed record Normal : PressureConfig
    {
        /// <summary>Pressure magnitude per face (positive = along outward normal).</summary>
        public required List<double> Pressures { get; init; }

        public override int GetContentHashCode()
        {
            var h = new HashCode();
            h.Add(base.GetContentHashCode());
            foreach (var pressure in Pressures)
                h.Add(pressure);
            return h.ToHashCode();
        }
    }

    /// <summary>Hydrostatic pressure varying linearly with depth (e.g. fabric formwork).</summary>
    public sealed record Hydrostatic : PressureConfig
    {
        /// <summary>Fluid density (e.g. 2400 kg/m^3 for wet concrete).</summary>
        public required double RhoFluid { get; init; }
        /// <summary>Gravitational acceleration magnitude (default 9.81).</summary>
        public double GMagnitude { get; init; } = 9.81;
        /// <summary>Free surface elevation (top of pour) along the up direction.</summary>
        public required double ZDatum { get; init; }
        /// <summary>Unit "up" direction (default Z-axis). Depth is measured opposite to this.</summary>
        public Vector3d UpDirection { get; init; } = new(0, 0, 1);

        public override int GetContentHashCode()
        {
            var h = new HashCode();
            h.Add(base.GetContentHashCode());
            h.Add(RhoFluid);
            h.Add(GMagnitude);
            h.Add(ZDatum);
            h.Add(UpDirection.X);
            h.Add(UpDirection.Y);
            h.Add(UpDirection.Z);
            return h.ToHashCode();
        }
    }

    /// <summary>Directional pressure proportional to projected face area (e.g. dead load, soil pressure).</summary>
    public sealed record Directional : PressureConfig
    {
        /// <summary>Pressure magnitude per face.</summary>
        public required List<double> Pressures { get; init; }
        /// <summary>Unit load direction (e.g. [0,0,-1] for gravity dead load).</summary>
        public required Vector3d Direction { get; init; }

        public override int GetContentHashCode()
        {
            var h = new HashCode();
            h.Add(base.GetContentHashCode());
            foreach (var pressure in Pressures)
                h.Add(pressure);
            h.Add(Direction.X);
            h.Add(Direction.Y);
            h.Add(Direction.Z);
            return h.ToHashCode();
        }
    }
}

/// <summary>
/// Inputs required for the solver. Bounds are nullable: null means
/// unconstrained (forward-only), non-null means optimization bounds.
/// </summary>
public sealed record SolverInputs
{
    /// <summary>Initial force densities (one per edge).</summary>
    public required List<double> QInit { get; init; }
    /// <summary>Load vectors on free nodes (one per free node).</summary>
    public required List<Vector3d> Loads { get; init; }
    /// <summary>Indices into the free-node list that should receive loads (null = all free nodes).</summary>
    public List<int>? LoadNodeIndices { get; init; }
    /// <summary>Lower bounds on q (null = forward-only).</summary>
    public List<double>? LowerBounds { get; init; }
    /// <summary>Upper bounds on q (null = forward-only).</summary>
    public List<double>? UpperBounds { get; init; }
    /// <summary>Objectives to minimize when optimizing.</summary>
    public List<Objective> Objectives { get; init; } = [];
    /// <summary>q optimization parameterization mode.</summary>
    public QParameterizationMode QParameterizationMode { get; init; } = QParameterizationMode.DirectSoftBounds;
    /// <summary>Self-weight configuration (null = no self-weight).</summary>
    public SelfWeightConfig? SelfWeight { get; init; }
    /// <summary>Pressure load configuration (null = no pressure loads).</summary>
    public PressureConfig? Pressure { get; init; }
    /// <summary>Variable support definitions (null/empty = no variable anchors).</summary>
    public List<VariableSupportConfig>? VariableSupports { get; init; }
}

/// <summary>
/// Results from a solve operation.
/// </summary>
public sealed record SolveResult
{
    /// <summary>Network with updated node positions (and optionally updated q).</summary>
    public required FDM_Network Network { get; init; }
    /// <summary>Force densities (q) per edge after solve.</summary>
    public required double[] ForceDensities { get; init; }
    /// <summary>Member forces (q * length) per edge.</summary>
    public required double[] MemberForces { get; init; }
    /// <summary>Member lengths per edge.</summary>
    public required double[] MemberLengths { get; init; }
    /// <summary>Reaction forces at fixed nodes.</summary>
    public required double[] Reactions { get; init; }
    /// <summary>Number of solver iterations (0 for forward-only).</summary>
    public required int Iterations { get; init; }
    /// <summary>True if the optimizer converged (or N/A for forward-only).</summary>
    public required bool Converged { get; init; }
    /// <summary>Native optimizer termination reason when available.</summary>
    public string TerminationReason { get; init; } = "";

    /// <summary>
    /// Node positions as Point3d list (convenience accessor).
    /// </summary>
    public List<Point3d> NodePositions => Network.Graph.Nodes.ConvertAll(n => n.Value);

    /// <summary>
    /// Edge curves as LineCurve list (convenience accessor).
    /// </summary>
    public List<LineCurve> EdgeCurves => Network.Graph.Edges.ConvertAll(e => (LineCurve)e.Value);
}

/// <summary>
/// Internal data structure for solver creation.
/// </summary>
internal sealed record SolverData(
    int NumEdges,
    int NumNodes,
    int NumFree,
    int[] CooRows,
    int[] CooCols,
    double[] CooVals,
    int[] FreeIndices,
    int[] FixedIndices,
    double[] Loads,
    double[] FixedPositions,
    double[] QInit,
    double[] LowerBounds,
    double[] UpperBounds,
    int[] VariableNodeIndices,
    int[] VariableSupportKinds,
    double[] SphereRadii,
    byte[] RollerEnabled,
    double[] RollerLower,
    double[] RollerUpper,
    double[] RailStart,
    double[] RailEnd);
