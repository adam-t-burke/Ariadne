namespace Ariadne.Solver;

using System.Collections.Generic;
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
    /// <summary>Invoke progress callback every N evaluations (0 = every evaluation).</summary>
    public int ReportFrequency { get; init; } = 10;
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
    /// <summary>Progress callback frequency (evaluations between callbacks).</summary>
    public int ReportFrequency { get; init; } = 10;
    /// <summary>When true, optimization runs (e.g. from a button or toggle).</summary>
    public bool Run { get; init; } = false;
    /// <summary>When true, stream intermediate results to outputs during optimization.</summary>
    public bool StreamPreview { get; init; } = true;
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
    double[] UpperBounds);
