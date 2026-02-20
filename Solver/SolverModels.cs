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
    public int MaxIterations { get; init; } = 500;
    public double AbsTol { get; init; } = 1e-6;
    public double RelTol { get; init; } = 1e-6;
    public double BarrierWeight { get; init; } = 1000.0;
    public double BarrierSharpness { get; init; } = 10.0;
}

/// <summary>
/// Inputs required for the solver.
/// </summary>
public sealed record SolverInputs
{
    public required List<double> QInit { get; init; }
    public required List<Vector3d> Loads { get; init; }
    public required List<double> LowerBounds { get; init; }
    public required List<double> UpperBounds { get; init; }
    public List<Objective> Objectives { get; init; } = [];
}

/// <summary>
/// Results from a solve operation.
/// </summary>
public sealed record SolveResult
{
    public required FDM_Network Network { get; init; }
    public required double[] ForceDensities { get; init; }
    public required double[] MemberForces { get; init; }
    public required double[] MemberLengths { get; init; }
    public required double[] Reactions { get; init; }
    public required int Iterations { get; init; }
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
