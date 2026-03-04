namespace Ariadne.FEA;

using System.Collections.Generic;
using Rhino.Geometry;

/// <summary>
/// FEA optimization configuration.
/// </summary>
public sealed record FeaOptimizationConfig
{
    public required IReadOnlyList<FeaObjective> Objectives { get; init; }
    public bool OptimizeNodePositions { get; init; } = true;
    public List<int>? NodeIndices { get; init; }
    public bool OptimizeAreas { get; init; } = false;
    public List<int>? AreaIndices { get; init; }
    public bool OptimizeSupportPositions { get; init; } = false;
    public List<int>? SupportNodeIndices { get; init; }
    public int MaxIterations { get; init; } = 500;
    public double AbsTol { get; init; } = 1e-6;
    public double RelTol { get; init; } = 1e-6;
    public double BarrierWeight { get; init; } = 10.0;
    public double BarrierSharpness { get; init; } = 10.0;
    public int ReportFrequency { get; init; } = 10;
    public bool Run { get; init; } = false;
    public bool StreamPreview { get; init; } = true;
}

/// <summary>
/// Results from an FEA solve.
/// </summary>
public sealed record FeaSolveResult
{
    public required FEA_Network Network { get; init; }
    public required double[] Displacements { get; init; }
    public required double[] Reactions { get; init; }
    public required double[] AxialForces { get; init; }
    public required double[] Stresses { get; init; }
    public required double[] Strains { get; init; }
    public required double[] Utilization { get; init; }
    public required Point3d[] DeformedNodes { get; init; }
    public required LineCurve[] DeformedEdges { get; init; }
    public required int Iterations { get; init; }
    public required bool Converged { get; init; }
}

/// <summary>
/// Internal data for FEA solver creation.
/// </summary>
internal sealed record FeaSolverData(
    int NumNodes,
    int NumElements,
    double[] NodePositions,
    int[] EdgeNodes,
    int NumMaterials,
    double[] Materials,
    int NumSections,
    double[] Sections,
    int[] ElementProps,
    int NumSupports,
    int[] Supports,
    int NumLoads,
    double[] LoadForces,
    int[] LoadNodes,
    bool IncludeSelfWeight,
    double[] Gravity);
