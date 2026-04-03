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

internal sealed record BarSolveData
{
    public required Vector3d[] Displacements { get; init; }
    public required Vector3d[] Reactions { get; init; }
    public required int[] ReactionNodeIndices { get; init; }
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

internal sealed record ShellSolveData
{
    public required Vector3d[] Displacements { get; init; }
    public required Vector3d[] Rotations { get; init; }
    public required Vector3d[] ReactionForces { get; init; }
    public required Vector3d[] ReactionMoments { get; init; }
    public required Point3d[] DeformedNodes { get; init; }
    public required double[] S1 { get; init; }
    public required double[] S2 { get; init; }
    public required double[] TopStresses { get; init; }
    public required double[] BottomStresses { get; init; }
    public required double[] VonMises { get; init; }
}

internal sealed record ShellSolverData(
    int NumNodes,
    int NumElements,
    double[] NodePositions,
    double[] NodeThicknesses,
    int[] Elements,
    int[] NumNodesPerElement,
    int NumMaterials,
    double[] Materials,
    int NumSections,
    double[] Sections,
    int[] ElementProps,
    int NumSupports,
    int[] Supports,
    int NumLoads,
    double[] LoadForces,
    double[] LoadMoments,
    int[] LoadNodes);
