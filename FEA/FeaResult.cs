namespace Ariadne.FEA;

using Rhino.Geometry;

public sealed record FeaResult
{
    public required FEA_Model Model { get; init; }
    public required FEA_Model DeformedModel { get; init; }
    public required Vector3d[] Displacements { get; init; }
    public required Vector3d[] Reactions { get; init; }
    public required int[] ReactionNodeIndices { get; init; }
    public required Point3d[] DeformedNodes { get; init; }
    public double[]? AxialForces { get; init; }
    public double[]? BarStresses { get; init; }
    public double[]? BarStrains { get; init; }
    public double[,]? SolidStresses { get; init; }
    public double[]? VonMises { get; init; }
    public Vector3d[]? ShellRotations { get; init; }
    public Vector3d[]? ShellReactionMoments { get; init; }
    public double[]? ShellS1 { get; init; }
    public double[]? ShellS2 { get; init; }
    public double[]? ShellVonMises { get; init; }
    public double[]? ShellTopStresses { get; init; }
    public double[]? ShellBottomStresses { get; init; }
    public required double[] Utilization { get; init; }
    public required int Iterations { get; init; }
    public required bool Converged { get; init; }
}
