namespace Ariadne.FEA;

/// <summary>
/// Cross-section properties for FEA elements.
/// For trusses: only Area is needed.
/// For beams (frames): Iy, Iz, J are required; Asy/Asz for Timoshenko only.
/// For shells: Offset shifts the mid-surface reference plane.
/// </summary>
public sealed class FeaSection
{
    public double Area { get; init; } = 0.01;
    public double Iy { get; init; } = 0.0;
    public double Iz { get; init; } = 0.0;
    public double J { get; init; } = 0.0;
    public double Asy { get; init; } = 0.0;
    public double Asz { get; init; } = 0.0;
    public double Offset { get; init; } = 0.0;
}
