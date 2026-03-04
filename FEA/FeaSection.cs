namespace Ariadne.FEA;

/// <summary>
/// Cross-section properties for FEA elements.
/// Phase 1: area only. Phase 2 adds Ix, Iy, J for beams.
/// </summary>
public sealed class FeaSection
{
    public double Area { get; init; } = 0.01;
}
