namespace Ariadne.FEA;

/// <summary>
/// Material properties for FEA elements.
/// </summary>
public sealed class FeaMaterial
{
    public double E { get; init; }
    public double Density { get; init; }
    public double YieldStress { get; init; }

    public static FeaMaterial Steel() => new()
    {
        E = 210e9,
        Density = 7850.0,
        YieldStress = 250e6,
    };
}
