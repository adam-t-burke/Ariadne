namespace Ariadne.FEA;

/// <summary>
/// Material properties for FEA elements (bar and solid).
/// </summary>
public sealed class FeaMaterial
{
    public double E { get; init; }
    public double Nu { get; init; } = 0.3;
    public double Density { get; init; }
    public double YieldStress { get; init; }

    public static FeaMaterial Steel() => new()
    {
        E = 210e9,
        Nu = 0.3,
        Density = 7850.0,
        YieldStress = 250e6,
    };

    public static FeaMaterial Concrete() => new()
    {
        E = 30e9,
        Nu = 0.2,
        Density = 2400.0,
        YieldStress = 30e6,
    };

    public static FeaMaterial Aluminum() => new()
    {
        E = 70e9,
        Nu = 0.33,
        Density = 2700.0,
        YieldStress = 270e6,
    };

    public static FeaMaterial Timber() => new()
    {
        E = 12e9,
        Nu = 0.35,
        Density = 500.0,
        YieldStress = 40e6,
    };
}
