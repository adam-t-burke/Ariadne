namespace Ariadne.FEA;

using Rhino.Geometry;

/// <summary>
/// Support boundary condition for FEA nodes.
/// Per-DOF constraint: fix any combination of X, Y, Z translations.
/// Phase 2 adds FixRx, FixRy, FixRz for beam elements.
/// </summary>
public sealed class FeaSupport
{
    public Point3d Location { get; init; }
    public int? NodeIndex { get; init; }
    public bool FixX { get; init; } = true;
    public bool FixY { get; init; } = true;
    public bool FixZ { get; init; } = true;
    public bool FixRX { get; init; } = false; // 6-DOF support
    public bool FixRY { get; init; } = false;
    public bool FixRZ { get; init; } = false;

    public static FeaSupport Fixed(Point3d location) => new()
    {
        Location = location,
        FixX = true,
        FixY = true,
        FixZ = true,
        FixRX = true,
        FixRY = true,
        FixRZ = true,
    };

    public static FeaSupport Pinned(Point3d location) => new()
    {
        Location = location,
        FixX = true,
        FixY = true,
        FixZ = true,
    };

    public static FeaSupport RollerX(Point3d location) => new()
    {
        Location = location,
        FixX = true,
        FixY = false,
        FixZ = false,
    };

    public static FeaSupport RollerY(Point3d location) => new()
    {
        Location = location,
        FixX = false,
        FixY = true,
        FixZ = false,
    };

    public static FeaSupport RollerZ(Point3d location) => new()
    {
        Location = location,
        FixX = false,
        FixY = false,
        FixZ = true,
    };
}
