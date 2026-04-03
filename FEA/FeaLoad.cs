namespace Ariadne.FEA;

using Rhino.Geometry;

/// <summary>
/// Canonical nodal load definition for a unified FEA model.
/// NodeIndex references the global node numbering in <see cref="FEA_Model"/>.
/// </summary>
public sealed record FeaLoad
{
    public required int NodeIndex { get; init; }
    public required Vector3d Force { get; init; }
    public Vector3d Moment { get; init; } = Vector3d.Zero;
}
