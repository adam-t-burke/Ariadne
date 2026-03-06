namespace Ariadne.FEA;

using Rhino.Geometry;

/// <summary>
/// Canonical nodal load definition for a unified FEA model.
/// Node indices use the same global numbering as <see cref="FEA_Model"/>.
/// </summary>
public sealed record FeaLoad
{
    public required int NodeIndex { get; init; }
    public required Vector3d Force { get; init; }
}
