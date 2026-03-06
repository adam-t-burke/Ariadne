namespace Ariadne.FEA;

using System.Collections.Generic;

/// <summary>
/// Exact support region for solid meshes, defined on explicit solid node sets.
/// </summary>
public sealed record FeaSolidSupportRegion
{
    public required List<int> NodeIndices { get; init; }
    public List<int> FaceIndices { get; init; } = [];
    public bool FixX { get; init; } = true;
    public bool FixY { get; init; } = true;
    public bool FixZ { get; init; } = true;
}
