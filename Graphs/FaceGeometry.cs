using System.Collections.Generic;
using Rhino.Geometry;

namespace Ariadne.Graphs;

/// <summary>
/// Shared face geometry helpers used by face-finding and pressure-load components.
/// </summary>
internal static class FaceGeometry
{
    /// <summary>
    /// Computes a unit face normal via Newell's method from an ordered vertex ring.
    /// </summary>
    public static Vector3d ComputeNewellNormal(IReadOnlyList<Point3d> vertices)
    {
        if (vertices.Count < 3)
            return Vector3d.Zero;

        var normal = Vector3d.Zero;
        for (int i = 0; i < vertices.Count; i++)
        {
            var vi = vertices[i];
            var vj = vertices[(i + 1) % vertices.Count];
            normal.X += (vi.Y - vj.Y) * (vi.Z + vj.Z);
            normal.Y += (vi.Z - vj.Z) * (vi.X + vj.X);
            normal.Z += (vi.X - vj.X) * (vi.Y + vj.Y);
        }

        normal *= 0.5;
        if (normal.Length > 1e-12)
            normal.Unitize();
        return normal;
    }

    /// <summary>
    /// Builds a closed polyline curve for a face vertex ring.
    /// </summary>
    public static PolylineCurve BuildFacePolyline(IReadOnlyList<Point3d> vertices)
    {
        var pts = new List<Point3d>(vertices.Count + 1);
        foreach (var pt in vertices)
            pts.Add(pt);
        pts.Add(pts[0]);
        return new PolylineCurve(new Polyline(pts));
    }

    /// <summary>
    /// Validates face vertex indices against a graph node count.
    /// Returns null when valid; otherwise an error message.
    /// </summary>
    public static string? ValidateFaceIndices(IReadOnlyList<int> face, int faceIndex, int nodeCount)
    {
        if (face.Count < 3)
            return $"Face {faceIndex}: fewer than 3 vertices.";

        for (int v = 0; v < face.Count; v++)
        {
            int idx = face[v];
            if (idx < 0 || idx >= nodeCount)
                return $"Face {faceIndex}: vertex index {idx} out of range [0, {nodeCount - 1}].";
        }

        return null;
    }
}
