namespace Ariadne.FEA;

using System.Collections.Generic;
using System.Linq;
using Rhino.Geometry;

/// <summary>
/// Geometry + topology container for a triangulated shell mesh, analogous to
/// <see cref="SolidMesh"/> for tets. Stores pre-triangulated elements, per-node
/// thickness, and boundary edges so downstream components can discover nodes
/// before FEA model construction.
/// </summary>
public class ShellMesh
{
    public List<Point3d> Nodes { get; }
    public List<int[]> Elements { get; }
    public double[] NodeThicknesses { get; }
    public List<(int, int)> BoundaryEdges { get; }
    public FeaMaterial Material { get; }
    public Mesh OriginalMesh { get; }
    public int OriginalFaceCount { get; }
    public int[] TriToFaceMap { get; }

    public int Nn => Nodes.Count;
    public int Ne => Elements.Count;

    public ShellMesh(Mesh mesh, double[] thicknesses, FeaMaterial material)
    {
        OriginalMesh = mesh;
        OriginalFaceCount = mesh.Faces.Count;
        Material = material;

        Nodes = mesh.Vertices.Select(v => (Point3d)v).ToList();
        NodeThicknesses = thicknesses;

        var elements = new List<int[]>();
        var triToFace = new List<int>();

        for (int i = 0; i < mesh.Faces.Count; i++)
        {
            var f = mesh.Faces[i];
            elements.Add([f.A, f.B, f.C]);
            triToFace.Add(i);

            if (f.IsQuad)
            {
                elements.Add([f.A, f.C, f.D]);
                triToFace.Add(i);
            }
        }

        Elements = elements;
        TriToFaceMap = triToFace.ToArray();

        BoundaryEdges = ComputeBoundaryEdges();
    }

    private List<(int, int)> ComputeBoundaryEdges()
    {
        var edgeCounts = new Dictionary<(int, int), int>();

        foreach (var tri in Elements)
        {
            for (int i = 0; i < 3; i++)
            {
                int a = tri[i];
                int b = tri[(i + 1) % 3];
                var edge = a < b ? (a, b) : (b, a);
                edgeCounts.TryGetValue(edge, out int count);
                edgeCounts[edge] = count + 1;
            }
        }

        return edgeCounts
            .Where(kvp => kvp.Value == 1)
            .Select(kvp => kvp.Key)
            .OrderBy(e => e.Item1)
            .ThenBy(e => e.Item2)
            .ToList();
    }

    /// <summary>
    /// Build a Rhino Mesh from the triangulated elements for visualization.
    /// </summary>
    public Mesh BuildTriangulatedMesh()
    {
        var rhinoMesh = new Mesh();
        rhinoMesh.Vertices.AddVertices(Nodes);
        foreach (var tri in Elements)
            rhinoMesh.Faces.AddFace(tri[0], tri[1], tri[2]);
        rhinoMesh.Normals.ComputeNormals();
        return rhinoMesh;
    }

    public override string ToString() => $"ShellMesh ({Nn} nodes, {Ne} triangles, {OriginalFaceCount} original faces)";
}
