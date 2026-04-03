namespace Ariadne.FEA;

using System.Collections.Generic;
using System.Linq;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

/// <summary>
/// Pure geometry + topology container for a tetrahedral volume mesh.
/// Inherits from <see cref="GH_Mesh"/> so the boundary surface mesh is
/// available via the inherited Value property for Grasshopper preview/piping.
/// The tet4 element connectivity and full node list live as additional properties.
/// </summary>
public class SolidMesh : GH_Mesh
{
    /// <summary>Full nodal coordinates (all tet nodes, not just boundary).</summary>
    public List<Point3d> TetNodes { get; set; } = [];

    /// <summary>Tet4 element connectivity. Each int[4] contains the four node indices.</summary>
    public List<int[]> Elements { get; set; } = [];

    /// <summary>Boundary face connectivity. Each int[3] contains the three node indices of a boundary triangle.</summary>
    public List<int[]> BoundaryFaces { get; set; } = [];

    public FeaMaterial? Material { get; set; }

    private bool _boundaryMeshDirty = true;

    public int Nn => TetNodes.Count;
    public int Ne => Elements.Count;

    public SolidMesh() : base() { }

    public SolidMesh(List<Point3d> nodes, List<int[]> elements, List<int[]> boundaryFaces)
        : base()
    {
        TetNodes = nodes;
        Elements = elements;
        BoundaryFaces = boundaryFaces;
        _boundaryMeshDirty = true;
    }

    public SolidMesh(SolidMesh other) : base()
    {
        TetNodes = new List<Point3d>(other.TetNodes);
        Elements = other.Elements.Select(e => (int[])e.Clone()).ToList();
        BoundaryFaces = other.BoundaryFaces.Select(f => (int[])f.Clone()).ToList();
        _boundaryMeshDirty = true;
    }

    public override Mesh Value
    {
        get
        {
            if (_boundaryMeshDirty)
            {
                base.Value = BuildBoundaryMesh();
                _boundaryMeshDirty = false;
            }
            return base.Value;
        }
        set
        {
            base.Value = value;
            _boundaryMeshDirty = false;
        }
    }

    public override bool IsValid => Nn > 0 && Ne > 0;

    public override string ToString() => $"SolidMesh ({Nn} nodes, {Ne} elements)";

    /// <summary>
    /// Enables SolidMesh to convert to Rhino.Geometry.Mesh for native Grasshopper components (e.g. Deconstruct Mesh).
    /// </summary>
    public override bool CastTo<T>(out T target)
    {
        if (typeof(T) == typeof(Mesh))
        {
            target = (T)(object)Value;
            return true;
        }
        target = default!;
        return false;
    }

    /// <summary>
    /// Compute the volume of a single tet4 element.
    /// V = (1/6) |det([x1-x0, x2-x0, x3-x0])|
    /// </summary>
    public double ElementVolume(int elementIndex)
    {
        var e = Elements[elementIndex];
        var p0 = TetNodes[e[0]];
        var p1 = TetNodes[e[1]];
        var p2 = TetNodes[e[2]];
        var p3 = TetNodes[e[3]];

        var a = p1 - p0;
        var b = p2 - p0;
        var c = p3 - p0;

        double det = a.X * (b.Y * c.Z - b.Z * c.Y)
                   - a.Y * (b.X * c.Z - b.Z * c.X)
                   + a.Z * (b.X * c.Y - b.Y * c.X);

        return System.Math.Abs(det) / 6.0;
    }

    /// <summary>Compute the centroid of a single tet4 element.</summary>
    public Point3d ElementCentroid(int elementIndex)
    {
        var e = Elements[elementIndex];
        return new Point3d(
            (TetNodes[e[0]].X + TetNodes[e[1]].X + TetNodes[e[2]].X + TetNodes[e[3]].X) * 0.25,
            (TetNodes[e[0]].Y + TetNodes[e[1]].Y + TetNodes[e[2]].Y + TetNodes[e[3]].Y) * 0.25,
            (TetNodes[e[0]].Z + TetNodes[e[1]].Z + TetNodes[e[2]].Z + TetNodes[e[3]].Z) * 0.25);
    }

    public List<int> GetBoundaryFaceNodeIndices(IEnumerable<int> faceIndices)
    {
        var nodes = new HashSet<int>();
        foreach (int faceIndex in faceIndices)
        {
            if (faceIndex < 0 || faceIndex >= BoundaryFaces.Count)
                continue;

            foreach (int nodeIndex in BoundaryFaces[faceIndex])
                nodes.Add(nodeIndex);
        }

        return nodes.OrderBy(i => i).ToList();
    }

    public bool TryGetBoundaryFaceLoadData(int boundaryFaceIndex, out Vector3d outwardNormal, out double area)
    {
        outwardNormal = Vector3d.Unset;
        area = 0.0;

        if (boundaryFaceIndex < 0 || boundaryFaceIndex >= BoundaryFaces.Count)
            return false;

        var face = BoundaryFaces[boundaryFaceIndex];
        var p0 = TetNodes[face[0]];
        var p1 = TetNodes[face[1]];
        var p2 = TetNodes[face[2]];

        var rawNormal = Vector3d.CrossProduct(p1 - p0, p2 - p0);
        area = rawNormal.Length * 0.5;
        if (area < 1e-30)
            return false;

        rawNormal.Unitize();

        int oppositeNodeIndex = -1;
        foreach (var element in Elements)
        {
            if (!element.Contains(face[0]) || !element.Contains(face[1]) || !element.Contains(face[2]))
                continue;

            oppositeNodeIndex = element.First(nodeIndex => nodeIndex != face[0] && nodeIndex != face[1] && nodeIndex != face[2]);
            break;
        }

        if (oppositeNodeIndex >= 0)
        {
            var faceCentroid = new Point3d(
                (p0.X + p1.X + p2.X) / 3.0,
                (p0.Y + p1.Y + p2.Y) / 3.0,
                (p0.Z + p1.Z + p2.Z) / 3.0);
            var toInterior = TetNodes[oppositeNodeIndex] - faceCentroid;
            if (Vector3d.Multiply(rawNormal, toInterior) > 0.0)
                rawNormal = -rawNormal;
        }

        outwardNormal = rawNormal;
        return true;
    }

    private Mesh BuildBoundaryMesh()
    {
        var rhinoMesh = new Mesh();
        if (TetNodes.Count == 0 || BoundaryFaces.Count == 0)
            return rhinoMesh;

        rhinoMesh.Vertices.AddVertices(TetNodes);
        foreach (var face in BoundaryFaces)
            rhinoMesh.Faces.AddFace(face[0], face[1], face[2]);

        rhinoMesh.Normals.ComputeNormals();
        rhinoMesh.Compact();
        return rhinoMesh;
    }
}
