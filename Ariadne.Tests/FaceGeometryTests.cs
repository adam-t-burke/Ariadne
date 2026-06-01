using System.Collections.Generic;
using Ariadne.Graphs;
using Rhino.Geometry;
using Xunit;

namespace Ariadne.Tests;

public class FaceGeometryTests
{
    [Fact]
    public void ComputeNewellNormal_UnitSquareInXY_IsPositiveZ()
    {
        var verts = new List<Point3d>
        {
            new(0, 0, 0),
            new(1, 0, 0),
            new(1, 1, 0),
            new(0, 1, 0),
        };

        var normal = FaceGeometry.ComputeNewellNormal(verts);
        Assert.True(normal.Z > 0.9);
    }

    [Fact]
    public void ValidateFaceIndices_RejectsOutOfRange()
    {
        var error = FaceGeometry.ValidateFaceIndices([0, 1, 99], faceIndex: 2, nodeCount: 4);
        Assert.Contains("out of range", error);
    }
}
