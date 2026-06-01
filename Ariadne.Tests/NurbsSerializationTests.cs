using System.Collections.Generic;
using Ariadne.Utilities;
using Rhino.Geometry;
using Xunit;

namespace Ariadne.Tests;

public class NurbsSerializationTests
{
    [Fact]
    public void FromCurve_ExportsFullBookKnotVectorAndHomogeneousPoints()
    {
        var curve = NurbsCurve.Create(
            false,
            2,
            new List<Point3d>
            {
                new(0, 0, 0),
                new(1, 1, 0),
                new(2, 0, 0),
            });

        var data = NurbsSerialization.FromNurbsCurve(curve);

        Assert.Equal(curve.Points.Count + curve.Degree + 1, data.Knots.Length);
        Assert.Equal(curve.Points.Count * 4, data.ControlPoints.Length);
        Assert.Equal(curve.Domain.T0, data.Domain[0]);
        Assert.Equal(curve.Domain.T1, data.Domain[1]);
    }

    [Fact]
    public void FromSurface_ExportsControlPointGrid()
    {
        using var surface = NurbsSurface.CreateFromCorners(
            new Point3d(0, 0, 0),
            new Point3d(1, 0, 0),
            new Point3d(1, 1, 0),
            new Point3d(0, 1, 0));

        var data = NurbsSerialization.FromNurbsSurface(surface);

        Assert.Equal(surface.Points.CountU, data.CountU);
        Assert.Equal(surface.Points.CountV, data.CountV);
        Assert.Equal(data.CountU * data.CountV * 4, data.ControlPoints.Length);
        Assert.Equal(surface.Points.CountU + surface.Degree(0) + 1, data.KnotsU.Length);
        Assert.Equal(surface.Points.CountV + surface.Degree(1) + 1, data.KnotsV.Length);
    }
}
