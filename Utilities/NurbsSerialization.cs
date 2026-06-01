using System;
using System.Collections.Generic;
using Rhino.Geometry;

namespace Ariadne.Utilities;

/// <summary>
/// Converts Rhino NURBS geometry to the flat, book-style knot/control-point
/// layout consumed by the native nurbsbook crate.
/// </summary>
public static class NurbsSerialization
{
    public sealed record CurveData(
        int Degree,
        double[] Knots,
        double[] ControlPoints,
        double[] Domain);

    public sealed record SurfaceData(
        int DegreeU,
        int DegreeV,
        int CountU,
        int CountV,
        double[] KnotsU,
        double[] KnotsV,
        double[] ControlPoints,
        double[] DomainU,
        double[] DomainV);

    public static CurveData FromCurve(Curve curve)
    {
        if (curve == null)
            throw new ArgumentNullException(nameof(curve));

        using var nurbs = curve.ToNurbsCurve();
        if (nurbs == null || !nurbs.IsValid)
            throw new ArgumentException("Curve could not be converted to a valid NURBS curve.", nameof(curve));

        return FromNurbsCurve(nurbs);
    }

    public static CurveData FromNurbsCurve(NurbsCurve curve)
    {
        if (curve == null)
            throw new ArgumentNullException(nameof(curve));
        if (!curve.IsValid)
            throw new ArgumentException("NURBS curve is invalid.", nameof(curve));

        int count = curve.Points.Count;
        var cp = new double[count * 4];
        for (int i = 0; i < count; i++)
        {
            var p = curve.Points[i];
            double w = p.Weight;
            cp[i * 4 + 0] = p.Location.X * w;
            cp[i * 4 + 1] = p.Location.Y * w;
            cp[i * 4 + 2] = p.Location.Z * w;
            cp[i * 4 + 3] = w;
        }

        return new CurveData(
            curve.Degree,
            ExpandRhinoKnots(curve.Knots),
            cp,
            [curve.Domain.T0, curve.Domain.T1]);
    }

    public static SurfaceData FromSurface(Surface surface)
    {
        if (surface == null)
            throw new ArgumentNullException(nameof(surface));

        using var nurbs = surface.ToNurbsSurface();
        if (nurbs == null || !nurbs.IsValid)
            throw new ArgumentException("Surface could not be converted to a valid NURBS surface.", nameof(surface));

        return FromNurbsSurface(nurbs);
    }

    public static SurfaceData FromNurbsSurface(NurbsSurface surface)
    {
        if (surface == null)
            throw new ArgumentNullException(nameof(surface));
        if (!surface.IsValid)
            throw new ArgumentException("NURBS surface is invalid.", nameof(surface));

        int countU = surface.Points.CountU;
        int countV = surface.Points.CountV;
        var cp = new double[countU * countV * 4];
        for (int u = 0; u < countU; u++)
        {
            for (int v = 0; v < countV; v++)
            {
                var p = surface.Points.GetControlPoint(u, v);
                double w = p.Weight;
                int idx = (u * countV + v) * 4;
                cp[idx + 0] = p.Location.X * w;
                cp[idx + 1] = p.Location.Y * w;
                cp[idx + 2] = p.Location.Z * w;
                cp[idx + 3] = w;
            }
        }

        var domainU = surface.Domain(0);
        var domainV = surface.Domain(1);
        return new SurfaceData(
            surface.Degree(0),
            surface.Degree(1),
            countU,
            countV,
            ExpandRhinoKnots(surface.KnotsU),
            ExpandRhinoKnots(surface.KnotsV),
            cp,
            [domainU.T0, domainU.T1],
            [domainV.T0, domainV.T1]);
    }

    public static void AddToHash(HashCode hash, CurveData data)
    {
        hash.Add(data.Degree);
        AddRange(hash, data.Knots);
        AddRange(hash, data.ControlPoints);
        AddRange(hash, data.Domain);
    }

    public static void AddToHash(HashCode hash, SurfaceData data)
    {
        hash.Add(data.DegreeU);
        hash.Add(data.DegreeV);
        hash.Add(data.CountU);
        hash.Add(data.CountV);
        AddRange(hash, data.KnotsU);
        AddRange(hash, data.KnotsV);
        AddRange(hash, data.ControlPoints);
        AddRange(hash, data.DomainU);
        AddRange(hash, data.DomainV);
    }

    private static void AddRange(HashCode hash, IReadOnlyList<double> values)
    {
        for (int i = 0; i < values.Count; i++)
            hash.Add(values[i]);
    }

    /// <summary>
    /// Rhino omits the two superfluous endpoint knots. The NURBS Book algorithms
    /// use the full vector, so duplicate the first and last Rhino knots.
    /// </summary>
    private static double[] ExpandRhinoKnots(dynamic knots)
    {
        if (knots.Count == 0)
            throw new ArgumentException("Knot list is empty.", nameof(knots));

        var result = new double[knots.Count + 2];
        result[0] = knots[0];
        for (int i = 0; i < knots.Count; i++)
            result[i + 1] = knots[i];
        result[^1] = knots[knots.Count - 1];
        return result;
    }
}
