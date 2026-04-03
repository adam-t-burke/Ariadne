using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Text;
using Rhino.Geometry;

namespace Ariadne.TetGen;

internal static class TetGenNative
{
    const string DLL = "tetgen_wrapper";

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int tetgen_mesh(
        int num_points, double[] points,
        int num_facets, int[] facets,
        double max_volume, double min_ratio,
        int quadratic,
        out int out_num_points, out IntPtr out_points,
        out int out_num_tets, out IntPtr out_tets,
        out int out_nodes_per_tet,
        out int out_num_faces, out IntPtr out_faces);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void tetgen_free(IntPtr ptr);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int tetgen_last_error(byte[] buf, int buf_len);
}

/// <summary>
/// Managed wrapper around the native TetGen library.
/// Converts Rhino Mesh to/from flat arrays for the C FFI.
/// </summary>
public static class TetGenMesher
{
    public static string GetLastError()
    {
        var buf = new byte[2048];
        int n = TetGenNative.tetgen_last_error(buf, buf.Length);
        if (n <= 0) return string.Empty;
        return Encoding.UTF8.GetString(buf, 0, n);
    }

    /// <summary>
    /// Generate a tetrahedral volume mesh from a closed surface mesh.
    /// </summary>
    /// <param name="mesh">Closed triangulated surface mesh.</param>
    /// <param name="maxVolume">Maximum tet volume (0 = no constraint).</param>
    /// <param name="minRatio">Minimum radius-edge quality ratio (0 = no quality, typical 1.2-2.0).</param>
    /// <param name="nodes">Output node positions.</param>
    /// <param name="elements">Output tet4 connectivity (zero-based).</param>
    /// <param name="boundaryFaces">Output boundary triangle connectivity (zero-based).</param>
    /// <returns>Null on success, error message on failure.</returns>
    public static string? Tetrahedralize(
        Mesh mesh,
        double maxVolume,
        double minRatio,
        bool quadratic,
        out List<Point3d> nodes,
        out List<int[]> elements,
        out List<int[]> boundaryFaces)
    {
        nodes = [];
        elements = [];
        boundaryFaces = [];

        mesh.Faces.ConvertQuadsToTriangles();

        int nv = mesh.Vertices.Count;
        int nf = mesh.Faces.Count;

        var points = new double[nv * 3];
        for (int i = 0; i < nv; i++)
        {
            var v = mesh.Vertices[i];
            points[i * 3 + 0] = v.X;
            points[i * 3 + 1] = v.Y;
            points[i * 3 + 2] = v.Z;
        }

        var facets = new int[nf * 3];
        for (int i = 0; i < nf; i++)
        {
            var f = mesh.Faces[i];
            facets[i * 3 + 0] = f.A;
            facets[i * 3 + 1] = f.B;
            facets[i * 3 + 2] = f.C;
        }

        int rc = TetGenNative.tetgen_mesh(
            nv, points,
            nf, facets,
            maxVolume, minRatio,
            quadratic ? 1 : 0,
            out int outNp, out IntPtr outPointsPtr,
            out int outNt, out IntPtr outTetsPtr,
            out int outNodesPerTet,
            out int outNf, out IntPtr outFacesPtr);

        if (rc != 0)
            return GetLastError();

        try
        {
            var outPoints = new double[outNp * 3];
            Marshal.Copy(outPointsPtr, outPoints, 0, outNp * 3);

            var outTets = new int[outNt * outNodesPerTet];
            Marshal.Copy(outTetsPtr, outTets, 0, outNt * outNodesPerTet);

            var outFaces = new int[outNf * 3];
            Marshal.Copy(outFacesPtr, outFaces, 0, outNf * 3);

            nodes = new List<Point3d>(outNp);
            for (int i = 0; i < outNp; i++)
                nodes.Add(new Point3d(outPoints[i * 3], outPoints[i * 3 + 1], outPoints[i * 3 + 2]));

            elements = new List<int[]>(outNt);
            for (int i = 0; i < outNt; i++)
            {
                var el = new int[outNodesPerTet];
                for (int j = 0; j < outNodesPerTet; j++)
                    el[j] = outTets[i * outNodesPerTet + j];
                elements.Add(el);
            }

            boundaryFaces = new List<int[]>(outNf);
            for (int i = 0; i < outNf; i++)
                boundaryFaces.Add(new int[] { outFaces[i * 3], outFaces[i * 3 + 1], outFaces[i * 3 + 2] });

            return null;
        }
        finally
        {
            TetGenNative.tetgen_free(outPointsPtr);
            TetGenNative.tetgen_free(outTetsPtr);
            TetGenNative.tetgen_free(outFacesPtr);
        }
    }
}
