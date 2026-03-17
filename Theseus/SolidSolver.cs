using System;
using System.Runtime.InteropServices;

namespace Theseus.Interop;

/// <summary>
/// Number of Gauss points per tet4 element. Must match NUM_GP in the Rust core.
/// </summary>
public static class SolidConstants
{
    public const int NumGP = 4;
}

/// <summary>
/// Result of a solid FEA forward solve.
/// All arrays use row-major layout.
/// </summary>
public sealed class SolidSolverResult
{
    /// <summary>Flat displacement array, length num_nodes * 3.</summary>
    public double[] Displacements { get; init; } = [];

    /// <summary>Flat deformed node positions, length num_nodes * 3.</summary>
    public double[] DeformedXyz { get; init; } = [];

    /// <summary>Flat reaction forces at all nodes, length num_nodes * 3.</summary>
    public double[] Reactions { get; init; } = [];

    /// <summary>
    /// Per-GP stress tensors: [elem0_gp0(6), elem0_gp1(6), ..., elem0_gp3(6), elem1_gp0(6), ...].
    /// Length = num_elements * NumGP * 6.
    /// </summary>
    public double[] Stresses { get; init; } = [];

    /// <summary>
    /// Per-GP von Mises: [elem0_gp0, elem0_gp1, ..., elem0_gp3, elem1_gp0, ...].
    /// Length = num_elements * NumGP.
    /// </summary>
    public double[] VonMises { get; init; } = [];
}

/// <summary>
/// Raw P/Invoke declarations for the solid FEA functions in theseus.dll.
/// Maps 1:1 to <c>extern "C"</c> symbols in <c>rust/src/solid_ffi.rs</c>.
/// </summary>
internal static class SolidInterop
{
    const string DLL = "theseus";

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern IntPtr theseus_solid_create(
        nuint num_nodes, nuint num_elements,
        double[] node_positions, int[] elements,
        nuint num_materials, double[] materials,
        int[] element_props,
        nuint num_supports, int[] supports,
        nuint num_loads, double[] load_forces, int[] load_nodes,
        byte include_self_weight, double[] gravity);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_solid_solve_forward(
        IntPtr handle,
        double[] out_displacements, double[] out_deformed_xyz,
        double[] out_reactions, double[] out_stresses,
        double[] out_von_mises);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void theseus_solid_destroy(IntPtr handle);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_solid_last_error(byte[] buf, nuint buf_len);
}

/// <summary>
/// Managed wrapper around the native solid FEA solver in theseus.dll.
/// </summary>
public sealed class SolidSolver : IDisposable
{
    private IntPtr _handle;
    private readonly int _numNodes;
    private readonly int _numElements;
    private bool _disposed;

    private SolidSolver(IntPtr handle, int numNodes, int numElements)
    {
        _handle = handle;
        _numNodes = numNodes;
        _numElements = numElements;
    }

    ~SolidSolver() { Dispose(); }

    private void ThrowIfDisposed()
    {
        if (_disposed) throw new ObjectDisposedException(nameof(SolidSolver));
    }

    private static string GetLastError()
    {
        var buf = new byte[2048];
        int n = SolidInterop.theseus_solid_last_error(buf, (nuint)buf.Length);
        if (n <= 0) return string.Empty;
        return System.Text.Encoding.UTF8.GetString(buf, 0, n);
    }

    private static void Check(int rc)
    {
        if (rc != 0)
            throw new TheseusException(GetLastError(), rc);
    }

    // ── Construction ─────────────────────────────────────────

    /// <summary>
    /// Create a solid FEA solver from flat arrays.
    /// </summary>
    /// <param name="numNodes">Number of nodes in the mesh.</param>
    /// <param name="numElements">Number of tet4 elements.</param>
    /// <param name="nodePositions">Flat [x0,y0,z0, x1,...], length numNodes*3.</param>
    /// <param name="elements">Flat [n0,n1,n2,n3, ...], length numElements*4.</param>
    /// <param name="numMaterials">Number of materials.</param>
    /// <param name="materials">Flat [E,nu,density,yieldStress, ...], length numMaterials*4.</param>
    /// <param name="elementProps">Material index per element, length numElements.</param>
    /// <param name="numSupports">Number of support conditions.</param>
    /// <param name="supports">Flat [nodeIdx,fixX,fixY,fixZ, ...], length numSupports*4.</param>
    /// <param name="numLoads">Number of point loads.</param>
    /// <param name="loadForces">Flat [fx,fy,fz, ...], length numLoads*3.</param>
    /// <param name="loadNodes">Node index per load, length numLoads.</param>
    /// <param name="includeSelfWeight">Whether to include self-weight.</param>
    /// <param name="gravity">Gravity vector [gx,gy,gz], length 3.</param>
    public static SolidSolver Create(
        int numNodes, int numElements,
        double[] nodePositions, int[] elements,
        int numMaterials, double[] materials,
        int[] elementProps,
        int numSupports, int[] supports,
        int numLoads, double[] loadForces, int[] loadNodes,
        bool includeSelfWeight, double[] gravity)
    {
        var handle = SolidInterop.theseus_solid_create(
            (nuint)numNodes, (nuint)numElements,
            nodePositions, elements,
            (nuint)numMaterials, materials,
            elementProps,
            (nuint)numSupports, supports,
            (nuint)numLoads, loadForces, loadNodes,
            (byte)(includeSelfWeight ? 1 : 0), gravity);

        if (handle == IntPtr.Zero)
            throw new TheseusException(GetLastError(), -1);

        return new SolidSolver(handle, numNodes, numElements);
    }

    // ── Forward solve ────────────────────────────────────────

    public SolidSolverResult SolveForward()
    {
        ThrowIfDisposed();
        int nn = _numNodes;
        int ne = _numElements;

        int ngp = SolidConstants.NumGP;
        var displacements = new double[nn * 3];
        var deformedXyz = new double[nn * 3];
        var reactions = new double[nn * 3];
        var stresses = new double[ne * ngp * 6];
        var vonMises = new double[ne * ngp];

        Check(SolidInterop.theseus_solid_solve_forward(
            _handle, displacements, deformedXyz,
            reactions, stresses, vonMises));

        return new SolidSolverResult
        {
            Displacements = displacements,
            DeformedXyz = deformedXyz,
            Reactions = reactions,
            Stresses = stresses,
            VonMises = vonMises,
        };
    }

    // ── Disposal ─────────────────────────────────────────────

    public void Dispose()
    {
        if (_disposed) return;
        _disposed = true;
        if (_handle != IntPtr.Zero)
        {
            SolidInterop.theseus_solid_destroy(_handle);
            _handle = IntPtr.Zero;
        }
        GC.SuppressFinalize(this);
    }
}
