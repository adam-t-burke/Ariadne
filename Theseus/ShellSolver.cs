using System;
using System.Runtime.InteropServices;

namespace Theseus.Interop;

/// <summary>
/// Result of a shell FEA forward solve.
/// </summary>
public sealed class ShellSolverResult
{
    /// <summary>Flat displacement array (translation + rotation), length num_nodes * 6.</summary>
    public double[] Displacements { get; init; } = [];

    /// <summary>Flat reaction forces + moments at all nodes, length num_nodes * 6.</summary>
    public double[] Reactions { get; init; } = [];

    /// <summary>Flat principal stress array [S1, S2] per element, length num_elements * 2.</summary>
    public double[] PrincipalStresses { get; init; } = [];

    /// <summary>Per-element top-fiber global stress tensors, length num_elements * 6.</summary>
    public double[] TopStresses { get; init; } = [];

    /// <summary>Per-element bottom-fiber global stress tensors, length num_elements * 6.</summary>
    public double[] BottomStresses { get; init; } = [];

    /// <summary>Per-element von Mises (max of top/bottom), length num_elements.</summary>
    public double[] VonMises { get; init; } = [];
}

/// <summary>
/// Managed wrapper around the native shell FEA solver in theseus.dll.
/// </summary>
public sealed class ShellSolver : IDisposable
{
    private IntPtr _handle;
    private readonly int _numNodes;
    private readonly int _numElements;
    private bool _disposed;

    private ShellSolver(IntPtr handle, int numNodes, int numElements)
    {
        _handle = handle;
        _numNodes = numNodes;
        _numElements = numElements;
    }

    ~ShellSolver() { Dispose(); }

    private void ThrowIfDisposed()
    {
        if (_disposed) throw new ObjectDisposedException(nameof(ShellSolver));
    }

    private static string GetLastError()
    {
        var buf = new byte[2048];
        int n = ShellInterop.theseus_shell_last_error(buf, (nuint)buf.Length);
        if (n <= 0) return string.Empty;
        return System.Text.Encoding.UTF8.GetString(buf, 0, n);
    }

    private static void Check(int rc)
    {
        if (rc != 0)
            throw new Exception($"Shell solver error: {rc}. {GetLastError()}");
    }

    public static ShellSolver Create(
        int numNodes, int numElements,
        double[] nodePositions, double[] nodeThicknesses, int[] elements, int[] numNodesPerElement,
        int numMaterials, double[] materials,
        int numSections, double[] sections,
        int[] elementProps,
        int numSupports, int[] supports,
        int numLoads, double[] loadForces, double[] loadMoments, int[] loadNodes,
        bool includeSelfWeight, double[] gravity)
    {
        var handle = ShellInterop.theseus_shell_create(
            (nuint)numNodes, (nuint)numElements,
            nodePositions, nodeThicknesses, elements, numNodesPerElement,
            (nuint)numMaterials, materials,
            (nuint)numSections, sections,
            elementProps,
            (nuint)numSupports, supports,
            (nuint)numLoads, loadForces, loadMoments, loadNodes,
            (byte)(includeSelfWeight ? 1 : 0), gravity);

        if (handle == IntPtr.Zero)
            throw new Exception($"Failed to create shell solver handle. {GetLastError()}");

        return new ShellSolver(handle, numNodes, numElements);
    }

    public ShellSolverResult SolveForward()
    {
        ThrowIfDisposed();
        int nn = _numNodes;
        int ne = _numElements;
        var displacements = new double[nn * 6];
        var reactions = new double[nn * 6];
        var principalStresses = new double[ne * 2];
        var topStresses = new double[ne * 6];
        var bottomStresses = new double[ne * 6];
        var vonMises = new double[ne];

        Check(ShellInterop.theseus_shell_solve_forward(
            _handle, displacements, reactions, principalStresses,
            topStresses, bottomStresses, vonMises));

        return new ShellSolverResult
        {
            Displacements = displacements,
            Reactions = reactions,
            PrincipalStresses = principalStresses,
            TopStresses = topStresses,
            BottomStresses = bottomStresses,
            VonMises = vonMises,
        };
    }

    public void Dispose()
    {
        if (_disposed) return;
        _disposed = true;
        if (_handle != IntPtr.Zero)
        {
            ShellInterop.theseus_shell_destroy(_handle);
            _handle = IntPtr.Zero;
        }
        GC.SuppressFinalize(this);
    }
}
