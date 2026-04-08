using System;
using System.Runtime.InteropServices;
using System.Text;

namespace Theseus.Interop;

/// <summary>
/// Result of an optimisation or forward solve.
/// All arrays use row-major layout: xyz[node * 3 + dim].
/// </summary>
public sealed class SolverResult
{
    public double[] Xyz { get; init; } = [];
    public double[] MemberLengths { get; init; } = [];
    public double[] MemberForces { get; init; } = [];
    public double[] ForceDensities { get; init; } = [];
    public double[] Reactions { get; init; } = [];
    public int Iterations { get; init; }
    public bool Converged { get; init; }
}

/// <summary>
/// Managed wrapper around the native Theseus solver (theseus.dll).
///
/// Implements <see cref="IDisposable"/> to ensure the native handle is freed.
/// A destructor is provided as a safety net for cases where Dispose is not called.
/// </summary>
public sealed class TheseusSolver : IDisposable
{
    private IntPtr _handle;
    private readonly int _numNodes;
    private readonly int _numEdges;
    private bool _disposed;
    private TheseusInterop.NativeProgressCallback? _pinnedCallback;

    private TheseusSolver(IntPtr handle, int numNodes, int numEdges)
    {
        _handle = handle;
        _numNodes = numNodes;
        _numEdges = numEdges;
    }

    ~TheseusSolver()
    {
        Dispose();
    }

    private void ThrowIfDisposed()
    {
        if (_disposed)
            throw new ObjectDisposedException(nameof(TheseusSolver));
    }

    public static string GetLastError()
    {
        var buf = new byte[2048];
        int n = TheseusInterop.theseus_last_error(buf, (nuint)buf.Length);
        if (n <= 0) return string.Empty;
        return Encoding.UTF8.GetString(buf, 0, n);
    }

    private static void Check(int rc)
    {
        if (rc != 0)
            throw new TheseusException(GetLastError(), rc);
    }

    // ── Construction ─────────────────────────────────────────

    public static TheseusSolver Create(
        int numEdges, int numNodes, int numFree,
        int[] cooRows, int[] cooCols, double[] cooVals,
        int[] freeNodeIndices, int[] fixedNodeIndices,
        double[] loads, double[] fixedPositions,
        double[] qInit, double[] lowerBounds, double[] upperBounds)
    {
        int numFixed = fixedNodeIndices.Length;

        var handle = TheseusInterop.theseus_create(
            (nuint)numEdges, (nuint)numNodes, (nuint)numFree,
            ToNuint(cooRows), ToNuint(cooCols), cooVals, (nuint)cooRows.Length,
            ToNuint(freeNodeIndices), ToNuint(fixedNodeIndices), (nuint)numFixed,
            loads, fixedPositions,
            qInit, lowerBounds, upperBounds);

        if (handle == IntPtr.Zero)
            throw new TheseusException(GetLastError(), -1);

        return new TheseusSolver(handle, numNodes, numEdges);
    }

    // ── Objectives ───────────────────────────────────────────

    public void AddTargetXyz(double weight, int[] nodeIndices, double[] targetXyz)
    {
        ThrowIfDisposed();
        Check(TheseusInterop.theseus_add_target_xyz(
            _handle, weight, ToNuint(nodeIndices), (nuint)nodeIndices.Length, targetXyz));
    }

    public void AddTargetXy(double weight, int[] nodeIndices, double[] targetXy)
    {
        ThrowIfDisposed();
        Check(TheseusInterop.theseus_add_target_xy(
            _handle, weight, ToNuint(nodeIndices), (nuint)nodeIndices.Length, targetXy));
    }

    public void AddTargetPlane(double weight, int[] nodeIndices, double[] targetXyz, double[] origin, double[] xAxis, double[] yAxis)
    {
        ThrowIfDisposed();
        if (targetXyz.Length != nodeIndices.Length * 3)
            throw new ArgumentException("targetXyz length must be nodeIndices.Length * 3.", nameof(targetXyz));
        if (origin == null || origin.Length != 3)
            throw new ArgumentException("origin must have length 3.", nameof(origin));
        if (xAxis == null || xAxis.Length != 3)
            throw new ArgumentException("xAxis must have length 3.", nameof(xAxis));
        if (yAxis == null || yAxis.Length != 3)
            throw new ArgumentException("yAxis must have length 3.", nameof(yAxis));
        Check(TheseusInterop.theseus_add_target_plane(
            _handle, weight, ToNuint(nodeIndices), (nuint)nodeIndices.Length,
            targetXyz, origin, xAxis, yAxis));
    }

    public void AddPlanarConstraintAlongDirection(double weight, int[] nodeIndices, double[] origin, double[] xAxis, double[] yAxis, double[] direction)
    {
        ThrowIfDisposed();
        if (origin == null || origin.Length != 3)
            throw new ArgumentException("origin must have length 3.", nameof(origin));
        if (xAxis == null || xAxis.Length != 3)
            throw new ArgumentException("xAxis must have length 3.", nameof(xAxis));
        if (yAxis == null || yAxis.Length != 3)
            throw new ArgumentException("yAxis must have length 3.", nameof(yAxis));
        if (direction == null || direction.Length != 3)
            throw new ArgumentException("direction must have length 3.", nameof(direction));
        Check(TheseusInterop.theseus_add_planar_constraint_along_direction(
            _handle, weight, ToNuint(nodeIndices), (nuint)nodeIndices.Length,
            origin, xAxis, yAxis, direction));
    }

    public void AddTargetLength(double weight, int[] edgeIndices, double[] targets)
    {
        ThrowIfDisposed();
        Check(TheseusInterop.theseus_add_target_length(
            _handle, weight, ToNuint(edgeIndices), (nuint)edgeIndices.Length, targets));
    }

    public void AddLengthVariation(double weight, int[] edgeIndices, double sharpness)
    {
        ThrowIfDisposed();
        Check(TheseusInterop.theseus_add_length_variation(
            _handle, weight, ToNuint(edgeIndices), (nuint)edgeIndices.Length, sharpness));
    }

    public void AddForceVariation(double weight, int[] edgeIndices, double sharpness)
    {
        ThrowIfDisposed();
        Check(TheseusInterop.theseus_add_force_variation(
            _handle, weight, ToNuint(edgeIndices), (nuint)edgeIndices.Length, sharpness));
    }

    public void AddSumForceLength(double weight, int[] edgeIndices)
    {
        ThrowIfDisposed();
        Check(TheseusInterop.theseus_add_sum_force_length(
            _handle, weight, ToNuint(edgeIndices), (nuint)edgeIndices.Length));
    }

    public void AddMinLength(double weight, int[] edgeIndices, double[] thresholds, double sharpness)
    {
        ThrowIfDisposed();
        Check(TheseusInterop.theseus_add_min_length(
            _handle, weight, ToNuint(edgeIndices), (nuint)edgeIndices.Length, thresholds, sharpness));
    }

    public void AddMaxLength(double weight, int[] edgeIndices, double[] thresholds, double sharpness)
    {
        ThrowIfDisposed();
        Check(TheseusInterop.theseus_add_max_length(
            _handle, weight, ToNuint(edgeIndices), (nuint)edgeIndices.Length, thresholds, sharpness));
    }

    public void AddMinForce(double weight, int[] edgeIndices, double[] thresholds, double sharpness)
    {
        ThrowIfDisposed();
        Check(TheseusInterop.theseus_add_min_force(
            _handle, weight, ToNuint(edgeIndices), (nuint)edgeIndices.Length, thresholds, sharpness));
    }

    public void AddMaxForce(double weight, int[] edgeIndices, double[] thresholds, double sharpness)
    {
        ThrowIfDisposed();
        Check(TheseusInterop.theseus_add_max_force(
            _handle, weight, ToNuint(edgeIndices), (nuint)edgeIndices.Length, thresholds, sharpness));
    }

    public void AddRigidSetCompare(double weight, int[] nodeIndices, double[] targetXyz)
    {
        ThrowIfDisposed();
        Check(TheseusInterop.theseus_add_rigid_set_compare(
            _handle, weight, ToNuint(nodeIndices), (nuint)nodeIndices.Length, targetXyz));
    }

    public void AddReactionDirection(double weight, int[] anchorIndices, double[] targetDirs)
    {
        ThrowIfDisposed();
        Check(TheseusInterop.theseus_add_reaction_direction(
            _handle, weight, ToNuint(anchorIndices), (nuint)anchorIndices.Length, targetDirs));
    }

    public void AddReactionDirectionMagnitude(double weight, int[] anchorIndices, double[] targetDirs, double[] targetMags)
    {
        ThrowIfDisposed();
        Check(TheseusInterop.theseus_add_reaction_direction_magnitude(
            _handle, weight, ToNuint(anchorIndices), (nuint)anchorIndices.Length, targetDirs, targetMags));
    }

    // ── Self-weight ───────────────────────────────────────────

    public void SetSelfWeightPrescribed(
        double[] linearDensities, double[] gravity,
        int maxIters = 50, double tolerance = 1e-6, double relaxation = 1.0)
    {
        ThrowIfDisposed();
        Check(TheseusInterop.theseus_set_self_weight_prescribed(
            _handle, linearDensities, gravity,
            (nuint)maxIters, tolerance, relaxation));
    }

    public void SetSelfWeightSizing(
        double rho, double sigma, double[] gravity,
        int maxIters = 50, double tolerance = 1e-6, double relaxation = 1.0)
    {
        ThrowIfDisposed();
        Check(TheseusInterop.theseus_set_self_weight_sizing(
            _handle, rho, sigma, gravity,
            (nuint)maxIters, tolerance, relaxation));
    }

    public void ClearSelfWeight()
    {
        ThrowIfDisposed();
        Check(TheseusInterop.theseus_clear_self_weight(_handle));
    }

    // ── Pressure loads ──────────────────────────────────────

    public void SetPressure(
        int[][] faces, double[] pressures,
        int maxIters = 50, double tolerance = 1e-6, double relaxation = 1.0)
    {
        ThrowIfDisposed();
        int numFaces = faces.Length;
        var offsets = new nuint[numFaces + 1];
        int totalVerts = 0;
        for (int f = 0; f < numFaces; f++)
        {
            offsets[f] = (nuint)totalVerts;
            totalVerts += faces[f].Length;
        }
        offsets[numFaces] = (nuint)totalVerts;

        var verts = new nuint[totalVerts];
        int idx = 0;
        for (int f = 0; f < numFaces; f++)
            foreach (int v in faces[f])
                verts[idx++] = (nuint)v;

        Check(TheseusInterop.theseus_set_pressure(
            _handle, (nuint)numFaces, offsets, verts, pressures,
            (nuint)maxIters, tolerance, relaxation));
    }

    public void SetPressureHydrostatic(
        int[][] faces, double rhoFluid, double gMagnitude, double zDatum,
        double[] upDirection,
        int maxIters = 50, double tolerance = 1e-6, double relaxation = 1.0)
    {
        ThrowIfDisposed();
        int numFaces = faces.Length;
        var offsets = new nuint[numFaces + 1];
        int totalVerts = 0;
        for (int f = 0; f < numFaces; f++)
        {
            offsets[f] = (nuint)totalVerts;
            totalVerts += faces[f].Length;
        }
        offsets[numFaces] = (nuint)totalVerts;

        var verts = new nuint[totalVerts];
        int idx = 0;
        for (int f = 0; f < numFaces; f++)
            foreach (int v in faces[f])
                verts[idx++] = (nuint)v;

        Check(TheseusInterop.theseus_set_pressure_hydrostatic(
            _handle, (nuint)numFaces, offsets, verts,
            rhoFluid, gMagnitude, zDatum, upDirection,
            (nuint)maxIters, tolerance, relaxation));
    }

    public void SetPressureDirectional(
        int[][] faces, double[] pressures, double[] direction,
        int maxIters = 50, double tolerance = 1e-6, double relaxation = 1.0)
    {
        ThrowIfDisposed();
        int numFaces = faces.Length;
        var offsets = new nuint[numFaces + 1];
        int totalVerts = 0;
        for (int f = 0; f < numFaces; f++)
        {
            offsets[f] = (nuint)totalVerts;
            totalVerts += faces[f].Length;
        }
        offsets[numFaces] = (nuint)totalVerts;

        var verts = new nuint[totalVerts];
        int idx = 0;
        for (int f = 0; f < numFaces; f++)
            foreach (int v in faces[f])
                verts[idx++] = (nuint)v;

        Check(TheseusInterop.theseus_set_pressure_directional(
            _handle, (nuint)numFaces, offsets, verts,
            pressures, direction,
            (nuint)maxIters, tolerance, relaxation));
    }

    public void ClearPressure()
    {
        ThrowIfDisposed();
        Check(TheseusInterop.theseus_clear_pressure(_handle));
    }

    // ── Solver options ───────────────────────────────────────

    public void SetSolverOptions(
        int maxIterations = 500,
        double absTol = 1e-6,
        double relTol = 1e-6,
        double barrierWeight = 10.0,
        double barrierSharpness = 10.0)
    {
        ThrowIfDisposed();
        Check(TheseusInterop.theseus_set_solver_options(
            _handle, (nuint)maxIterations, absTol, relTol, barrierWeight, barrierSharpness));
    }

    // ── Progress callback ─────────────────────────────────────

    /// <summary>
    /// Register a managed callback invoked every <paramref name="frequency"/>
    /// evaluations with (iteration, loss, xyz[numNodes*3], q[numEdges]).
    /// Return <c>true</c> to continue, <c>false</c> to cancel.
    /// Pass null to clear.  The delegate is pinned for the lifetime of this solver.
    /// </summary>
    public void SetProgressCallback(Func<int, double, double[], double[], bool>? callback, int frequency)
    {
        ThrowIfDisposed();
        if (callback == null)
        {
            _pinnedCallback = null;
            Check(TheseusInterop.theseus_set_progress_callback(_handle, null, (nuint)1));
            return;
        }

        int nn = _numNodes;
        int ne = _numEdges;
        _pinnedCallback = (nuint iteration, double loss, IntPtr xyzPtr, nuint numNodes, IntPtr qPtr, nuint numEdges) =>
        {
            var xyz = new double[nn * 3];
            Marshal.Copy(xyzPtr, xyz, 0, nn * 3);
            var q = new double[ne];
            Marshal.Copy(qPtr, q, 0, ne);
            bool shouldContinue = callback((int)iteration, loss, xyz, q);
            return shouldContinue ? (byte)1 : (byte)0;
        };

        Check(TheseusInterop.theseus_set_progress_callback(
            _handle, _pinnedCallback, (nuint)Math.Max(1, frequency)));
    }

    // ── Solve ────────────────────────────────────────────────

    public SolverResult Optimize()
    {
        ThrowIfDisposed();
        var xyz = new double[_numNodes * 3];
        var lengths = new double[_numEdges];
        var forces = new double[_numEdges];
        var q = new double[_numEdges];
        var reactions = new double[_numNodes * 3];
        nuint iterations = 0;
        byte converged = 0;

        Check(TheseusInterop.theseus_optimize(
            _handle, xyz, lengths, forces, q, reactions,
            ref iterations, ref converged));

        return new SolverResult
        {
            Xyz = xyz,
            MemberLengths = lengths,
            MemberForces = forces,
            ForceDensities = q,
            Reactions = reactions,
            Iterations = (int)iterations,
            Converged = converged != 0,
        };
    }

    public SolverResult SolveForward()
    {
        ThrowIfDisposed();
        var xyz = new double[_numNodes * 3];
        var lengths = new double[_numEdges];
        var forces = new double[_numEdges];
        var q = new double[_numEdges];
        var reactions = new double[_numNodes * 3];

        Check(TheseusInterop.theseus_solve_forward(_handle, xyz, lengths, forces, q, reactions));

        return new SolverResult
        {
            Xyz = xyz,
            MemberLengths = lengths,
            MemberForces = forces,
            ForceDensities = q,
            Reactions = reactions,
            Iterations = 1,
            Converged = true,
        };
    }

    // ── Inverse solvers (experimental) ──────────────────────

    public SolverResult SolvePseudoinverse(
        double[] targetFreeXyz, double regularization,
        bool useL2 = true, int maxL1Iter = 20, bool useAugmented = false,
        bool enforceZeroRx = false, bool enforceZeroRy = false,
        bool enforceZeroRz = false, bool solveForQ = true)
    {
        ThrowIfDisposed();
        var q = new double[_numEdges];
        var xyz = new double[_numNodes * 3];
        var lengths = new double[_numEdges];
        var forces = new double[_numEdges];
        var reactions = new double[_numNodes * 3];

        Check(TheseusInterop.theseus_solve_pseudoinverse(
            _handle, targetFreeXyz, regularization,
            useL2 ? 1 : 0, (nuint)maxL1Iter, useAugmented ? 1 : 0,
            enforceZeroRx ? 1 : 0, enforceZeroRy ? 1 : 0,
            enforceZeroRz ? 1 : 0, solveForQ ? 1 : 0,
            q, xyz, lengths, forces, reactions));

        return new SolverResult
        {
            Xyz = xyz,
            MemberLengths = lengths,
            MemberForces = forces,
            ForceDensities = q,
            Reactions = reactions,
            Iterations = 1,
            Converged = true,
        };
    }

    public SolverResult SolveNnls(double[] targetFreeXyz, int maxIter, double tol)
    {
        ThrowIfDisposed();
        var q = new double[_numEdges];
        var xyz = new double[_numNodes * 3];
        var lengths = new double[_numEdges];
        var forces = new double[_numEdges];
        var reactions = new double[_numNodes * 3];

        Check(TheseusInterop.theseus_solve_nnls(
            _handle, targetFreeXyz, (nuint)maxIter, tol,
            q, xyz, lengths, forces, reactions));

        return new SolverResult
        {
            Xyz = xyz,
            MemberLengths = lengths,
            MemberForces = forces,
            ForceDensities = q,
            Reactions = reactions,
            Iterations = maxIter,
            Converged = true,
        };
    }

    // ── IDisposable ──────────────────────────────────────────

    public void Dispose()
    {
        if (!_disposed && _handle != IntPtr.Zero)
        {
            TheseusInterop.theseus_free(_handle);
            _handle = IntPtr.Zero;
            _disposed = true;
        }
        GC.SuppressFinalize(this);
    }

    // ── Helpers ──────────────────────────────────────────────

    private static nuint[] ToNuint(int[] arr)
    {
        var result = new nuint[arr.Length];
        for (int i = 0; i < arr.Length; i++)
            result[i] = (nuint)arr[i];
        return result;
    }
}

/// <summary>
/// Exception thrown when the native Theseus library returns an error.
/// </summary>
public class TheseusException : Exception
{
    public int NativeCode { get; }

    public TheseusException(string message, int nativeCode)
        : base(message)
    {
        NativeCode = nativeCode;
    }
}
