using System;
using System.Runtime.InteropServices;
using System.Text;

namespace Theseus.Interop;

/// <summary>
/// Result of an FEA solve (forward or optimization).
/// All arrays use row-major layout: displacements[node * 3 + dim].
/// </summary>
public sealed class FeaSolverResult
{
    public double[] Displacements { get; init; } = [];
    public double[] Reactions { get; init; } = [];
    public double[] AxialForces { get; init; } = [];
    public double[] Stresses { get; init; } = [];
    public double[] Strains { get; init; } = [];
    public double[] DeformedXyz { get; init; } = [];
    public double[] Utilization { get; init; } = [];
    public double[] NodePositions { get; init; } = [];
    public double[] Areas { get; init; } = [];
    public int Iterations { get; init; }
    public bool Converged { get; init; }
}

/// <summary>
/// Managed wrapper around the native FEA solver in theseus.dll.
/// </summary>
public sealed class FeaSolver : IDisposable
{
    private IntPtr _handle;
    private readonly int _numNodes;
    private readonly int _numElements;
    private readonly int _numSections;
    private bool _disposed;
    private FeaInterop.FeaNativeProgressCallback? _pinnedCallback;

    private FeaSolver(IntPtr handle, int numNodes, int numElements, int numSections)
    {
        _handle = handle;
        _numNodes = numNodes;
        _numElements = numElements;
        _numSections = numSections;
    }

    ~FeaSolver() { Dispose(); }

    private void ThrowIfDisposed()
    {
        if (_disposed) throw new ObjectDisposedException(nameof(FeaSolver));
    }

    private static void Check(int rc)
    {
        if (rc != 0)
            throw new TheseusException(TheseusSolver.GetLastError(), rc);
    }

    private static nuint[] ToNuint(int[] arr) => Array.ConvertAll(arr, i => (nuint)i);

    // ── Construction ─────────────────────────────────────────

    public static FeaSolver Create(
        int numNodes, int numElements,
        double[] nodePositions, int[] edgeNodes,
        int numMaterials, double[] materials,
        int numSections, double[] sections,
        int[] elementProps,
        int numSupports, int[] supports,
        int numLoads, double[] loadForces, int[] loadNodes,
        bool includeSelfWeight, double[] gravity)
    {
        var handle = FeaInterop.theseus_fea_create(
            (nuint)numNodes, (nuint)numElements,
            nodePositions, ToNuint(edgeNodes),
            (nuint)numMaterials, materials,
            (nuint)numSections, sections,
            ToNuint(elementProps),
            (nuint)numSupports, ToNuint(supports),
            (nuint)numLoads, loadForces, ToNuint(loadNodes),
            includeSelfWeight ? 1 : 0, gravity);

        if (handle == IntPtr.Zero)
            throw new TheseusException(TheseusSolver.GetLastError(), -1);

        return new FeaSolver(handle, numNodes, numElements, numSections);
    }

    // ── Objectives ───────────────────────────────────────────

    public void AddCompliance(double weight)
    {
        ThrowIfDisposed();
        Check(FeaInterop.theseus_fea_add_compliance(_handle, weight));
    }

    public void AddMaxDisplacement(double weight, int[] nodeIndices)
    {
        ThrowIfDisposed();
        Check(FeaInterop.theseus_fea_add_max_displacement(
            _handle, weight, ToNuint(nodeIndices), (nuint)nodeIndices.Length));
    }

    public void AddTargetDisplacement(double weight, int[] nodeIndices, double[] targets)
    {
        ThrowIfDisposed();
        Check(FeaInterop.theseus_fea_add_target_displacement(
            _handle, weight, ToNuint(nodeIndices), (nuint)nodeIndices.Length, targets));
    }

    public void AddMinWeight(double weight)
    {
        ThrowIfDisposed();
        Check(FeaInterop.theseus_fea_add_min_weight(_handle, weight));
    }

    public void AddMaxStress(double weight, int[] edgeIndices, double[] thresholds, double sharpness)
    {
        ThrowIfDisposed();
        Check(FeaInterop.theseus_fea_add_max_stress(
            _handle, weight, ToNuint(edgeIndices), (nuint)edgeIndices.Length, thresholds, sharpness));
    }

    public void AddTargetGeometry(double weight, int[] nodeIndices, double[] targets)
    {
        ThrowIfDisposed();
        Check(FeaInterop.theseus_fea_add_target_geometry(
            _handle, weight, ToNuint(nodeIndices), (nuint)nodeIndices.Length, targets));
    }

    // ── Solver options ───────────────────────────────────────

    public void SetSolverOptions(int maxIterations, double absTol, double relTol,
        double barrierWeight, double barrierSharpness)
    {
        ThrowIfDisposed();
        Check(FeaInterop.theseus_fea_set_solver_options(
            _handle, (nuint)maxIterations, absTol, relTol, barrierWeight, barrierSharpness));
    }

    // ── Variable mask ────────────────────────────────────────

    public void SetVariableMask(
        bool optNodePositions, int[] nodeIndices,
        bool optAreas, int[] areaIndices,
        bool optSupportPositions, int[] supportNodeIndices)
    {
        ThrowIfDisposed();
        Check(FeaInterop.theseus_fea_set_variable_mask(
            _handle,
            optNodePositions ? 1 : 0, ToNuint(nodeIndices), (nuint)nodeIndices.Length,
            optAreas ? 1 : 0, ToNuint(areaIndices), (nuint)areaIndices.Length,
            optSupportPositions ? 1 : 0, ToNuint(supportNodeIndices), (nuint)supportNodeIndices.Length));
    }

    // ── Progress callback ────────────────────────────────────

    internal void SetProgressCallback(FeaInterop.FeaNativeProgressCallback? callback, int frequency)
    {
        ThrowIfDisposed();
        _pinnedCallback = callback;
        Check(FeaInterop.theseus_fea_set_progress_callback(
            _handle, callback, (nuint)Math.Max(1, frequency)));
    }

    // ── Forward solve ────────────────────────────────────────

    public FeaSolverResult SolveForward()
    {
        ThrowIfDisposed();
        int nn = _numNodes;
        int ne = _numElements;
        int nTotal = nn * 3;

        var displacements = new double[nTotal];
        var reactions = new double[nTotal];
        var axialForces = new double[ne];
        var stresses = new double[ne];
        var strains = new double[ne];
        var deformedXyz = new double[nTotal];
        var utilization = new double[ne];

        Check(FeaInterop.theseus_fea_solve_forward(
            _handle, displacements, reactions, axialForces,
            stresses, strains, deformedXyz, utilization));

        return new FeaSolverResult
        {
            Displacements = displacements,
            Reactions = reactions,
            AxialForces = axialForces,
            Stresses = stresses,
            Strains = strains,
            DeformedXyz = deformedXyz,
            Utilization = utilization,
        };
    }

    // ── Optimization ─────────────────────────────────────────

    public FeaSolverResult Optimize()
    {
        ThrowIfDisposed();
        int nn = _numNodes;
        int ne = _numElements;
        int ns = _numSections;
        int nTotal = nn * 3;

        var displacements = new double[nTotal];
        var reactions = new double[nTotal];
        var axialForces = new double[ne];
        var stresses = new double[ne];
        var strains = new double[ne];
        var deformedXyz = new double[nTotal];
        var utilization = new double[ne];
        var nodePositions = new double[nTotal];
        var areas = new double[ns];
        nuint iterations = 0;
        byte converged = 0;

        Check(FeaInterop.theseus_fea_optimize(
            _handle, displacements, reactions, axialForces,
            stresses, strains, deformedXyz, utilization,
            nodePositions, areas, ref iterations, ref converged));

        return new FeaSolverResult
        {
            Displacements = displacements,
            Reactions = reactions,
            AxialForces = axialForces,
            Stresses = stresses,
            Strains = strains,
            DeformedXyz = deformedXyz,
            Utilization = utilization,
            NodePositions = nodePositions,
            Areas = areas,
            Iterations = (int)iterations,
            Converged = converged != 0,
        };
    }

    // ── Disposal ─────────────────────────────────────────────

    public void Dispose()
    {
        if (_disposed) return;
        _disposed = true;
        if (_handle != IntPtr.Zero)
        {
            FeaInterop.theseus_fea_free(_handle);
            _handle = IntPtr.Zero;
        }
        _pinnedCallback = null;
        GC.SuppressFinalize(this);
    }
}
