namespace Ariadne.Tests;

using System;
using Ariadne.Solver;
using Theseus.Interop;
using Xunit;

public sealed class LinearSolverKindTests
{
    // ── Enum values (must match the native i32 encoding) ─────────

    [Fact]
    public void EnumValuesMatchNativeEncoding()
    {
        Assert.Equal(0, (int)LinearSolverKind.Direct);
        Assert.Equal(1, (int)LinearSolverKind.IterativeCpu);
        Assert.Equal(2, (int)LinearSolverKind.IterativeGpu);
        Assert.Equal(3, Enum.GetValues<LinearSolverKind>().Length);
    }

    [Fact]
    public void IterativeOptionEnumsMatchNativeEncoding()
    {
        Assert.Equal(0, (int)IterativeToleranceMode.Fixed);
        Assert.Equal(1, (int)IterativeToleranceMode.Adaptive);

        Assert.Equal(0, (int)MultigridCycle.V);
        Assert.Equal(1, (int)MultigridCycle.K);

        Assert.Equal(0, (int)PreconditionerPrecision.F64);
        Assert.Equal(1, (int)PreconditionerPrecision.F32);

        Assert.Equal(0, (int)GpuOuterLoop.Auto);
        Assert.Equal(1, (int)GpuOuterLoop.Device);
        Assert.Equal(2, (int)GpuOuterLoop.Host);

        Assert.Equal(0, (int)GpuAdapterPreference.Discrete);
        Assert.Equal(1, (int)GpuAdapterPreference.Integrated);
        Assert.Equal(2, (int)GpuAdapterPreference.Any);
    }

    // ── Persistence ──────────────────────────────────────────────

    [Fact]
    public void PersistenceKeyIsStable()
    {
        Assert.Equal("LinearSolverKind", LinearSolverSelection.PersistenceKey);
    }

    [Fact]
    public void MissingKeyLoadsAsDirect()
    {
        Assert.Equal(LinearSolverKind.Direct, LinearSolverSelection.Default);
        Assert.Equal(LinearSolverKind.Direct, LinearSolverSelection.FromPersisted(null));
    }

    [Theory]
    [InlineData(0, LinearSolverKind.Direct)]
    [InlineData(1, LinearSolverKind.IterativeCpu)]
    [InlineData(2, LinearSolverKind.IterativeGpu)]
    public void KnownPersistedValuesRoundTrip(int stored, LinearSolverKind expected)
    {
        Assert.Equal(expected, LinearSolverSelection.FromPersisted(stored));
        Assert.Equal(stored, (int)expected);
    }

    [Theory]
    [InlineData(-1)]
    [InlineData(3)]
    [InlineData(int.MaxValue)]
    public void UnknownPersistedValuesFallBackToDirect(int stored)
    {
        Assert.Equal(LinearSolverKind.Direct, LinearSolverSelection.FromPersisted(stored));
    }

    [Fact]
    public void MenuOrderAndLabels()
    {
        Assert.Equal(
            [LinearSolverKind.Direct, LinearSolverKind.IterativeCpu, LinearSolverKind.IterativeGpu],
            LinearSolverSelection.All);
        Assert.Equal("Direct", LinearSolverSelection.MenuLabel(LinearSolverKind.Direct));
        Assert.Equal("Iterative (CPU)", LinearSolverSelection.MenuLabel(LinearSolverKind.IterativeCpu));
        Assert.Equal("Iterative (GPU)", LinearSolverSelection.MenuLabel(LinearSolverKind.IterativeGpu));
    }

    [Fact]
    public void MessageSuffixes()
    {
        Assert.Equal("lin: Direct", LinearSolverSelection.MessageSuffix(LinearSolverKind.Direct));
        Assert.Equal("lin: Iter-CPU", LinearSolverSelection.MessageSuffix(LinearSolverKind.IterativeCpu));
        Assert.Equal("lin: Iter-GPU", LinearSolverSelection.MessageSuffix(LinearSolverKind.IterativeGpu));
    }

    [Fact]
    public void ConfigDefaultsToDirectWithoutOptions()
    {
        var config = new OptimizationConfig { Objectives = [] };
        Assert.Equal(LinearSolverKind.Direct, config.LinearSolver);
        Assert.Null(config.IterativeOptions);
    }

    // ── IterativeSolverOptions defaults (match IterativeSolverOptions::default()) ──

    [Fact]
    public void IterativeOptionsDefaultsMatchNative()
    {
        var o = new IterativeSolverOptions();

        Assert.Equal(IterativeToleranceMode.Adaptive, o.ToleranceMode);
        Assert.Equal(1e-8, o.Tolerance);
        Assert.Equal(1e-10, o.ToleranceFloor);
        Assert.Equal(1e-6, o.ToleranceCeiling);
        Assert.Equal(1e-2, o.ToleranceFactor);
        Assert.Equal(200u, o.MaxIterations);
        Assert.Equal(MultigridCycle.K, o.Cycle);
        Assert.Equal(2u, o.SmootherDegree);
        Assert.Equal(2u, o.AggregationPasses);
        Assert.Equal(2000u, o.CoarsestSize);
        Assert.Null(o.PreconditionPrecision);
        Assert.Equal(-1, o.NativePreconditionPrecision);
        Assert.Equal(GpuOuterLoop.Auto, o.GpuOuterLoop);
        Assert.Equal(GpuAdapterPreference.Discrete, o.AdapterPreference);
    }

    [Fact]
    public void ExplicitPrecisionMapsToNativeCode()
    {
        Assert.Equal(0, new IterativeSolverOptions { PreconditionPrecision = PreconditionerPrecision.F64 }.NativePreconditionPrecision);
        Assert.Equal(1, new IterativeSolverOptions { PreconditionPrecision = PreconditionerPrecision.F32 }.NativePreconditionPrecision);
    }

    [Fact]
    public void ContentHashDistinguishesOptions()
    {
        var a = new IterativeSolverOptions();
        var b = new IterativeSolverOptions();
        var c = new IterativeSolverOptions { MaxIterations = 201 };
        var d = new IterativeSolverOptions { PreconditionPrecision = PreconditionerPrecision.F32 };

        Assert.Equal(a.GetContentHashCode(), b.GetContentHashCode());
        Assert.NotEqual(a.GetContentHashCode(), c.GetContentHashCode());
        Assert.NotEqual(a.GetContentHashCode(), d.GetContentHashCode());
    }

    [Fact]
    public void IterativeOptionsInputValidationMirrorsNativeRanges()
    {
        Assert.Null(IterativeSolverOptionsInput.Validate(false, double.NaN, 1e-10, 1e-6, 1e-2, 200, 2, 2, 2000));
        Assert.Null(IterativeSolverOptionsInput.Validate(true, 1e-8, 0, 0, 0, 1, 1, 1, 1));

        Assert.NotNull(IterativeSolverOptionsInput.Validate(true, 0.0, 1e-10, 1e-6, 1e-2, 200, 2, 2, 2000));
        Assert.NotNull(IterativeSolverOptionsInput.Validate(true, double.NaN, 1e-10, 1e-6, 1e-2, 200, 2, 2, 2000));
        Assert.NotNull(IterativeSolverOptionsInput.Validate(false, 0, 1e-6, 1e-10, 1e-2, 200, 2, 2, 2000)); // floor > ceiling
        Assert.NotNull(IterativeSolverOptionsInput.Validate(false, 0, -1e-10, 1e-6, 1e-2, 200, 2, 2, 2000));
        Assert.NotNull(IterativeSolverOptionsInput.Validate(false, 0, 1e-10, 1e-6, 1e-2, 0, 2, 2, 2000));
        Assert.NotNull(IterativeSolverOptionsInput.Validate(false, 0, 1e-10, 1e-6, 1e-2, 200, 0, 2, 2000));
        Assert.NotNull(IterativeSolverOptionsInput.Validate(false, 0, 1e-10, 1e-6, 1e-2, 200, 256, 2, 2000));
        Assert.NotNull(IterativeSolverOptionsInput.Validate(false, 0, 1e-10, 1e-6, 1e-2, 200, 2, 0, 2000));
        Assert.NotNull(IterativeSolverOptionsInput.Validate(false, 0, 1e-10, 1e-6, 1e-2, 200, 2, 2, 0));
    }

    // ── GPU probe JSON ───────────────────────────────────────────

    [Fact]
    public void ProbeParsesUnavailableDocument()
    {
        const string json = "{\"available\":false,\"adapters\":[],\"chosen\":null,\"reason\":\"GPU backend not yet built\"}";

        var probe = GpuProbe.Parse(json);

        Assert.False(probe.Available);
        Assert.Empty(probe.Adapters);
        Assert.Null(probe.ChosenIndex);
        Assert.Null(probe.ChosenAdapter);
        Assert.Equal("GPU backend not yet built", probe.Reason);
        Assert.Equal(json, probe.Json);
        Assert.Contains("GPU backend not yet built", probe.Describe());
    }

    [Fact]
    public void ProbeParsesAvailableDocumentWithAdapters()
    {
        const string json =
            "{\"available\":true,\"adapters\":[" +
            "{\"name\":\"Intel(R) UHD\",\"backend\":\"vulkan\",\"device_type\":\"integrated\",\"shader_f64\":false,\"max_storage_buffer_binding_size\":1073741824}," +
            "{\"name\":\"NVIDIA \\\"RTX\\\"\",\"backend\":\"vulkan\",\"device_type\":\"discrete\",\"shader_f64\":true,\"max_storage_buffer_binding_size\":4294967296}" +
            "],\"chosen\":1,\"reason\":null}";

        var probe = GpuProbe.Parse(json);

        Assert.True(probe.Available);
        Assert.Equal(2, probe.Adapters.Count);
        Assert.Equal(1, probe.ChosenIndex);
        Assert.Null(probe.Reason);

        var chosen = probe.ChosenAdapter;
        Assert.NotNull(chosen);
        Assert.Equal("NVIDIA \"RTX\"", chosen!.Name);
        Assert.Equal("vulkan", chosen.Backend);
        Assert.Equal("discrete", chosen.DeviceType);
        Assert.True(chosen.ShaderF64);
        Assert.Equal(4294967296UL, chosen.MaxStorageBufferBindingSize);

        Assert.False(probe.Adapters[0].ShaderF64);
        Assert.Equal(1073741824UL, probe.Adapters[0].MaxStorageBufferBindingSize);
        Assert.Contains("NVIDIA", probe.Describe());
    }

    [Fact]
    public void ProbeToleratesMissingOptionalFields()
    {
        var probe = GpuProbe.Parse("{\"available\":false}");

        Assert.False(probe.Available);
        Assert.Empty(probe.Adapters);
        Assert.Null(probe.ChosenIndex);
        Assert.Null(probe.Reason);
    }

    [Fact]
    public void ProbeRejectsNonObjectDocuments()
    {
        Assert.Throws<FormatException>(() => GpuProbe.Parse("[]"));
    }

    [Fact]
    public void OutOfRangeChosenIndexYieldsNoAdapter()
    {
        var probe = GpuProbe.Parse("{\"available\":true,\"adapters\":[],\"chosen\":0,\"reason\":null}");
        Assert.Null(probe.ChosenAdapter);
    }

    // ── Error classification ─────────────────────────────────────

    [Theory]
    [InlineData(-3, true)]  // IterativeSolverDidNotConverge
    [InlineData(-4, true)]  // IterativeSolverUnsupported
    [InlineData(-5, true)]  // GpuUnavailable
    [InlineData(-6, true)]  // GpuOutOfMemory
    [InlineData(-1, false)]
    [InlineData(-2, false)]
    [InlineData(-7, false)]
    [InlineData(0, false)]
    public void LinearSolverErrorCodesAreClassified(int code, bool expected)
    {
        Assert.Equal(expected, TheseusSolverService.IsLinearSolverErrorCode(code));
    }

    // ── Native round trips (need the native library, like PerSolveCancellationTests) ──

    [Fact]
    public void NativeProbeReportsUnavailableGpuInThisBuild()
    {
        var probe = GpuProbe.Query();

        Assert.False(probe.Available);
        Assert.Empty(probe.Adapters);
        Assert.Null(probe.ChosenIndex);
        Assert.False(string.IsNullOrEmpty(probe.Reason));
        Assert.StartsWith("{", probe.Json.TrimStart());
    }

    [Fact]
    public void NativeDirectForwardSolveReportsDirectStats()
    {
        using var solver = CreateForwardSolver();

        var result = solver.SolveForward();

        Assert.NotNull(result.LinearSolver);
        Assert.Equal(LinearSolverKind.Direct, result.LinearSolver!.Backend);
        Assert.True(result.LinearSolver.ConvergedAll);
        Assert.Equal(0UL, result.LinearSolver.IterationsTotal);

        var again = solver.GetLinearSolverStats();
        Assert.Equal(result.LinearSolver, again);
    }

    [Theory]
    [InlineData(LinearSolverKind.IterativeCpu)]
    [InlineData(LinearSolverKind.IterativeGpu)]
    public void NativeIterativeKindsFailWithUnsupportedCodeAndNeverFallBack(LinearSolverKind kind)
    {
        using var solver = CreateForwardSolver();
        solver.SetLinearSolver(kind);
        solver.SetIterativeOptions(new IterativeSolverOptions());

        var forward = Assert.Throws<TheseusException>(() => solver.SolveForward());
        Assert.Equal(-4, forward.NativeCode);
        Assert.True(TheseusSolverService.IsLinearSolverErrorCode(forward.NativeCode));

        var optimize = Assert.Throws<TheseusException>(() => solver.Optimize());
        Assert.Equal(-4, optimize.NativeCode);

        Assert.Equal(kind, solver.GetLinearSolverStats().Backend);

        solver.SetLinearSolver(LinearSolverKind.Direct);
        var result = solver.SolveForward();
        Assert.Equal(LinearSolverKind.Direct, result.LinearSolver!.Backend);
    }

    [Fact]
    public void NativeSetIterativeOptionsRejectsBadValues()
    {
        using var solver = CreateForwardSolver();

        var ex = Assert.Throws<TheseusException>(() =>
            solver.SetIterativeOptions(new IterativeSolverOptions { SmootherDegree = 0 }));
        Assert.NotEqual(0, ex.NativeCode);

        Assert.Throws<ArgumentOutOfRangeException>(() => solver.SetLinearSolver((LinearSolverKind)7));
    }

    private static TheseusSolver CreateForwardSolver()
    {
        (int Start, int End)[] edges =
        [
            (0, 1), (1, 2), (2, 3), (3, 4),
            (4, 5), (5, 6), (1, 5), (2, 4),
        ];
        var rows = new int[edges.Length * 2];
        var cols = new int[edges.Length * 2];
        var vals = new double[edges.Length * 2];
        for (int i = 0; i < edges.Length; i++)
        {
            rows[i * 2] = i;
            cols[i * 2] = edges[i].Start;
            vals[i * 2] = -1.0;
            rows[i * 2 + 1] = i;
            cols[i * 2 + 1] = edges[i].End;
            vals[i * 2 + 1] = 1.0;
        }

        return TheseusSolver.Create(
            edges.Length,
            numNodes: 7,
            numFree: 5,
            rows,
            cols,
            vals,
            freeNodeIndices: [1, 2, 3, 4, 5],
            fixedNodeIndices: [0, 6],
            loads:
            [
                0.0, 0.0, -1.0,
                0.0, 0.0, -1.0,
                0.0, 0.0, -2.0,
                0.0, 0.0, -1.0,
                0.0, 0.0, -1.0,
            ],
            fixedPositions: [0.0, 0.0, 0.0, 6.0, 0.0, 0.0],
            qInit: [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
            lowerBounds: [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1],
            upperBounds: [100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0]);
    }

    [Fact]
    public void StatsDescribeMentionsBackendAndIterations()
    {
        var stats = new LinearSolverStats
        {
            Backend = LinearSolverKind.IterativeCpu,
            Solves = 12,
            IterationsTotal = 340,
            IterationsMax = 41,
            SolveMsTotal = 12.5,
            SetupMsTotal = 2.5,
            ConvergedAll = false,
        };

        string text = stats.Describe();
        Assert.Contains("IterativeCpu", text);
        Assert.Contains("12 solve(s)", text);
        Assert.Contains("340 iteration(s)", text);
        Assert.Contains("15.0 ms", text);
        Assert.Contains("NOT all converged", text);
    }
}
