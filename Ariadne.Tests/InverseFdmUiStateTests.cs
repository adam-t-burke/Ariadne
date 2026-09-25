using System.Linq;
using Ariadne.Solver.Components.Experimental;
using Xunit;

namespace Ariadne.Tests;

public sealed class InverseFdmUiStateTests
{
    [Fact]
    public void DefaultsMatchTheBenchmarkPipeline()
    {
        // Unboxed: one sparse LDL of the augmented saddle. A box uses the projected quadratic.
        Assert.Equal(ParticularMode.Tikhonov, InverseFdmUiState.DefaultParticular);
        Assert.Equal(MetricMode.Geometric, InverseFdmUiState.DefaultMetric);
        Assert.Equal(1, InverseFdmUiState.DefaultFrozenIterations);
        Assert.Equal(2, InverseFdmUiState.DefaultGnIterations);
        Assert.False(InverseFdmUiState.DefaultSolveForQ);
        Assert.Equal(Stage2Mode.ActiveSet, InverseFdmUiState.DefaultStage2);
        Assert.True(InverseFdmUiState.DefaultNondimensionalize);
        Assert.Equal(3.0, InverseFdmUiState.DefaultSeedGuardMargin);
        Assert.Equal(0.0, InverseFdmUiState.DefaultLmDamping);
        Assert.Equal(1.0, InverseFdmUiState.DefaultReactionWeight);
        Assert.Equal(
            ActiveInverseEngine.Tikhonov,
            InverseFdmUiState.ResolveEngine(
                LinearAlgebraMode.Direct,
                InverseFdmUiState.DefaultParticular,
                hasEffectiveBounds: false));
        Assert.Equal(
            ActiveInverseEngine.Projected,
            InverseFdmUiState.ResolveEngine(
                LinearAlgebraMode.Direct,
                InverseFdmUiState.DefaultParticular,
                hasEffectiveBounds: true));
    }

    [Theory]
    [InlineData(3, 0)] // Gram (sparse)
    [InlineData(0, 1)] // Moore–Penrose (augmented)
    [InlineData(1, 1)] // Tikhonov (augmented)
    [InlineData(2, 2)] // QR
    [InlineData(4, 3)] // Clarabel
    [InlineData(5, 4)] // Gram (dense)
    public void NativeParticularMethodMatchesTheRustAbi(int particularValue, int expected)
    {
        Assert.Equal(expected, InverseFdmUiState.NativeParticularMethod((ParticularMode)particularValue));
    }

    [Fact]
    public void Stage2MappingAndLabelFollowTheMenu()
    {
        Assert.Equal(
            Theseus.Interop.InverseStage2Method.ActiveSet,
            InverseFdmUiState.NativeStage2Method(Stage2Mode.ActiveSet));
        Assert.Equal(
            Theseus.Interop.InverseStage2Method.Clarabel,
            InverseFdmUiState.NativeStage2Method(Stage2Mode.Clarabel));

        Assert.Equal("QP", InverseFdmUiState.Stage2Label(Stage2Mode.ActiveSet, 3.0, true));
        Assert.Equal("IP · guard off · dim", InverseFdmUiState.Stage2Label(Stage2Mode.Clarabel, 0.0, false));
    }

    [Fact]
    public void DiagnosticWarningsFireOnlyOnRecordedEvents()
    {
        var quiet = new Theseus.Interop.InverseFdmDiagnostics
        {
            Stage1Error = 1e-3,
            UniformSeedError = 2e-3,
        };
        Assert.Empty(InverseFdmUiState.DiagnosticWarnings(quiet, loadNorm: 10.0));

        var noisy = new Theseus.Interop.InverseFdmDiagnostics
        {
            Stage1Error = 5.0,
            UniformSeedError = 0.5,
            UsedUniformSeed = true,
            ActiveSetCapped = 1,
            DegenerateLinearizations = 2,
            ClarabelFallbacks = 1,
            ReactionResidual = 4.0,
        };
        var warnings = InverseFdmUiState.DiagnosticWarnings(noisy, loadNorm: 10.0).ToList();
        Assert.Equal(5, warnings.Count);
        Assert.Contains(warnings, w => w.StartsWith("Seed guard"));
        Assert.Contains(warnings, w => w.Contains("pass limit"));
        Assert.Contains(warnings, w => w.Contains("collapsed edge"));
        Assert.Contains(warnings, w => w.Contains("fell back"));
        Assert.Contains(warnings, w => w.Contains("Realised reaction"));

        // A small realised reaction relative to the load is not a warning.
        var smallReaction = new Theseus.Interop.InverseFdmDiagnostics { ReactionResidual = 1e-6 };
        Assert.Empty(InverseFdmUiState.DiagnosticWarnings(smallReaction, loadNorm: 10.0));
    }

    [Fact]
    public void DiagnosticLinesCoverEveryField()
    {
        var lines = InverseFdmUiState.DiagnosticLines(
            new Theseus.Interop.InverseFdmDiagnostics { FrozenSteps = 1, NewtonSteps = 2 },
            iterations: 3, converged: true, geometricError: 0.5);
        Assert.Contains("frozen_steps = 1", lines);
        Assert.Contains("newton_steps = 2", lines);
        Assert.Contains("iterations = 3", lines);
        Assert.Contains("converged = True", lines);
        Assert.Equal(13, lines.Count);
    }

    [Fact]
    public void EffectiveBoundsMatchNativeFiniteBoxRules()
    {
        Assert.False(InverseFdmUiState.HasEffectiveBounds(
            [0, 0],
            [double.NegativeInfinity, double.NegativeInfinity],
            [double.PositiveInfinity, double.PositiveInfinity]));

        Assert.True(InverseFdmUiState.HasEffectiveBounds([0, 1], [], []));
        Assert.True(InverseFdmUiState.HasEffectiveBounds([], [-2.0], []));
        Assert.True(InverseFdmUiState.HasEffectiveBounds([], [], [3.0]));
    }

    [Fact]
    public void StrictSignDefiniteBoundsSuppressOnlyTheSingularityWarning()
    {
        Assert.True(InverseFdmUiState.HasStrictSignDefiniteBounds([0.1, 2.0], []));
        Assert.True(InverseFdmUiState.HasStrictSignDefiniteBounds([], [-0.1, -2.0]));
        Assert.False(InverseFdmUiState.HasStrictSignDefiniteBounds([0.0], []));
        Assert.False(InverseFdmUiState.HasStrictSignDefiniteBounds([], [0.0]));
        Assert.False(InverseFdmUiState.HasStrictSignDefiniteBounds([-2.0], [3.0]));
    }

    [Theory]
    [InlineData(4, 0)]
    [InlineData(0, 1)]
    [InlineData(1, 2)]
    [InlineData(2, 3)]
    [InlineData(3, 4)]
    [InlineData(5, 5)]
    public void UnconstrainedDirectUsesSelectedEngine(int particularValue, int expectedValue)
    {
        var particular = (ParticularMode)particularValue;
        var expected = (ActiveInverseEngine)expectedValue;
        Assert.Equal(
            expected,
            InverseFdmUiState.ResolveEngine(
                LinearAlgebraMode.Direct,
                particular,
                hasEffectiveBounds: false));
    }

    [Fact]
    public void DirectBoxUsesTheProjectedQuadraticAndKeepsTheParticular()
    {
        ParticularMode selected = InverseFdmUiState.UpdateParticular(
            LinearAlgebraMode.Direct,
            ParticularMode.QrLeastSquares,
            hasEffectiveBounds: true);

        Assert.Equal(ParticularMode.QrLeastSquares, selected);
        Assert.Equal(
            ActiveInverseEngine.Projected,
            InverseFdmUiState.ResolveEngine(
                LinearAlgebraMode.Direct,
                selected,
                hasEffectiveBounds: true));
        Assert.Equal(
            ActiveInverseEngine.Clarabel,
            InverseFdmUiState.ResolveEngine(
                LinearAlgebraMode.Direct,
                selected,
                hasEffectiveBounds: true,
                Stage2Mode.Clarabel));

        Assert.Equal(
            ActiveInverseEngine.QrLeastSquares,
            InverseFdmUiState.ResolveEngine(
                LinearAlgebraMode.Direct,
                selected,
                hasEffectiveBounds: false));
    }

    [Fact]
    public void IterativeEngineIgnoresDirectSelection()
    {
        Assert.Equal(
            ActiveInverseEngine.Lsqr,
            InverseFdmUiState.ResolveEngine(
                LinearAlgebraMode.Iterative,
                ParticularMode.MoorePenrose,
                hasEffectiveBounds: false));
        Assert.Equal(
            ActiveInverseEngine.Spg,
            InverseFdmUiState.ResolveEngine(
                LinearAlgebraMode.Iterative,
                ParticularMode.MoorePenrose,
                hasEffectiveBounds: true));
        Assert.Equal(
            ParticularMode.MoorePenrose,
            InverseFdmUiState.UpdateParticular(
                LinearAlgebraMode.Iterative,
                ParticularMode.MoorePenrose,
                hasEffectiveBounds: true));
    }

    [Theory]
    [InlineData(3)] // Gram
    [InlineData(2)] // QrLeastSquares
    public void GramAndQrCanInitializeIndependentGeometricStage(int particularValue)
    {
        Assert.True(InverseFdmUiState.SupportsGeometricMetric(
            LinearAlgebraMode.Direct,
            (ParticularMode)particularValue,
            hasEffectiveBounds: false));
    }

    [Theory]
    [InlineData(4)] // Clarabel
    [InlineData(0)] // MoorePenrose
    [InlineData(1)] // Tikhonov
    public void GeometricMetricAcceptsLeftWeightableDirectSolvers(int particularValue)
    {
        Assert.True(InverseFdmUiState.SupportsGeometricMetric(
            LinearAlgebraMode.Direct,
            (ParticularMode)particularValue,
            hasEffectiveBounds: false));
    }

    [Theory]
    [InlineData(3)] // Gram
    [InlineData(2)] // QrLeastSquares
    public void BoundsAndIterativeModesRouteAwayFromDensifyingSolvers(int particularValue)
    {
        // A finite box routes Direct to Clarabel and Iterative to SPG, so the
        // geometric metric is available even when the menu still shows Gram/QR.
        var particular = (ParticularMode)particularValue;
        Assert.True(InverseFdmUiState.SupportsGeometricMetric(
            LinearAlgebraMode.Direct,
            particular,
            hasEffectiveBounds: true));
        Assert.True(InverseFdmUiState.SupportsGeometricMetric(
            LinearAlgebraMode.Iterative,
            particular,
            hasEffectiveBounds: false));
    }

    [Fact]
    public void GeometricPhaseBudgetsAreIndependentAndNonnegative()
    {
        Assert.Equal(0, InverseFdmUiState.NativeMetric(MetricMode.Force));
        Assert.Equal(2, InverseFdmUiState.NativeMetric(MetricMode.Geometric));

        Assert.Equal(0, InverseFdmUiState.FrozenIterationBudget(MetricMode.Force, 12));
        Assert.Equal(0, InverseFdmUiState.GaussNewtonIterationBudget(MetricMode.Force, 12));
        Assert.Equal(0, InverseFdmUiState.FrozenIterationBudget(MetricMode.Geometric, -1));
        Assert.Equal(0, InverseFdmUiState.GaussNewtonIterationBudget(MetricMode.Geometric, -1));
        Assert.Equal(7, InverseFdmUiState.FrozenIterationBudget(MetricMode.Geometric, 7));
        Assert.Equal(500, InverseFdmUiState.GaussNewtonIterationBudget(MetricMode.Geometric, 500));
    }

    [Theory]
    [InlineData(0, 0, "Stage 1 only")]
    [InlineData(3, 0, "Frozen×3")]
    [InlineData(0, 3, "GN×3")]
    [InlineData(3, 3, "Frozen×3 → GN×3")]
    public void PhaseLabelReportsActualPipeline(int frozen, int gn, string expected)
    {
        Assert.Equal(expected, InverseFdmUiState.PhaseLabel(frozen, gn));
    }
}
