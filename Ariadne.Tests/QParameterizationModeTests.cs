namespace Ariadne.Tests;

using System;
using Ariadne.Solver;
using Xunit;

public sealed class QParameterizationModeTests
{
    [Fact]
    public void LegacyImplicitModeMigratesToBoxBoundsWhenAllBoundsAreFinite()
    {
        var mode = TheseusSolverService.ResolveQParameterizationMode(
            (QParameterizationMode)1,
            [0.1, 0.2],
            [10.0, 20.0]);

        Assert.Equal(QParameterizationMode.DirectBoxBounds, mode);
    }

    [Fact]
    public void LegacyImplicitModeMigratesToSoftBoundsWhenAnyBoundIsInfinite()
    {
        var mode = TheseusSolverService.ResolveQParameterizationMode(
            (QParameterizationMode)1,
            [0.1, 0.2],
            [10.0, double.PositiveInfinity]);

        Assert.Equal(QParameterizationMode.DirectSoftBounds, mode);
    }

    [Fact]
    public void CurrentModesPassThroughUnchanged()
    {
        var lower = new[] { 0.1 };
        var upper = new[] { 10.0 };

        Assert.Equal(
            QParameterizationMode.DirectSoftBounds,
            TheseusSolverService.ResolveQParameterizationMode(
                QParameterizationMode.DirectSoftBounds,
                lower,
                upper));
        Assert.Equal(
            QParameterizationMode.DirectBoxBounds,
            TheseusSolverService.ResolveQParameterizationMode(
                QParameterizationMode.DirectBoxBounds,
                lower,
                upper));
    }
}
