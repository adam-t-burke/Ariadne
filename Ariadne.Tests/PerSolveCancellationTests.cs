namespace Ariadne.Tests;

using System;
using System.Threading;
using System.Threading.Tasks;
using Ariadne.Solver;
using Theseus.Interop;
using Xunit;

public sealed class PerSolveCancellationTests
{
    [Fact]
    public void CancellingOneConcurrentSolveDoesNotCancelTheOther()
    {
        using var firstSource = new CancellationTokenSource();
        using var secondSource = new CancellationTokenSource();
        var first = new CancellationTarget();
        var second = new CancellationTarget();

        using var firstRegistration =
            PerSolveCancellation.Register(firstSource.Token, first.RequestCancel);
        using var secondRegistration =
            PerSolveCancellation.Register(secondSource.Token, second.RequestCancel);

        firstSource.Cancel();

        Assert.Equal(1, first.RequestCount);
        Assert.Equal(0, second.RequestCount);
    }

    [Fact]
    public void CancellationAfterSolveRegistrationIsDisposedIsHarmless()
    {
        using var source = new CancellationTokenSource();
        var target = new CancellationTarget();
        var registration = PerSolveCancellation.Register(source.Token, target.RequestCancel);

        registration.Dispose();
        source.Cancel();

        Assert.Equal(0, target.RequestCount);
    }

    [Fact]
    public void AlreadyCancelledSolveIsCancelledWhenRegistered()
    {
        using var source = new CancellationTokenSource();
        source.Cancel();
        var target = new CancellationTarget();

        using var registration =
            PerSolveCancellation.Register(source.Token, target.RequestCancel);

        Assert.Equal(1, target.RequestCount);
    }

    [Fact]
    public void AlreadyCancelledForwardSolveDoesNotEnterNativeWork()
    {
        using var solver = CreateForwardSolver();
        using var source = new CancellationTokenSource();
        source.Cancel();

        Assert.Throws<OperationCanceledException>(
            () => solver.SolveForward(source.Token));
    }

    [Fact]
    public void RunningForwardSolveObservesItsCancellationToken()
    {
        using var solver = CreateForwardSolver();
        solver.SetSelfWeightPrescribed(
            new double[8],
            [0.0, 0.0, 0.0],
            maxIters: 100_000,
            tolerance: 0.0,
            relaxation: 0.5);
        using var source = new CancellationTokenSource();
        source.CancelAfter(TimeSpan.FromMilliseconds(10));

        Assert.Throws<OperationCanceledException>(
            () => solver.SolveForward(source.Token));
    }

    [Fact]
    public void ForwardCancellationAfterManagedCheckIsLatchedAndNextRunIsClean()
    {
        using var solver = CreateForwardSolver();
        using var source = new CancellationTokenSource();

        Assert.Throws<OperationCanceledException>(
            () => solver.SolveForwardWithScopeReadyHook(source.Token, source.Cancel));

        SolverResult cleanResult = solver.SolveForward();
        Assert.True(cleanResult.Converged);
    }

    [Fact]
    public void OptimizeCancellationAfterManagedCheckIsLatchedAndNextRunIsClean()
    {
        using var solver = CreateForwardSolver();
        solver.AddTargetXyz(1.0, [3], [3.0, 0.0, -1.0]);
        using var source = new CancellationTokenSource();

        Assert.Throws<OperationCanceledException>(
            () => solver.OptimizeWithScopeReadyHook(source.Token, source.Cancel));

        SolverResult cleanResult = solver.SolveForward();
        Assert.True(cleanResult.Converged);
    }

    [Fact]
    public void CancellationBeforePublicationIsLatchedForThatScope()
    {
        using var source = new CancellationTokenSource();
        bool pending = false;
        bool activeCancelled = false;
        using var scope = PerSolveCancellation.Begin(
            source.Token,
            () => pending = true,
            () => pending = false);

        source.Cancel();
        Assert.True(pending);

        activeCancelled = pending;
        pending = false;
        Assert.True(activeCancelled);
        scope.Complete();
        Assert.False(pending);
    }

    [Fact]
    public void CancellationDuringActiveScopeReachesActiveRun()
    {
        using var source = new CancellationTokenSource();
        bool activeCancelled = false;
        using var scope = PerSolveCancellation.Begin(
            source.Token,
            () => activeCancelled = true,
            () => { });

        source.Cancel();

        Assert.True(activeCancelled);
    }

    [Fact]
    public async Task CompletionWaitsForInFlightCancellationCallback()
    {
        using var source = new CancellationTokenSource();
        using var callbackEntered = new ManualResetEventSlim();
        using var releaseCallback = new ManualResetEventSlim();
        using var completionStarted = new ManualResetEventSlim();
        int completed = 0;
        using var scope = PerSolveCancellation.Begin(
            source.Token,
            () =>
            {
                callbackEntered.Set();
                releaseCallback.Wait();
            },
            () => Interlocked.Increment(ref completed));

        Task cancellation = Task.Run(source.Cancel);
        Assert.True(callbackEntered.Wait(TimeSpan.FromSeconds(5)));
        Task completion = Task.Run(() =>
        {
            completionStarted.Set();
            scope.Complete();
        });
        Assert.True(completionStarted.Wait(TimeSpan.FromSeconds(5)));
        Assert.False(completion.IsCompleted);

        releaseCallback.Set();
        await Task.WhenAll(cancellation, completion);
        Assert.Equal(1, completed);
    }

    [Fact]
    public void CancellationAfterCompletionCannotPoisonSubsequentScope()
    {
        using var firstSource = new CancellationTokenSource();
        int firstRequests = 0;
        int firstCompletions = 0;
        using (var firstScope = PerSolveCancellation.Begin(
            firstSource.Token,
            () => Interlocked.Increment(ref firstRequests),
            () => Interlocked.Increment(ref firstCompletions)))
        {
            firstScope.Complete();
            firstSource.Cancel();
        }

        using var secondSource = new CancellationTokenSource();
        int secondRequests = 0;
        using var secondScope = PerSolveCancellation.Begin(
            secondSource.Token,
            () => Interlocked.Increment(ref secondRequests),
            () => { });

        Assert.Equal(0, firstRequests);
        Assert.Equal(1, firstCompletions);
        Assert.Equal(0, secondRequests);
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

    private sealed class CancellationTarget
    {
        private int _requestCount;

        public int RequestCount => Volatile.Read(ref _requestCount);

        public void RequestCancel()
        {
            Interlocked.Increment(ref _requestCount);
        }
    }
}
