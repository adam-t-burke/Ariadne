namespace Ariadne.Solver;

using System;
using System.Threading;

internal sealed class PerSolveCancellation : IDisposable
{
    private readonly object _gate = new();
    private readonly Action _requestCancellation;
    private readonly Action _completeScope;
    private CancellationTokenRegistration _registration;
    private bool _completed;

    private PerSolveCancellation(
        CancellationToken cancellationToken,
        Action requestCancellation,
        Action completeScope)
    {
        _requestCancellation = requestCancellation;
        _completeScope = completeScope;
        _registration = cancellationToken.Register(
            static state => ((PerSolveCancellation)state!).RequestCancellation(),
            this);
    }

    public static PerSolveCancellation Begin(
        CancellationToken cancellationToken,
        Action requestCancellation,
        Action completeScope)
    {
        ArgumentNullException.ThrowIfNull(requestCancellation);
        ArgumentNullException.ThrowIfNull(completeScope);
        return new PerSolveCancellation(
            cancellationToken,
            requestCancellation,
            completeScope);
    }

    public static CancellationTokenRegistration Register(
        CancellationToken cancellationToken,
        Action requestCancellation)
    {
        ArgumentNullException.ThrowIfNull(requestCancellation);
        return cancellationToken.Register(
            static state => ((Action)state!).Invoke(),
            requestCancellation);
    }

    private void RequestCancellation()
    {
        lock (_gate)
        {
            if (!_completed)
                _requestCancellation();
        }
    }

    public void Complete()
    {
        lock (_gate)
        {
            if (_completed)
                return;
            _completed = true;
            _completeScope();
        }
    }

    public void Dispose()
    {
        Complete();
        _registration.Dispose();
    }
}
