# Future Issues

## Warm-starting optimization

The solver accepts initial force densities (q) as input. Currently, users can
manually wire the Q output back to the q input for warm-starting, but this
creates a feedback loop in the Grasshopper graph and can cause confusing
behavior (stale values, infinite re-triggering).

**Design options to explore:**
- A dedicated "Q Cache" component that stores the final q values from a
  completed optimization run and outputs them on demand (e.g. via a button).
  This avoids feedback loops and gives explicit user control.
- An auto-warm-start toggle on OptConfig that internally caches the last
  optimized q and uses it as the initial guess for the next run.
- A "snapshot" component that captures the full solver state (q + geometry)
  for reproducibility and warm-starting.

**Key consideration:** warm-starting only helps when the problem is similar
between runs (e.g. small parameter changes). For large topology changes, a
fresh start from uniform q is often better.

## Rust-side cancellation without callback

Currently, cancellation only works at callback invocation points. If
StreamPreview is on and ReportFrequency is high, there can be a delay of up
to ReportFrequency evaluations before cancellation takes effect. For
immediate cancellation regardless of callback, consider:
- A shared atomic flag checked by the Rust solver after each evaluation
- A separate `theseus_cancel(handle)` FFI function

## CancellationToken for forward solves

Forward solves are typically fast enough that cancellation is unnecessary.
If very large networks make forward solves slow, consider adding cancellation
support via a diagnostic perturbation callback or a pre-solve size check.
