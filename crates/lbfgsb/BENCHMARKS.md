# Performance baselines

The lean `benches/performance.rs` harness has no benchmark-framework
dependency. Its default sweep covers every objective family once; `--full`
runs every combination of:

- objectives: diagonal and neighbor-rotated quadratics, extended Rosenbrock,
  active corner, mixed bounds, and tied breakpoints;
- dimensions: 25, 180, 512, 1000, 1984, and 4096;
- history sizes: 5, 10, and 20.

Run a representative deterministic release sweep with:

```powershell
$env:LBFGSB_BENCH_MILLIS = '1000'
cargo bench -p ariadne-lbfgsb --bench performance -- --backend=deterministic
```

Use `--full` for the 108-case matrix. Normal `cargo test` does not compile or
run this harness. `LBFGSB_BENCH_MILLIS` is the minimum sample time per case.
Use `--objective=<name>`, `--n=<dimension>`, and `--m=<history>` to benchmark
an exact case.
The harness prints the resolved backend plus termination, final objective,
projected-gradient norm, and weighted `x` checksum. Before timing, it compares
those outputs against the deterministic backend.

## Baseline recorded 2026-08-28 (rerun)

Environment:

- Windows 10.0.26200, Intel Core Ultra 9 285K, 63.5 GiB RAM
- rustc 1.89.0 (LLVM 20.1.7), Cargo 1.89.0
- release profile: `opt-level=3`, LTO, one codegen unit
- one-second minimum per case, warm reusable `Solver`; values are arithmetic
  mean wall-clock microseconds per complete solve

```text
engine objective          n    m     us/solve  iterations evaluations
solver diagonal          25    5       36.376          28          32
solver rotated          180   10        6.309           2           3
solver rosenbrock       512   20      903.289          22          27
solver active-corner   1000    5       24.145           1           2
solver mixed-bounds    1984   10     2026.600          20          22
solver tied-breakpoints 4096  20       90.152           1           2
```

## Official Fortran comparison

`tests/reference/benchmark.ps1` downloads the official L-BFGS-B 3.0 archive,
verifies the pinned SHA-256 documented in `tests/reference/PROVENANCE.md`,
and compiles the unmodified upstream routines with the local
`benchmark_driver.f90`. Sources and executables remain temporary.

```powershell
.\crates\lbfgsb\tests\reference\benchmark.ps1 -Repeats 1000
```

With GNU Fortran 15.2.0, `-O3 -fno-fast-math`, the same machine measured:

```text
engine             n    m     us/solve  iterations evaluations
official Fortran  25    5       51.957          23          28
official Fortran 1000  10     1149.447          24          29
Rust scalar       25    5       29.608          23          28
Rust scalar      1000  10     1266.347          24          29
Rust faer        1000  10     1154.362          24          29
```

The official comparison is the extended Rosenbrock driver with `factr=1e7`
and `pgtol=1e-5`. The Rust and Fortran executables use different compilers and
timers, and runs were not pinned to a CPU core, so small differences should not
be treated as statistically significant. At `n=25,m=5`, scalar Rust was 43.0%
faster than Fortran in this run. At `n=1000,m=10`, scalar Rust was 10.2% slower
and faer Rust was 0.4% slower. Faer improved the Rust result by 8.8% for that
case. The large-n scalar dense-history path remains the clearest solver
bottleneck.

## Dense-history backend recorded 2026-08-28

Profiling showed that `matupd` recomputed every retained `S^T Y` and `S^T S`
entry after each accepted update. The shared kernel path now shifts retained
logical entries when the physical ring wraps and computes only the new row and
column. The scalar kernel preserves the original accumulation order. With the
optional `faer-backend` feature, faer 0.17 performs the two batched
matrix-vector products (`Y^T s_new` and `S^T s_new`) sequentially. A wrapped
ring is passed as at most two native column-major views; no history is repacked
and no temporary allocation is made.

The following extended-Rosenbrock results used explicit deterministic and faer
options, a 100 ms minimum per cell, and otherwise the environment and release
settings above. Times are microseconds per solve.

```text
n     m   scalar      faer    faer improvement
180   5    157.768   148.890        5.6%
180  10    265.833   258.924        2.6%
180  20    436.783   438.633       -0.4%
512   5    442.912   417.936        5.6%
512  10    675.442   624.235        7.6%
512  20    910.973   849.202        6.8%
1000  5    900.686   843.339        6.4%
1000 10   1258.693  1175.790        6.6%
1000 20   1697.336  1633.326        3.8%
1984  5   1957.113  1774.849        9.3%
1984 10   3416.273  3062.012       10.4%
1984 20   4142.388  3905.500        5.7%
4096  5   4073.516  3794.711        6.8%
4096 10   5452.363  5349.500        1.9%
4096 20   8344.338  8111.354        2.8%
```

Faer improved 14 of the 15 extended-Rosenbrock cells in this rerun. The gains
ranged from 1.9% to 10.4%; `n=180,m=20` was 0.4% slower. At the directly
comparable `n=1000,m=10` case, faer reduced Rust solve time from 1258.693 us to
1175.790 us, close to the separately measured 1149.447 us Fortran result.

`Auto` intentionally resolves to the deterministic scalar backend. This
Rosenbrock-only sweep is insufficient evidence for a general size crossover;
callers may explicitly select faer when representative measurements for their
objective justify it. `crates/theseus` does not enable the feature because its
recorded end-to-end runs are objective/gradient dominated and no material
application-level benefit has been measured.

Remaining solver costs are gather-heavy free-variable compact updates,
generalized-Cauchy breakpoint traversal, tiny compact factorizations, and
line-search/objective work. They remain scalar because they are non-contiguous
or too small for this backend.

## Large-n cache-layout improvements recorded 2026-09-02

The benchmark range now extends to 8,192 variables. Two arithmetic-preserving
loop transformations stream complete history columns in `cauchy` and `subsm`
instead of repeatedly loading them with an `n`-element stride. Internal
projected-gradient and active-variable bookkeeping now share one validated
pass. Deterministic objective values, checksums, termination, iteration counts,
and evaluation counts remained unchanged; all reference and allocation tests
passed.

Each value below is the arithmetic mean of three independent one-second
process runs with `m=10`. Baselines were recorded from the parent commit before
the transformations.

```text
backend        objective       n    baseline us   final us   improvement
deterministic  rosenbrock   1000       1232.556     992.171       19.5%
deterministic  rosenbrock   4096       5502.381    4385.175       20.3%
deterministic  rosenbrock   8192      12618.035    9854.659       21.9%
deterministic  mixed-bounds 4096       4382.217    3338.271       23.8%
faer           rosenbrock   1000       1131.143     888.934       21.4%
faer           rosenbrock   4096       5290.724    3943.184       25.5%
faer           rosenbrock   8192      11787.944    8937.428       24.2%
faer           mixed-bounds 4096       4108.152    3008.746       26.8%
```

An attempted extension of faer into subspace products was rejected because it
changed the accepted trajectory from 24 to 25 iterations in the backend
equivalence test. The retained transformations preserve scalar accumulation
order and improve both backends.

## Theseus end-to-end context

The existing ignored release benchmark was run with:

```powershell
cargo test --release -p theseus --test bench_release `
  bench_full_optimize_scaling -- --ignored --nocapture
```

The current rerun measured 95.43 ms/run for the 10x10 grid (180 edges, 20
runs) and 5.54 s/run for the 32x32 grid (1,984 edges, 5 runs).

This existing benchmark uses `DirectSoftBounds`, so it measures Theseus's
argmin path and is not a valid scalar-versus-faer comparison for
`ariadne-lbfgsb`. The previous backend comparison has therefore been removed.
The standalone matrix above remains useful solver-level evidence, but backend
selection in Theseus is based on a dedicated DirectBoxBounds measurement.

A dedicated 48x48 DirectBoxBounds benchmark (4,512 edges) now provides that
measurement. Three scalar process runs averaged 176.690 ms/run. Five confirming
faer runs averaged 163.657 ms/run, a 7.4% end-to-end improvement. Theseus
therefore compiles the existing `faer-backend` feature and explicitly selects
faer for DirectBoxBounds; `Backend::Auto` remains deterministic.

For coarse phase attribution, enable `benchmark-instrumentation` and inspect
`Solver::benchmark_timings()`. It reports total solve time, measured callback
and kernel phases, and a residual containing validation, line-search
bookkeeping, vector operations, report construction, and timer overhead. Timer
reads are compiled out when the feature is disabled.
