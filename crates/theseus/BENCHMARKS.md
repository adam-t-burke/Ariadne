# Theseus benchmarks

1. [Basin optimizer comparison](#basin-optimizer-comparison)
2. [Null-space paper benchmarks](#null-space-paper-benchmarks)

# Basin optimizer comparison

This compares the integration for [issue #12](https://github.com/adam-t-burke/Ariadne/issues/12)
with commit `c8e957cfe58648f07ce4222e4f6e925563d834b2`, which uses argmin
0.10 for soft bounds and ariadne-lbfgsb 0.1.0 with its faer backend for box
bounds. The replacement uses Basin 1.13.0 with `Vec<f64>` parameters and no
optional features. The standalone ariadne-lbfgsb crate is preserved.

## Reproduce

From the repository root, with Python 3, Git, and Rust available:

```sh
python scripts/compare-optimizers.py
```

The script creates a temporary worktree at the baseline commit, copies the
same benchmark cases into it, and builds both revisions before timing either
one. It runs the bounded reference checks on both revisions and removes its
temporary worktree afterward. `--baseline <revision>` selects another
compatible pre-Basin revision. The current checkout may contain uncommitted
implementation changes.

Results below were measured on September 18, 2026, on NixOS with an AMD
Ryzen 9 7900 (12 cores, 24 threads), Rust 1.98.1, and Cargo 1.98.1. Both
revisions use the workspace release profile (`opt-level = 3`, LTO, one
codegen unit) and `RAYON_NUM_THREADS=1`. Each case has one warmup and ten
timed samples; the reported time is their median.
Reference samples batch repeated solves for at least 50 ms and report the
time per solve. Grid samples time a single solve, including optimization
state/cache initialization and result construction, but exclude network
construction and the separate final objective/gradient check. Baseline cases
run before Basin cases. An initial run overlapped unrelated Nix builds and
was discarded. The reported run started after those builds stopped, with
total CPU utilization below 1% before compilation. These are local
measurements without CPU pinning; small timing differences should be
treated cautiously.

## Full Theseus solves

The grid cases use the existing `bench_release.rs` topology: four fixed
corners, uniform downward loads, and a target grid at `z = -0.2`. Both
revisions start every force density at 1, use history length 10, allow 200
accepted iterations, and set both tolerances to `1e-6`. Soft bounds use
`q >= 0.1`; box bounds use `0.1 <= q <= 10`.

| Mode and grid | Edges | Baseline (ms) | Basin (ms) | Time change |
| --- | ---: | ---: | ---: | ---: |
| Soft, 10 x 10 | 180 | 21.069 | 20.783 | -1.4% |
| Soft, 32 x 32 | 1,984 | 194.317 | 143.492 | -26.2% |
| Box, 48 x 48 | 4,512 | 13.042 | 12.500 | -4.2% |

Each cell below lists **baseline / Basin**. Loss and gradient are recomputed
at the returned parameters, outside the timed run. Evaluations count fused
objective/gradient calls made during optimization, including line-search
trials. The residual is the Euclidean gradient norm for soft bounds and the
projected-gradient infinity norm for box bounds.

| Case | Final loss | Gradient residual | Iterations | Evaluations | Stop reason (both) |
| --- | ---: | ---: | ---: | ---: | --- |
| Soft, 10 x 10 | 12.98921036 / 15.37482678 | 0.09670735 / 0.01876123 | 200 / 200 | 720 / 740 | Iteration limit |
| Soft, 32 x 32 | 19,788.11560 / 2,568.73200 | 4.544078 / 11.88551 | 200 / 200 | 633 / 471 | Iteration limit |
| Box, 48 x 48 | 26,427,795.5224 / 26,427,820.9137 | 9.833693 / 9.819662 | 12 / 11 | 14 / 12 | Relative cost tolerance |

The soft-bound implementations take different paths. At the iteration limit,
Basin has a higher loss on the smaller grid and a lower loss on the larger
grid; neither run meets its stopping tolerances. Faster execution alone does
not establish a better solution. On the bounded grid, both satisfy the
relative cost criterion, with final losses differing by about `9.6e-7`
relative to the baseline. Their projected-gradient residuals remain well
above `1e-6`, so this is cost convergence, not gradient convergence.

## Bounded reference problems

The reference suite uses the existing, independently generated L-BFGS-B 3.0
[fixtures and provenance](../lbfgsb/tests/reference/PROVENANCE.md). It does
not replace or regenerate them. The cases cover the three upstream drivers,
mixed finite/one-sided/unbounded variables, a fixed variable, and history
rollover. Every correctness run checks feasibility at each evaluation as
well as the final objective and projected residual. Those feasibility
assertions are disabled for both timed backends.

Driver 1 uses history length 5, projected-gradient tolerance `1e-5`, and
relative cost tolerance `1e7 * f64::EPSILON`. Drivers 2 and 3 use history
lengths 5 and 10 and the upstream custom stopping rule. The mixed quadratic
uses history length 5 and projected-gradient tolerance `1e-12`; the audit
case uses history length 2, projected-gradient tolerance `1e-10`, and
relative cost tolerance `1e7 * f64::EPSILON`.

| Case | Variables | Baseline (us) | Basin (us) | Basin / baseline | Iterations / evaluations (both) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Driver 1 | 25 | 34.267 | 33.782 | 0.99 | 23 / 28 |
| Driver 2 | 25 | 71.489 | 69.755 | 0.98 | 46 / 53 |
| Driver 3 | 1,000 | 2,864.619 | 2,546.342 | 0.89 | 49 / 58 |
| Mixed bounds | 4 | 0.795 | 0.971 | 1.22 | 2 / 3 |
| Fixed variable/history rollover | 8 | 3.523 | 4.118 | 1.17 | 8 / 11 |

The numerical results below again list **baseline / Basin**. Both revisions
pass the reference checks. Exact iteration trajectories are not required.

| Case | Final loss | Projected-gradient infinity norm |
| --- | ---: | ---: |
| Driver 1 | 1.083490083505e-9 / 1.083490083469e-9 | 1.720523e-4 / 1.720523e-4 |
| Driver 2 | 5.807023130099e-15 / 5.807023130104e-15 | 6.619508e-11 / 6.619686e-11 |
| Driver 3 | 5.352273181607e-22 / 5.349971425016e-22 | 9.745380e-11 / 9.730987e-11 |
| Mixed bounds | 8.184431891668e-30 / 8.135128085092e-30 | 4.440892e-15 / 3.552714e-15 |
| Fixed variable/history rollover | 128.7 / 128.7 | 0 / 0 |

Basin 1.13.0 is within a few percent of the baseline on the 25-variable
drivers and takes 11.1% less time on the 1,000-variable driver with the
selected `Vec<f64>` backend. The four- and eight-variable cases still take
22.1% and 16.9% more time, respectively, adding about 0.18 and 0.60 us per
solve. The matching iteration/evaluation counts make these reference cases
an optimizer overhead comparison; they do not support a general speedup
claim. The full Theseus grid measurements also include the cost of the
forward and adjoint solves and changes in the number of evaluations.

## Validation

With Basin 1.13.0, `cargo test --locked --workspace --release` passes all
192 tests; six manual benchmarks are ignored by that command. The comparison
script separately runs its grid and reference benchmarks and verifies the
bounded reference solutions on both revisions. The preserved standalone
solver also passes all 49 tests with
`cargo test --locked -p ariadne-lbfgsb --features faer-backend`.

The workspace suite includes progress reporting, callback and atomic-flag
cancellation, iteration limits, failed line searches, variable supports,
and FFI round trips. The distributed Basin license and attribution files
match the published 1.13.0 crate byte for byte.

Windows and macOS distribution checks were not run on this NixOS host.
`dotnet build Ariadne.csproj -c Release` fails here because the C# compiler
cannot resolve `ToolStripDropDown` in the Grasshopper components.

# Null-space paper benchmarks

`tests/bench_nullspace.rs` is an ignored, dependency-free harness for the
Pellegrino--Calladine implementation. It compares the dense referee with the
sparse production path and keeps every method/fixture cell in a fresh process.
Consequently, an out-of-memory dense SVD is recorded as `process-failed` and
does not invalidate later projector cells.

Run the default release sweep from the repository root:

```powershell
cargo test --release -p theseus --test bench_nullspace `
  bench_nullspace_fresh_processes -- --ignored --nocapture
```

The default fixtures are the reduced 4-bar, a prestressed three-spoke triangle,
a planar grid under normal loads, and corner-anchored warped grids of sizes 4,
8, and 12. Override the matrix without changing the harness:

```powershell
$env:THESEUS_NULLSPACE_FIXTURES = 'four-bar,prestressed-triangle,planar-z-grid'
$env:THESEUS_NULLSPACE_GRIDS = '8,16,24,32,48,64,96,128'
$env:THESEUS_NULLSPACE_METHODS = 'dense,projector,qr,angles,saddle,gram'
$env:THESEUS_NULLSPACE_LAMBDAS = '0,1e-12,1e-10,1e-8,1e-6'
$env:THESEUS_NULLSPACE_MAX_MODES = '32'
cargo test --release -p theseus --test bench_nullspace `
  bench_nullspace_fresh_processes -- --ignored --nocapture
```

Increase the grid list until the dense child is killed or cannot allocate.
Do not run an intentionally memory-exhausting sweep alongside Grasshopper or
other unsaved work. The projector currently obtains exact nullities from
projector traces, requiring one pseudoinverse action per row and column; this
is a correctness-first implementation, so its runtime can become the limiting
factor before its sparse memory scaling does.

## Methods and validity

- `dense` densifies the force-form equilibrium matrix \(A\) and computes a
  full `faer-svd`. This is the naive paper baseline and small-fixture referee.
- `projector` applies \(I-A^+A\) and \(I-AA^+\) using the production LSQR
  Moore--Penrose action, followed by SVD only on small probe panels. It never
  densifies \(A\).
- `qr` uses faer sparse QR's built-in COLAMD ordering. It is available only
  when rows are at least columns. The basic least-squares solution is always
  residual-checked; rank and Calladine counts are deliberately never
  published from its diagonal. Rank-deficient non-finite solves are reported
  as `unsupported-rank-deficient`.
- `angles` runs dense SVD and projector together in a separate comparison
  process. It reports the largest principal angle from the singular values of
  \(N_\mathrm{svd}^T N_\mathrm{projector}\), and likewise for mechanisms.
  Its time and memory are not attributed to either primary method.
- `saddle` measures the augmented particular solve. At \(\lambda=0\), this is
  the Moore--Penrose particular with the LSQR fallback; positive values are
  Tikhonov solves.
- `gram` explicitly forms \(A^TA+\lambda I\), demonstrating normal-equation
  fill and condition squaring. A singular zero-\(\lambda\) cell is expected to
  report unsupported. The lambda sweep applies only to these particular
  solves, never to kernel rank.

All methods use the same assembled force-form \(A\) within a fixture.

## Output

Rows beginning with `RESULT` are tab-separated for direct import. The columns
are:

- `wall_ms`: complete child method time after fixture assembly;
- `peak_mib`: process peak working set (`GetProcessMemoryInfo` on Windows,
  `VmHWM` on Linux);
- `rank`, `s`, `m_raw`: SVD/projector counts only;
- `calladine`: integer residual
  \(s-m_\mathrm{raw}-(n_e-n_\mathrm{eq})\), which must be zero;
- `AN_F`, `ATPhi_F`: measured Frobenius kernel residuals
  \(\lVert AN\rVert_F\) and \(\lVert A^T\Phi\rVert_F\);
- `load_residual`: \(\lVert At^+-p\rVert_2/\max(\lVert p\rVert_2,\epsilon)\);
- `self_angle_deg`, `mechanism_angle_deg`: largest principal angles, present
  only in `angles` cells where both complete returned subspaces fit the mode
  cap;
- `status` and `note`: unsupported shapes/rank, allocation/process failures,
  rigid-body stripping, and returned-mode truncation.

Peak working set includes the Rust test process baseline and loaded libraries.
Compare cells from the same build and machine; do not subtract two peaks or
interpret small differences as allocator-exact live memory.

## Publication protocol

Use a release build, close unrelated high-memory applications, record OS, CPU,
RAM, Rust/Cargo versions, commit, thread settings, grid list, tolerance, mode
cap, and the complete `RESULT` stream. Repeat successful cells in independent
sweeps. Treat a dense process failure as a memory-wall observation only after
confirming the preceding grid succeeds; record the failing grid and available
physical memory. Never infer QR rank from this harness.

This document intentionally contains no paper baseline yet. A focused smoke
run validates compilation and process isolation, but it is not a controlled
end-to-end measurement suitable for publication.
