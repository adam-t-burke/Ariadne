# Theseus benchmarks

1. [Scaling the FDM solve to 100k edges](#scaling-the-fdm-solve-to-100k-edges)
2. [Basin optimizer comparison](#basin-optimizer-comparison)
3. [Null-space paper benchmarks](#null-space-paper-benchmarks)

# Scaling the FDM solve to 100k edges

Everything in a Theseus iteration other than the L-BFGS-B update is one fused
objective/gradient evaluation: assemble `A(q) = Cnᵀ diag(q) Cn`, factor it,
solve for the free-node positions, evaluate the objectives, solve the adjoint
system with the same factor, and accumulate the gradient. This section
profiles that evaluation on square cable-net grids up to 224×224 (99,904
edges, 50,172 free nodes), records the optimisations made on the strength of
the profile, and reports the matrix-free experiment.

## Reproduce

```sh
# Per-phase profile of one evaluation (grid sizes as arguments).
RAYON_NUM_THREADS=1 cargo run --release -p theseus --example profile_phases -- 72 160 224

# Full box-constrained L-BFGS-B solves with a fixed iteration budget.
RAYON_NUM_THREADS=1 cargo test --release -p theseus --test bench_scale -- --ignored --nocapture
THESEUS_SCALE_GRIDS=72,160,224,320 THESEUS_SCALE_ITERS=40 cargo test --release -p theseus --test bench_scale -- --ignored --nocapture

# Matrix-free PCG versus the direct solve.
RAYON_NUM_THREADS=1 cargo run --release -p theseus --example matrix_free_pcg -- 224
```

The fixture (`tests/support/grid.rs`) is an `n × n` grid with four corner
supports, unit downward loads and a `TargetXYZ` objective. `bench_scale` uses
the *recoverable* variant, whose target is the equilibrium shape of a smooth
interior force-density field, so the optimum is reachable inside the bounds
and the line search behaves as it does on real problems (about 1.2
evaluations per accepted iteration). With the unreachable flat target the
solver stalls on the upper bound and burns ~3 evaluations per iteration, which
says nothing about the linear algebra.

Environment: 4 vCPU Linux VM, Rust 1.89, faer-sparse 0.17.1, Basin 1.13.0.
Timings are medians of 3–5 repetitions and vary by roughly ±5% between runs.

## Where an evaluation spends its time

Milliseconds per phase of one warm evaluation (symbolic analysis already
done), single-threaded, before this work:

| grid | edges | free nodes | assemble A | assemble b | numeric factor | solve (3 rhs) | geometry | loss | explicit ∇ | adjoint solve | implicit ∇ | total |
|-----:|------:|-----------:|-----------:|-----------:|---------------:|--------------:|---------:|-----:|-----------:|--------------:|-----------:|------:|
|  72 | 10,224 | 5,180 | 0.06 | 0.12 | 2.06 | 0.36 | 0.34 | 0.02 | 0.05 | 0.36 | 0.07 | 3.5 |
| 160 | 50,880 | 25,596 | 0.43 | 0.58 | 25.7 | 6.5 | 1.48 | 0.10 | 0.33 | 6.7 | 0.39 | 42.4 |
| 224 | 99,904 | 50,172 | 1.37 | 1.27 | 51.1 | 12.4 | 3.3 | 0.22 | 0.79 | 11.7 | 0.62 | 83.1 |

With the default rayon pool (4 threads) the same evaluation took 94 ms at
224: the numeric factorization was not faster (`Parallelism::Rayon(0)` gave
no speedup at any size and was 6–15% slower than `Parallelism::None`), and
geometry went from 3.3 ms to 12.9 ms because the reaction fold allocated an
`nn × 3` buffer per rayon work split and reduced them.

After the changes below:

| grid | assemble A | assemble b | numeric factor | solve (3 rhs) | geometry | loss | explicit ∇ | adjoint solve | implicit ∇ | total |
|-----:|-----------:|-----------:|---------------:|--------------:|---------:|-----:|-----------:|--------------:|-----------:|------:|
|  72 | 0.04 | 0.01 | 1.56 | 0.20 | 0.04 | 0.02 | 0.05 | 0.20 | 0.06 | 2.2 |
| 160 | 0.19 | 0.04 | 18.7 | 1.86 | 0.20 | 0.04 | 0.32 | 1.87 | 0.29 | 23.6 |
| 224 | 0.38 | 0.11 | 38.9 | 3.87 | 0.40 | 0.10 | 0.74 | 4.05 | 0.55 | 49.3 |

The evaluation at 100k edges is 1.7× faster single-threaded (83 → 49 ms) and
1.9× faster with the default pool (94 → 51 ms), and is now bit-identical
across thread counts. The numeric Cholesky factorization is ~79% of what is
left.

## What changed

* **Triangular solves** (`src/factor_solve.rs`). faer's generic
  `solve_in_place_with_conj` was 12 ms per 3-column solve at 100k edges,
  roughly a quarter of the factorization cost, because it goes through
  entity-generic dense kernels per supernode and a column-major right-hand
  side. The new routines walk faer's public factor layout directly
  (simplicial CSC with the diagonal first; supernodal column-major blocks;
  LLᵀ and LDLᵀ) with the three coordinates of a row stored contiguously, so a
  factor column is streamed once and the inner loops vectorise over the three
  right-hand sides. 12.4 ms → 3.9 ms. Unit tests check the result against
  faer's own solver for every layout/kind combination.
* **Scratch buffers**: the faer stack buffer was reallocated on every
  factorization and every solve; it is now grown only when a request exceeds
  the current size.
* **Parallelism**: `Parallelism::None` for the numeric factorization and
  solves. On these factors (~29 nonzeros per column, supernodes of a few
  columns) faer's Rayon path has nothing to parallelise and only adds
  scheduling cost.
* **Assembly**: `assemble_a` gathers each nonzero from a flat CSR map of its
  contributing edges instead of zero-filling and scattering through nested
  `Vec`s; `assemble_rhs` accumulates only the edges with one fixed endpoint
  instead of two sparse-dense products over all edges. The `ne × 3` scratch
  arrays and the cached `Cn`, `Cnᵀ`, `Cf` copies were removed.
* **Geometry**: one sequential pass computes lengths, forces and reactions.
  The rayon fold is gone (3.3 ms → 0.4 ms single-threaded; 12.9 ms → 0.44 ms
  with 4 threads), and the reaction sums are now deterministic.

## Full solves

40 L-BFGS-B iterations on the recoverable grid, box bounds `[0.1, 10]`,
tolerances disabled so every run does the full budget. `evals` is the number
of fused evaluations; `non-eval` is the total minus `evals × eval`, i.e.
setup plus Basin's own bookkeeping (parameter copies, two-loop recursion,
line search).

Before (merge of the Basin PR onto the warm-start branch), 1 thread:

| grid | edges | setup ms | eval ms | iters | evals | total ms | ms/iter | non-eval ms | peak RSS MB |
|-----:|------:|---------:|--------:|------:|------:|---------:|--------:|------------:|------------:|
|  72 | 10,224 | 8.0 | 3.68 | 40 | 50 | 181 | 4.5 | −3 | 17 |
| 160 | 50,880 | 55.1 | 36.7 | 40 | 56 | 2,156 | 53.9 | 102 | 80 |
| 224 | 99,904 | 114.3 | 74.3 | 40 | 55 | 4,310 | 107.8 | 225 | 147 |

Before, default rayon pool: 204 / 2,253 / 5,034 ms, peak RSS 22 / 110 / 210 MB,
and the final loss differed from the single-threaded run (non-deterministic
reaction sums).

After, 1 thread:

| grid | edges | setup ms | eval ms | iters | evals | total ms | ms/iter | non-eval ms | peak RSS MB |
|-----:|------:|---------:|--------:|------:|------:|---------:|--------:|------------:|------------:|
|  72 | 10,224 | 7.3 | 2.51 | 40 | 50 | 145 | 3.6 | 19 | 14 |
| 160 | 50,880 | 46.8 | 25.5 | 40 | 53 | 1,456 | 36.4 | 106 | 61 |
| 224 | 99,904 | 100.5 | 55.4 | 40 | 53 | 3,195 | 79.9 | 261 | 110 |

After, default rayon pool: 189 / 1,414 / 3,097 ms, peak RSS 14 / 59 / 111 MB,
final loss identical to the single-threaded run.

At 100k edges a full iteration is now ~80 ms (from 108 ms single-threaded and
126 ms with the default pool), so a 200-iteration solve takes about 16 s in
~110 MB. Basin's bookkeeping is ~6 ms per iteration at that size (`Vec`
clones of the 100k-element parameter and gradient vectors); it is 8% of the
iteration and the cost of its owned-`Vec` API rather than anything in
Theseus.

## Matrix-free solves

`examples/matrix_free_pcg.rs` applies `A(q)` without assembling it — one pass
over the edges — and solves the forward system with conjugate gradients on
the three coordinate columns at once. Two preconditioners were tried: Jacobi
(pure matrix-free) and the existing Cholesky factor of `A(q_prev)` from the
previous evaluation ("frozen factor"). Each is started from zero and from the
previous solution `x_prev`, for relative force-density steps of 0, 1, 5 and
20%. 224×224 grid, `tol` is the relative residual, single thread:

| method | q step | iterations | time | direct refactor + solve |
|--------|-------:|-----------:|-----:|------------------------:|
| Jacobi-CG, cold, tol 1e-8 | any | 1,190–1,220 | 1.15–1.17 s | 44–47 ms |
| Jacobi-CG, warm, tol 1e-8 | 1% | 1,100 | 1.07 s | |
| Jacobi-CG, warm, tol 1e-8 | 20% | 1,238 | 1.20 s | |
| frozen factor, warm, tol 1e-8 | 1% | 3 | 30 ms | |
| frozen factor, warm, tol 1e-8 | 5% | 5 | 29 ms | |
| frozen factor, warm, tol 1e-8 | 20% | 7 | 39 ms | |

One matrix-free `A·x` for three right-hand sides costs 0.4 ms, so the
operator is not the problem; the iteration count is. `A` is a weighted grid
Laplacian whose condition number grows with the number of nodes, and Jacobi
does nothing about that: ~1,200 iterations regardless of the start vector
(the warm start only helps when `q` has not changed). A pure matrix-free path
would need a multilevel preconditioner (AMG or a geometric hierarchy of the
network) to be competitive, which is a project in itself.

Preconditioning with the previous factor converges in 3–7 iterations, and
with the new triangular solves each iteration costs ~5.5 ms, so the forward
solve is 30–39 ms against 44–47 ms for a refactor. But an evaluation also
needs the adjoint solve, which with the direct method reuses the fresh factor
for 4 ms and with PCG needs its own 3–7 iterations. Per evaluation that is
34–76 ms against 47 ms: a gain only when consecutive force densities differ
by a percent or so, a loss for the larger steps typical early in a solve, and
the gradient becomes inexact to the CG tolerance. The variant is therefore
left as an experiment. Making it pay would need an adaptive policy (attempt
`k ≤ 3` PCG iterations, refactor on failure) and tolerance handling in the
optimizer for inexact gradients, for a ceiling of roughly 1.3× on the
favourable iterations.

## Remaining floor

The numeric factorization is 39 ms of the 49 ms evaluation at 100k edges.
faer-sparse 0.17 only offers AMD ordering (fill is already good at ~29
nonzeros per factor column) and its supernodal kernels do not parallelise on
supernodes this small. Forcing the simplicial factorization cuts the factor's
fill by 25% and halves the solve, but the numeric factorization becomes 2.5×
slower, so the automatic choice stands. Further gains would come from a
nested-dissection ordering (not exposed by faer 0.17) or a different sparse
Cholesky implementation, not from anything around the factorization.

## Toward 1M edges

Same fixture at 320, 448 and 708 (204k, 401k and 1,001k edges), single
thread, after the changes above. Milliseconds per phase of one evaluation:

| grid | edges | free nodes | assemble | numeric factor | solve | geometry | loss + ∇ | adjoint solve | total |
|-----:|------:|-----------:|---------:|---------------:|------:|---------:|---------:|--------------:|------:|
| 320 | 204,160 | 102,396 | 1.3 | 87 | 9.0 | 0.8 | 4.3 | 9.3 | 113 |
| 448 | 400,512 | 200,700 | 3.7 | 200 | 20.7 | 2.0 | 8.3 | 19.4 | 256 |
| 708 | 1,001,112 | 501,260 | 9.4 | 661 | 85.8 | 5.8 | 24.0 | 87.0 | 877 |

Full solves (10 iterations): 1.95 s / 4.47 s / 14.8 s, i.e. 0.19 / 0.45 /
1.48 s per iteration at 1.5 evaluations per iteration; setup 0.21 / 0.47 /
1.51 s; peak RSS 214 / 438 / 1,123 MB. Basin's bookkeeping is ~20 ms per
iteration at 1M and irrelevant.

Factor properties (faer 0.17, AMD, automatic supernodal selection):

| grid | nnz(A) | true nnz(L) | per column | supernodal nnz(L) | supernodes | symbolic | numeric f64 | numeric f64 Rayon | numeric f32 | solve (3 rhs) |
|-----:|-------:|------------:|-----------:|------------------:|-----------:|---------:|------------:|------------------:|------------:|--------------:|
| 224 | 250k | 1.45M | 28.8 | 1.89M | 25k | 12 ms | 34 ms | 41 ms | 29 ms | 4.5 ms |
| 320 | 511k | 3.36M | 32.9 | 4.35M | 51k | 26 ms | 92 ms | 92 ms | 71 ms | 10.0 ms |
| 448 | 1.00M | 7.45M | 37.1 | 9.41M | 101k | 55 ms | 186 ms | 195 ms | 140 ms | 21.7 ms |
| 708 | 2.50M | 21.9M | 43.8 | 27.2M | 253k | 140 ms | 605 ms | 587 ms | 433 ms | 75.9 ms |

Things checked and their outcome:

* **Ordering.** Fill grows like `log n` (29 → 44 nonzeros per column from
  50k to 500k nodes), so AMD is not the problem. faer 0.17 has no
  nested-dissection ordering and its symbolic analysis does not accept a
  user permutation, so the flop count is what it is.
* **faer's parallel factorization.** `Parallelism::Rayon(0)` is within noise
  of `None` at 1M (587 vs 605 ms) and slower below that. Supernodes average
  two columns on these planar graphs, so there is nothing for its dense
  kernels to parallelise.
* **Supernode relaxation.** More aggressive amalgamation than faer's
  (CHOLMOD-default) cutoffs barely changes the supernode count (253k → 250k
  at 1M) and buys 2–10% on the factorization for 10–28% more fill, which the
  two solves pay back. Disabling relaxation is 20% slower. Leave the default.
* **Single precision.** An `f32` factorization is 1.4× faster (compute
  bound, not bandwidth bound), and would only be usable as a preconditioner
  followed by `f64` refinement; not worth the machinery on its own.
* **Frozen-factor PCG at scale.** The factor/solve ratio stays ~8 from 100k
  to 1M edges, so the break-even analysis from the 224 grid still holds:
  3–7 preconditioned iterations per solve, two solves per evaluation, no
  gain unless the preconditioner solve gets cheaper. The one thing that
  changes at 1M is that the solve is now memory bound (27M values streamed
  twice in 76 ms ≈ 5.7 GB/s), so storing the frozen factor in `f32` would
  roughly halve its cost and make the hybrid pay when consecutive force
  densities differ by a few percent. Estimated, not measured.
* **Setup.** `FdmCache::new` is 0.52 s at 1M: 0.10 s transposing `Cn`,
  0.21 s forming the pattern of `CnᵀCn` through sorted triplets, the rest
  in the `q → nonzero` map and per-node incidence lists. The pattern can be
  built directly from the edge list in one pass. Symbolic analysis and the
  first factorization add 0.75 s. Both are paid once per `optimize()` call,
  which for interactive re-solves on a fixed topology is once per re-solve.
* **Memory.** 1.1 GB peak at 1M edges: the factor is ~220 MB, the rest is
  construction transients (triplet buffers, per-edge `Vec`s) and
  `Vec<Vec<usize>>` incidence lists that could be CSR.

Conclusion: with faer 0.17's sparse Cholesky, 1M edges costs ~0.9 s per
evaluation and ~1.1 GB, and about 70% of the time is the numeric
factorization, which nothing around it can reduce. Getting substantially
below that requires either a different direct solver (nested dissection and
a parallel supernodal kernel) or a multilevel preconditioner that makes the
matrix-free path converge in tens rather than a thousand iterations. The
latter is the only option whose cost per evaluation and memory both scale
linearly.

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

# WS-B: parallel O(ne) loops

The O(ne)/O(nn) loops of the fused evaluation (`assemble_a`, `assemble_rhs`,
`compute_geometry`, the objective reductions, the explicit node-position
gradients and `accumulate_implicit_gradients`) now run as chunked rayon loops
over fixed 4,096-element pieces, with per-node work expressed as gathers over
the CSR adjacency `FdmCache::adjacency` (incident edges in ascending edge
order) instead of edge scatters. Chunk boundaries depend only on the data
length and reduction partials are combined in chunk order, so every result is
bitwise identical for any pool size (`tests/graph_loops_determinism.rs`
checks one evaluation and a full 8-iteration solve on the 160 grid under
pools of 1, 2 and 4 threads). Loops under 16,384 elements, and every loop on
a single-thread pool, run on the calling thread with the same chunking.

## Reproduce

```sh
RAYON_NUM_THREADS=1 cargo run --release -p theseus --example profile_phases -- 224 708
RAYON_NUM_THREADS=4 cargo run --release -p theseus --example profile_phases -- 224 708
```

## Per-phase times

Milliseconds per phase of one warm evaluation; `before` is the parent
commit built from the same tree, `after` this branch, run interleaved on the
same 4 vCPU VM and reported as the minimum over 8 runs (each run averages
3 repetitions). One core was occupied by an unrelated process for the whole
session, so the 4-thread rows effectively had three cores and run-to-run
variation was ±10% on the sub-millisecond phases; the factorization and the
two triangular solves (`solve`, `adjoint`) are unchanged by this work and
serve as the noise reference.

`RAYON_NUM_THREADS=1`:

| grid | edges | | assemble A | assemble b | numeric factor | solve | geometry | loss | explicit ∇ | adjoint solve | implicit ∇ | total |
|-----:|------:|:--|-----------:|-----------:|---------------:|------:|---------:|-----:|-----------:|--------------:|-----------:|------:|
| 224 | 99,904 | before | 0.42 | 0.08 | 38.1 | 4.05 | 0.39 | 0.14 | 0.41 | 4.03 | 0.57 | 48.6 |
| 224 | 99,904 | after | 0.41 | 0.07 | 38.2 | 3.92 | 0.40 | 0.13 | 0.40 | 3.91 | 0.37 | 48.1 |
| 708 | 1,001,112 | before | 8.08 | 1.62 | 668 | 93.6 | 6.00 | 2.23 | 7.44 | 95.4 | 6.72 | 900 |
| 708 | 1,001,112 | after | 8.23 | 1.49 | 677 | 96.1 | 6.09 | 2.31 | 7.55 | 98.5 | 5.91 | 908 |

`RAYON_NUM_THREADS=4`:

| grid | edges | | assemble A | assemble b | numeric factor | solve | geometry | loss | explicit ∇ | adjoint solve | implicit ∇ | total |
|-----:|------:|:--|-----------:|-----------:|---------------:|------:|---------:|-----:|-----------:|--------------:|-----------:|------:|
| 224 | 99,904 | before | 0.70 | 0.14 | 42.1 | 4.71 | 0.51 | 0.26 | 0.66 | 4.39 | 0.61 | 55.2 |
| 224 | 99,904 | after | 0.41 | 0.06 | 42.1 | 5.06 | 0.52 | 0.10 | 0.49 | 4.31 | 0.31 | 54.3 |
| 708 | 1,001,112 | before | 8.05 | 1.64 | 673 | 101.3 | 6.11 | 2.29 | 7.44 | 100.8 | 6.87 | 913 |
| 708 | 1,001,112 | after | 2.98 | 0.54 | 669 | 97.4 | 2.72 | 0.81 | 4.09 | 95.2 | 2.41 | 881 |

Sum of the six rewritten phases at 1M edges (assemble A, assemble b,
geometry, loss, explicit ∇, implicit ∇):

| threads | before | after | speed-up |
|--------:|-------:|------:|---------:|
| 1 | 32.1 ms | 31.6 ms | 1.02× |
| 4 | 32.4 ms | 13.6 ms | 2.4× |

A quieter run of the same protocol (load average 0.95, min over 6) gave
31.8 → 11.6 ms at 4 threads (2.7×) with the same 1-thread parity (31.8 →
31.4 ms); the loops are memory-bound gathers, so with three effective cores
2.4–2.7× is what the machine allows. At 1 thread no phase moved by more than
0.1 ms except the implicit gradient, which is 12% faster because the fixed
node contributions are gathered per node instead of scanning every edge for
a fixed endpoint. The triangular solves (~190 ms at 1M) are sequential and
outside this workstream, so the whole evaluation moves only from 913 to
881 ms; replacing them and the factorization with the AMG-preconditioned
iterative path (WS-C/WS-D) is what this infrastructure is for.

## Numerical changes

* `assemble_a`, `assemble_rhs`, `compute_geometry`, the explicit node
  gradients and `accumulate_implicit_gradients` reproduce the previous
  sequential loops **bitwise**: each output entry is a sum over the same
  terms in the same order (incident edges ascending, boundary edges sorted
  by free row then edge), and the parallel and sequential forms of each loop
  are identical.
* The objective loss reductions (`objectives.rs`) are the one place the
  order changed: sums are formed sequentially inside 4,096-element chunks
  and the chunk partials added in order, so losses over more than 4,096
  entries may differ from the old values in the last bits (tested at 1e-12
  relative against plain sequential sums; `total_loss` still adds the
  objectives in order). Sums of at most 4,096 entries are unchanged.
