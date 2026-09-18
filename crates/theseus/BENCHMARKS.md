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

## Phase 0: aggregation AMG prototype

Phase 0 of `AMG_PLAN.md` / workstream WS-A of `ITERATIVE_SOLVER_PROGRAM.md`:
a self-contained, single-threaded prototype of the matrix-free multilevel
solve, used to fix the algorithm choices before any production code. The
code is `examples/amg_prototype.rs`; it changes nothing under `src/`.

```sh
export RUSTUP_TOOLCHAIN=1.89.0
# Correctness: coarse operator = PᵀAP on a random 200-node graph, V-cycle
# preconditioner symmetric to 1e-10 and positive definite, solve matches faer.
cargo run -p theseus --release --example amg_prototype -- --self-test
# One sweep (plain aggregation), all sizes of one fixture:
for s in 64 128 224 448 708; do
  RAYON_NUM_THREADS=1 cargo run -p theseus --release --example amg_prototype -- \
    --fixture grid --size $s --cycle v,k --degree 1,2,3 --passes 1,2 --warm 2 --reps 3
done
# Smoothed aggregation, the recommended configuration:
RAYON_NUM_THREADS=1 cargo run -p theseus --release --example amg_prototype -- \
  --fixture irregular --size 708 --cycle v --degree 2 --passes 3 --alpha 10 --smoothed 1
```

What the prototype contains (all in the example):

* the free-node weighted graph Laplacian `A(q) x = anchor ⊙ x + Σ_e q_e
  (x_u − x_v)` in node-gather form (CSR incidence, each node sums its
  incident edges), on three right-hand sides at once;
* pairwise heavy-edge aggregation: strength `w_uv / sqrt(d_u d_v)`, nodes
  visited by decreasing degree (then index) and matched with their strongest
  unmatched neighbour, leftovers joining their strongest neighbour's pair
  (at most three members) or staying single; 1–3 matching passes per level
  (aggregates of ~2/4/8 nodes); coarsening stops below `coarsest` nodes or
  when a level keeps more than 70% of the nodes of the previous one; coarse
  operator = quotient graph (`w_c = Σ q`, `anchor_U = Σ anchor_u`), so every
  level reuses the edge kernel and a `q` change is an `O(nnz)` weight pass
  per level;
* Chebyshev smoother of degree 1–3 with Jacobi scaling on `[λ_max/α, λ_max]`,
  `λ_max` from 10 power iterations on `D⁻¹A` per level at setup, with a 10%
  margin, kept across the small `q` updates measured here;
* V-cycle, and Notay's K-cycle (at each coarse level up to two flexible-CG
  steps preconditioned by the next-coarser cycle, the second skipped when
  the first already reduced the residual by 4×);
* coarsest level solved with faer's sparse Cholesky through
  `theseus::types::Factorization` (≤ 2,000 nodes, well under 1 ms per solve);
* PCG (V-cycle) and flexible CG FCG(1) (K-cycle) on the block of three
  right-hand sides, per-column convergence to relative residual `1e-8`;
* the fallback of `AMG_PLAN.md`: smoothed aggregation, `P = (I − ω D⁻¹A) P₀`
  with `ω = 4/(3 λ_max)` and `A_c = Pᵀ A P`, whose coarse operators are
  general sparse matrices and whose `q` update is a sparse triple product
  per level rather than a weight pass.

Fixtures, each with a smooth non-constant force-density field
`q = exp(ln(ratio) · (½ + ½ sin 2πx cos 2πy))` (ratio 10 unless stated) and
unit downward loads:

* **grid**: the `n × n` cable net of `tests/support/grid.rs` (four corner
  supports, degree 4);
* **irregular**: an `m × m` lattice jittered by up to 0.4 of the spacing,
  every point joined to its 5–7 nearest neighbours, edge set symmetrised,
  boundary ring fixed (degree 5–9);
* **dome**: a spoke wheel / cable dome — hub, concentric rings whose node
  count doubles when the circumferential spacing exceeds 1.5× the radial
  one, radial spokes and diagonal bracing, outer ring fixed (hub degree
  12, degree 8 on doubling rings, 6 elsewhere).

Sizes are matched by edge count to the grids 64, 128, 224, 448 and 708
(8k, 33k, 100k, 400k and 1M edges). The direct reference is the library
path (`FdmCache::new` + `Factorization::update` + `solve_slices` on the same
three right-hand sides): one numeric refactorization plus one 3-rhs solve,
i.e. the per-`q` cost the optimizer pays today before the adjoint solve.

Method: single thread (`RAYON_NUM_THREADS=1`, and the prototype has no
parallelism), release build, Rust 1.89, faer-sparse 0.17.1. Iteration counts
are deterministic. Times are medians of 3 solves. "Cold" starts from `x = 0`;
"warm" perturbs every `q_e` by a uniform random ±2%, updates the hierarchy
(level weights and anchors and the coarsest factor; for smoothed aggregation
`P` and the coarse operators are rebuilt on the fixed aggregation) and
restarts from the cold solution against the new right-hand side; the update
is timed separately ("q-update"). The iterative solutions agree with
faer's to `5e-8` relative (max-norm) at the `1e-8` residual tolerance. The
4-vCPU VM was shared with other builds while the sweeps ran (load 2–6
during the first pass, 1–2 during the re-runs), which inflates the
memory-bound solves at 448 and 708 by up to 2×; every time reported is the
smallest median over the runs of that configuration, and the dedicated
1M-edge table was taken at load ≈ 1, where the grid-708 direct reference
reproduces the 757 ms of the section above. Iteration counts are unaffected
by load.

### Fixtures and direct reference

| fixture | grid-equivalent | edges | free nodes | direct refactor + 3-rhs solve (ms) |
|---|---:|---:|---:|---:|
| grid 64 | 64 | 8,064 | 4,092 | 1.4 |
| grid 128 | 128 | 32,512 | 16,380 | 12.6 |
| grid 224 | 224 | 99,904 | 50,172 | 46.0 |
| grid 448 | 448 | 400,512 | 200,700 | 226 |
| grid 708 | 708 | 1,001,112 | 501,260 | 756 |
| irregular 49 | 64 | 8,129 | 2,209 | 0.9 |
| irregular 98 | 128 | 32,376 | 9,216 | 6.1 |
| irregular 172 | 224 | 99,655 | 28,900 | 22.3 |
| irregular 345 | 448 | 400,434 | 117,649 | 118 |
| irregular 546 | 708 | 1,002,202 | 295,936 | 358 |
| dome 29 | 64 | 8,340 | 2,593 | 1.1 |
| dome 58 | 128 | 32,532 | 10,465 | 8.2 |
| dome 103 | 224 | 98,196 | 31,969 | 24.2 |
| dome 207 | 448 | 393,108 | 129,505 | 148 |
| dome 342 | 708 | 1,015,188 | 336,865 | 458 |

### Hierarchies

Levels, size of the coarse hierarchy relative to the fine edge count (at
708), and fine-operator-equivalent applications per outer iteration (all
levels' work converted to level-0 edge visits; the K-cycle's recursion
doubles it per level unless the early exit cuts it short). One matching
pass halves the node count per level and needs 8–9 levels at 500k nodes;
two passes quarter it (5 levels); three passes with smoothing divide by ~8
(4 levels).

| fixture | config | levels 64 → 708 | coarse edges / fine edges | A-applications per iteration 64 → 708 |
|---|---|---|---:|---|
| grid | plain V d1 p1 | 3 → 4 → 6 → 8 → 9 | 1.34 | 4.4 → 5.0 → 5.5 → 5.7 → 5.7 |
| grid | plain V d2 p2 | 2 → 3 → 4 → 5 → 5 | 0.45 | 5.1 → 6.4 → 6.7 → 6.8 → 6.8 |
| grid | plain K d2 p1 | 3 → 4 → 6 → 8 → 9 | 1.34 | 11.6 → 18.3 → 31.9 → 45.8 → 52.9 |
| grid | plain K d2 p2 | 2 → 3 → 4 → 5 → 5 | 0.45 | 5.1 → 8.3 → 10.1 → 10.9 → 11.0 |
| grid | smoothed V d2 p3 α10 | 2 → 2 → 3 → 4 → 4 | 0.45 | 5.1 → 5.1 → 6.5 → 6.8 → 6.8 |
| grid | smoothed K d2 p3 α10 | 2 → 2 → 3 → 4 → 4 | 0.45 | 5.1 → 5.1 → 7.4 → 9.2 → 8.9 |
| irregular | plain V d1 p1 | 2 → 4 → 5 → 7 → 8 | 0.85 | 3.1 → 4.5 → 4.6 → 4.7 → 4.7 |
| irregular | plain V d2 p2 | 2 → 3 → 3 → 4 → 5 | 0.25 | 5.1 → 5.9 → 5.9 → 6.0 → 6.1 |
| irregular | plain K d2 p1 | 2 → 4 → 5 → 7 → 8 | 0.85 | 5.1 → 12.0 → 15.3 → 21.2 → 23.1 |
| irregular | plain K d2 p2 | 2 → 3 → 3 → 4 → 5 | 0.25 | 5.1 → 7.1 → 7.1 → 8.0 → 8.3 |
| irregular | smoothed V d2 p3 α10 | 2 → 2 → 3 → 3 → 4 | 0.27 | 5.2 → 5.2 → 6.1 → 6.1 → 6.2 |
| irregular | smoothed K d2 p3 α10 | 2 → 2 → 3 → 3 → 4 | 0.27 | 5.2 → 5.2 → 6.4 → 6.4 → 6.6 |
| dome | plain V d1 p1 | 2 → 4 → 5 → 7 → 9 | 0.97 | 3.1 → 4.5 → 4.7 → 4.9 → 5.0 |
| dome | plain V d2 p2 | 2 → 3 → 3 → 4 → 5 | 0.32 | 5.1 → 6.0 → 6.0 → 6.2 → 6.3 |
| dome | plain K d2 p1 | 2 → 4 → 5 → 7 → 9 | 0.97 | 5.1 → 14.7 → 19.4 → 28.6 → 37.3 |
| dome | plain K d2 p2 | 2 → 3 → 3 → 4 → 5 | 0.32 | 5.1 → 7.4 → 7.5 → 8.7 → 9.3 |
| dome | smoothed V d2 p3 α10 | 2 → 2 → 3 → 3 → 4 | 0.34 | 5.2 → 5.2 → 6.2 → 6.2 → 6.5 |
| dome | smoothed K d2 p3 α10 | 2 → 2 → 3 → 3 → 4 | 0.34 | 5.2 → 5.2 → 6.6 → 6.6 → 7.2 |

### Iterations and wall time per configuration

`V`/`K` cycle, `d` Chebyshev degree, `p` matching passes per level, α = 30
unless shown, coarsest ≤ 2,000 nodes. Each cell is iterations to `1e-8`
for the x and z columns · solve time in ms (min over runs of the median of
3 solves). Setup is not included; it is 0.13–0.21 s (plain) and 0.46–0.64
s (smoothed) at 1M edges and paid once per topology.

**grid, cold start** — iterations (x column / z column) · solve ms

| method | config | 64 | 128 | 224 | 448 | 708 |
|---|---|---:|---:|---:|---:|---:|
| plain | V d1 p1 | 27/34 · 6.9 | 32/41 · 39.6 | 42/56 · 196 | 60/84 · 1,204 | 72/105 · 4,760 |
| plain | V d2 p1 | 20/24 · 8.3 | 24/31 · 46.0 | 34/45 · 246 | 51/69 · 1,554 | 62/88 · 5,702 |
| plain | V d3 p1 | 15/19 · 7.8 | 20/26 · 52.9 | 29/38 · 299 | 44/59 · 1,814 | 54/78 · 6,805 |
| plain | K d1 p1 | 26/32 · 10.8 | 28/37 · 79.7 | 29/39 · 594 | 30/42 · 2,576 | 30/43 · 8,957 |
| plain | K d2 p1 | 18/22 · 10.3 | 19/25 · 78.7 | 20/26 · 601 | 20/28 · 2,552 | 20/29 · 9,100 |
| plain | K d3 p1 | 13/16 · 10.2 | 14/18 · 72.9 | 15/20 · 429 | 15/21 · 2,559 | 15/22 · 8,555 |
| plain | V d1 p2 | 27/35 · 6.2 | 36/48 · 33.4 | 51/67 · 161 | 75/103 · 1,045 | 77/112 · 3,274 |
| plain | V d2 p2 | 23/29 · 6.6 | 33/41 · 43.1 | 45/58 · 210 | 65/88 · 1,348 | 68/100 · 4,304 |
| plain | V d3 p2 | 17/21 · 6.1 | 26/32 · 44.9 | 37/47 · 235 | 54/75 · 1,537 | 60/87 · 4,530 |
| plain | K d1 p2 | 27/35 · 5.9 | 31/41 · 40.9 | 33/45 · 162 | 35/49 · 808 | 35/50 · 2,668 |
| plain | K d2 p2 | 23/29 · 6.7 | 30/40 · 55.9 | 33/45 · 238 | 35/50 · 1,232 | 36/51 · 3,788 |
| plain | K d3 p2 | 17/21 · 6.4 | 20/26 · 47.8 | 21/29 · 201 | 23/33 · 1,045 | 23/33 · 2,806 |
| smoothed | V d2 p3 α10 | 12/15 · 4.4 | 13/16 · 17.5 | 14/19 · 74.0 | 16/22 · 353 | 15/22 · 920 |
| smoothed | V d3 p3 α10 | 10/12 · 3.6 | 11/14 · 19.7 | 12/16 · 80.6 | 13/18 · 381 | 13/18 · 982 |
| smoothed | K d2 p3 α10 | 12/15 · 3.7 | 13/16 · 19.3 | 14/19 · 84.4 | 15/21 · 440 | 15/21 · 1,144 |
| smoothed | K d3 p3 α10 | 10/12 · 3.7 | 11/14 · 22.9 | 12/16 · 86.7 | 12/17 · 393 | 12/17 · 1,009 |
| smoothed | V d2 p3 α30 | 16/20 · 4.6 | 17/22 · 25.0 | 19/25 · 97.0 | 21/29 · 475 | 21/29 · 1,302 |
| smoothed | V d3 p3 α30 | 12/14 · 5.6 | 12/16 · 22.7 | 14/18 · 89.3 | 15/20 · 429 | 15/21 · 1,192 |
| smoothed | K d2 p3 α30 | 16/20 · 4.8 | 17/22 · 25.5 | 19/26 · 125 | 21/29 · 639 | 21/29 · 1,754 |
| smoothed | K d3 p3 α30 | 12/14 · 4.7 | 12/16 · 23.2 | 13/18 · 104 | 14/18 · 490 | 14/19 · 1,468 |

**grid, warm start after a 2% `q` step** — iterations (x / z) · solve ms

| method | config | 64 | 128 | 224 | 448 | 708 |
|---|---|---:|---:|---:|---:|---:|
| plain | V d1 p1 | 19/26 · 5.3 | 21/32 · 32.1 | 27/41 · 141 | 40/65 · 928 | 41/77 · 3,282 |
| plain | V d2 p1 | 14/19 · 5.9 | 17/24 · 37.8 | 22/35 · 231 | 34/54 · 1,216 | 36/67 · 4,337 |
| plain | V d3 p1 | 11/15 · 6.2 | 14/20 · 41.4 | 19/29 · 300 | 29/47 · 1,437 | 31/60 · 5,347 |
| plain | K d1 p1 | 18/25 · 8.2 | 19/28 · 60.1 | 19/30 · 504 | 20/32 · 1,939 | 19/33 · 7,540 |
| plain | K d2 p1 | 13/18 · 9.2 | 13/19 · 61.0 | 13/20 · 314 | 13/21 · 1,915 | 13/22 · 6,727 |
| plain | K d3 p1 | 9/12 · 7.3 | 10/14 · 54.3 | 9/14 · 282 | 10/17 · 2,025 | 10/17 · 7,084 |
| plain | V d1 p2 | 19/27 · 4.7 | 25/37 · 26.8 | 32/50 · 121 | 49/80 · 814 | 45/83 · 2,386 |
| plain | V d2 p2 | 16/23 · 5.3 | 22/32 · 33.6 | 29/44 · 162 | 42/69 · 1,060 | 39/74 · 3,166 |
| plain | V d3 p2 | 12/16 · 5.0 | 17/26 · 39.0 | 23/36 · 176 | 36/59 · 1,238 | 34/66 · 3,626 |
| plain | K d1 p2 | 19/27 · 4.6 | 21/31 · 29.1 | 22/34 · 121 | 23/37 · 611 | 22/37 · 1,753 |
| plain | K d2 p2 | 16/23 · 5.5 | 20/30 · 42.0 | 21/32 · 171 | 23/38 · 918 | 21/38 · 2,463 |
| plain | K d3 p2 | 12/16 · 5.2 | 13/20 · 36.1 | 13/21 · 146 | 15/25 · 791 | 14/25 · 2,206 |
| smoothed | V d2 p3 α10 | 8/11 · 2.6 | 9/13 · 14.5 | 9/14 · 54.5 | 10/17 · 276 | 10/17 · 743 |
| smoothed | V d3 p3 α10 | 7/10 · 3.3 | 7/11 · 15.1 | 8/12 · 61.5 | 9/14 · 310 | 8/14 · 774 |
| smoothed | K d2 p3 α10 | 8/11 · 2.9 | 9/13 · 14.9 | 9/14 · 59.8 | 10/16 · 328 | 9/16 · 856 |
| smoothed | K d3 p3 α10 | 7/10 · 3.1 | 7/11 · 16.9 | 7/12 · 64.6 | 8/13 · 305 | 7/13 · 774 |
| smoothed | V d2 p3 α30 | 11/16 · 3.8 | 12/17 · 19.1 | 13/19 · 74.0 | 14/23 · 382 | 13/23 · 1,036 |
| smoothed | V d3 p3 α30 | 8/11 · 3.5 | 9/12 · 16.7 | 9/14 · 69.4 | 10/16 · 353 | 10/16 · 913 |
| smoothed | K d2 p3 α30 | 11/16 · 4.4 | 12/17 · 19.6 | 13/20 · 96.1 | 13/22 · 491 | 13/22 · 1,343 |
| smoothed | K d3 p3 α30 | 8/11 · 3.6 | 9/12 · 17.6 | 9/14 · 80.0 | 9/15 · 418 | 9/15 · 1,131 |

**irregular, cold start** — iterations (x column / z column) · solve ms

| method | config | 64 | 128 | 224 | 448 | 708 |
|---|---|---:|---:|---:|---:|---:|
| plain | V d1 p1 | 15/18 · 2.4 | 18/22 · 17.3 | 22/27 · 75.2 | 32/41 · 527 | 39/51 · 1,792 |
| plain | V d2 p1 | 12/14 · 2.7 | 16/19 · 24.4 | 20/24 · 108 | 30/36 · 758 | 36/46 · 2,590 |
| plain | V d3 p1 | 11/12 · 2.9 | 15/17 · 30.1 | 19/23 · 143 | 28/34 · 988 | 34/43 · 3,274 |
| plain | K d1 p1 | 15/18 · 2.5 | 16/21 · 32.5 | 18/24 · 158 | 17/23 · 896 | 17/24 · 2,780 |
| plain | K d2 p1 | 12/14 · 2.8 | 13/15 · 30.4 | 13/16 · 138 | 13/16 · 837 | 13/17 · 2,597 |
| plain | K d3 p1 | 11/12 · 3.0 | 11/13 · 40.2 | 11/14 · 186 | 11/14 · 1,117 | 11/15 · 3,587 |
| plain | V d1 p2 | 18/21 · 2.2 | 23/28 · 16.3 | 24/30 · 61.5 | 34/44 · 421 | 49/65 · 1,597 |
| plain | V d2 p2 | 15/17 · 2.8 | 21/24 · 22.6 | 22/26 · 84.3 | 32/39 · 586 | 45/58 · 2,311 |
| plain | V d3 p2 | 13/14 · 3.0 | 18/21 · 26.4 | 20/24 · 106 | 30/36 · 734 | 42/54 · 2,782 |
| plain | K d1 p2 | 18/21 · 2.4 | 20/24 · 17.8 | 21/27 · 71.7 | 20/28 · 371 | 22/30 · 1,073 |
| plain | K d2 p2 | 15/17 · 2.8 | 17/19 · 21.6 | 17/20 · 80.6 | 17/21 · 445 | 17/22 · 1,195 |
| plain | K d3 p2 | 13/14 · 3.2 | 15/17 · 25.5 | 15/18 · 96.7 | 15/19 · 526 | 15/19 · 1,368 |
| smoothed | V d2 p3 α10 | 9/10 · 1.9 | 10/11 · 10.8 | 10/12 · 43.2 | 10/12 · 199 | 10/13 · 556 |
| smoothed | V d3 p3 α10 | 8/9 · 2.1 | 8/10 · 12.5 | 8/10 · 46.9 | 9/11 · 234 | 9/11 · 606 |
| smoothed | K d2 p3 α10 | 9/10 · 1.8 | 10/11 · 11.0 | 9/12 · 45.6 | 10/12 · 207 | 10/13 · 634 |
| smoothed | K d3 p3 α10 | 8/9 · 2.2 | 8/10 · 12.7 | 8/10 · 49.0 | 9/10 · 226 | 8/10 · 598 |
| smoothed | V d2 p3 α30 | 13/15 · 2.7 | 13/16 · 15.6 | 14/16 · 57.6 | 14/17 · 278 | 14/18 · 766 |
| smoothed | V d3 p3 α30 | 9/10 · 2.4 | 10/11 · 13.8 | 10/12 · 56.1 | 10/13 · 278 | 10/13 · 710 |
| smoothed | K d2 p3 α30 | 13/15 · 2.8 | 13/16 · 15.9 | 14/17 · 70.0 | 14/17 · 338 | 14/18 · 906 |
| smoothed | K d3 p3 α30 | 9/10 · 2.5 | 10/11 · 14.1 | 10/12 · 58.1 | 10/12 · 279 | 10/13 · 744 |

**irregular, warm start after a 2% `q` step** — iterations (x / z) · solve ms

| method | config | 64 | 128 | 224 | 448 | 708 |
|---|---|---:|---:|---:|---:|---:|
| plain | V d1 p1 | 10/14 · 2.0 | 11/16 · 12.7 | 13/20 · 55.7 | 16/29 · 375 | 17/36 · 1,222 |
| plain | V d2 p1 | 7/10 · 2.0 | 9/14 · 18.1 | 10/18 · 81.9 | 13/27 · 570 | 15/33 · 1,824 |
| plain | V d3 p1 | 6/9 · 2.3 | 8/13 · 23.1 | 9/16 · 100 | 12/25 · 726 | 13/31 · 2,368 |
| plain | K d1 p1 | 10/14 · 2.1 | 10/15 · 23.3 | 11/18 · 118 | 10/17 · 660 | 11/19 · 2,259 |
| plain | K d2 p1 | 7/10 · 2.0 | 8/12 · 24.4 | 8/12 · 101 | 7/13 · 646 | 7/13 · 1,901 |
| plain | K d3 p1 | 6/9 · 2.3 | 7/10 · 28.9 | 7/11 · 134 | 7/11 · 799 | 6/11 · 2,366 |
| plain | V d1 p2 | 11/16 · 1.8 | 13/21 · 12.4 | 14/22 · 45.9 | 17/32 · 305 | 22/45 · 1,114 |
| plain | V d2 p2 | 9/13 · 2.2 | 11/18 · 16.9 | 11/19 · 62.2 | 15/29 · 436 | 18/41 · 1,620 |
| plain | V d3 p2 | 8/11 · 2.5 | 10/16 · 20.5 | 10/18 · 79.7 | 13/26 · 538 | 16/38 · 1,974 |
| plain | K d1 p2 | 11/16 · 1.9 | 12/18 · 13.4 | 13/20 · 53.4 | 12/21 · 287 | 12/23 · 810 |
| plain | K d2 p2 | 9/13 · 2.2 | 9/15 · 17.0 | 9/15 · 60.1 | 9/16 · 335 | 9/16 · 858 |
| plain | K d3 p2 | 8/11 · 2.5 | 8/13 · 19.7 | 8/13 · 69.4 | 8/14 · 377 | 8/15 · 1,056 |
| smoothed | V d2 p3 α10 | 5/8 · 1.5 | 6/8 · 8.0 | 6/9 · 32.8 | 5/9 · 150 | 5/10 · 434 |
| smoothed | V d3 p3 α10 | 5/6 · 1.5 | 5/7 · 9.0 | 5/7 · 33.5 | 5/8 · 174 | 5/8 · 454 |
| smoothed | K d2 p3 α10 | 5/8 · 1.6 | 6/8 · 8.2 | 5/9 · 34.2 | 5/9 · 161 | 5/10 · 467 |
| smoothed | K d3 p3 α10 | 5/6 · 1.5 | 5/7 · 9.1 | 5/7 · 34.5 | 4/8 · 185 | 4/8 · 484 |
| smoothed | V d2 p3 α30 | 8/12 · 2.2 | 8/12 · 11.8 | 8/13 · 46.5 | 8/13 · 217 | 8/14 · 602 |
| smoothed | V d3 p3 α30 | 6/8 · 1.9 | 6/9 · 11.3 | 6/9 · 42.6 | 6/10 · 222 | 6/10 · 567 |
| smoothed | K d2 p3 α30 | 8/12 · 2.3 | 8/12 · 11.8 | 8/13 · 52.3 | 8/13 · 240 | 8/15 · 685 |
| smoothed | K d3 p3 α30 | 6/8 · 1.9 | 6/9 · 11.8 | 6/9 · 43.4 | 6/10 · 236 | 6/10 · 586 |

**dome, cold start** — iterations (x column / z column) · solve ms

| method | config | 64 | 128 | 224 | 448 | 708 |
|---|---|---:|---:|---:|---:|---:|
| plain | V d1 p1 | 16/19 · 2.6 | 21/26 · 16.9 | 27/33 · 78.3 | 38/50 · 531 | 53/72 · 2,020 |
| plain | V d2 p1 | 13/14 · 2.7 | 17/20 · 20.6 | 22/26 · 93.8 | 31/41 · 687 | 46/59 · 2,582 |
| plain | V d3 p1 | 12/13 · 3.1 | 16/19 · 26.6 | 21/25 · 124 | 30/37 · 847 | 42/54 · 3,286 |
| plain | K d1 p1 | 16/19 · 2.7 | 17/21 · 28.9 | 18/23 · 135 | 18/24 · 920 | 18/25 · 3,366 |
| plain | K d2 p1 | 13/14 · 2.7 | 14/16 · 32.3 | 14/16 · 139 | 14/17 · 988 | 14/17 · 3,565 |
| plain | K d3 p1 | 12/13 · 3.1 | 13/15 · 40.4 | 13/16 · 187 | 14/17 · 1,318 | 14/18 · 5,003 |
| plain | V d1 p2 | 18/21 · 2.2 | 26/32 · 15.4 | 29/37 · 62.4 | 43/56 · 429 | 58/83 · 1,686 |
| plain | V d2 p2 | 16/17 · 2.7 | 21/24 · 17.7 | 24/30 · 78.0 | 36/45 · 530 | 50/66 · 2,132 |
| plain | V d3 p2 | 13/15 · 3.2 | 20/22 · 22.1 | 23/27 · 93.7 | 33/41 · 646 | 47/59 · 2,487 |
| plain | K d1 p2 | 18/21 · 2.4 | 22/27 · 17.2 | 23/30 · 67.7 | 24/32 · 365 | 25/35 · 1,233 |
| plain | K d2 p2 | 16/17 · 2.9 | 17/20 · 19.1 | 18/22 · 74.4 | 19/24 · 400 | 20/25 · 1,265 |
| plain | K d3 p2 | 13/15 · 3.2 | 16/19 · 24.2 | 17/20 · 88.7 | 18/22 · 484 | 18/24 · 1,534 |
| smoothed | V d2 p3 α10 | 9/11 · 1.9 | 9/11 · 8.9 | 11/13 · 38.4 | 11/14 · 174 | 11/14 · 463 |
| smoothed | V d3 p3 α10 | 8/9 · 2.0 | 8/10 · 10.2 | 10/12 · 46.7 | 10/12 · 193 | 10/13 · 559 |
| smoothed | K d2 p3 α10 | 9/11 · 2.0 | 9/11 · 9.2 | 11/13 · 42.3 | 11/13 · 182 | 11/14 · 517 |
| smoothed | K d3 p3 α10 | 8/9 · 2.1 | 8/10 · 10.3 | 9/11 · 45.4 | 10/12 · 209 | 10/12 · 573 |
| smoothed | V d2 p3 α30 | 13/15 · 2.5 | 14/16 · 12.6 | 14/17 · 50.3 | 15/19 · 236 | 18/22 · 723 |
| smoothed | V d3 p3 α30 | 10/11 · 2.5 | 10/12 · 12.1 | 14/16 · 61.5 | 18/21 · 336 | 25/30 · 1,290 |
| smoothed | K d2 p3 α30 | 13/15 · 2.7 | 14/16 · 13.2 | 14/17 · 62.8 | 14/18 · 281 | 15/19 · 842 |
| smoothed | K d3 p3 α30 | 10/11 · 2.5 | 10/12 · 12.3 | 12/14 · 65.6 | 12/15 · 298 | 14/20 · 1,116 |

**dome, warm start after a 2% `q` step** — iterations (x / z) · solve ms

| method | config | 64 | 128 | 224 | 448 | 708 |
|---|---|---:|---:|---:|---:|---:|
| plain | V d1 p1 | 10/14 · 1.9 | 13/19 · 12.4 | 15/24 · 55.9 | 19/35 · 371 | 25/50 · 1,397 |
| plain | V d2 p1 | 8/11 · 2.1 | 10/15 · 15.5 | 12/19 · 68.8 | 16/29 · 480 | 20/42 · 1,864 |
| plain | V d3 p1 | 7/10 · 2.4 | 9/14 · 20.0 | 10/18 · 90.2 | 15/27 · 611 | 18/38 · 2,359 |
| plain | K d1 p1 | 10/14 · 2.0 | 11/16 · 22.0 | 11/17 · 100 | 10/18 · 699 | 10/18 · 2,535 |
| plain | K d2 p1 | 8/11 · 2.1 | 8/12 · 21.4 | 8/13 · 97.7 | 8/13 · 666 | 8/14 · 2,524 |
| plain | K d3 p1 | 7/10 · 2.5 | 8/12 · 30.0 | 8/12 · 129 | 8/13 · 908 | 8/13 · 3,231 |
| plain | V d1 p2 | 12/16 · 1.7 | 16/24 · 11.5 | 17/27 · 45.5 | 21/40 · 307 | 27/56 · 1,166 |
| plain | V d2 p2 | 10/13 · 2.1 | 12/18 · 13.5 | 13/22 · 57.9 | 16/32 · 382 | 23/46 · 1,501 |
| plain | V d3 p2 | 9/11 · 2.3 | 11/17 · 17.2 | 11/19 · 65.7 | 15/29 · 460 | 20/42 · 1,785 |
| plain | K d1 p2 | 12/16 · 1.8 | 13/20 · 12.8 | 13/22 · 50.9 | 13/23 · 264 | 14/25 · 889 |
| plain | K d2 p2 | 10/13 · 2.2 | 10/15 · 14.2 | 10/16 · 53.2 | 10/18 · 305 | 10/19 · 942 |
| plain | K d3 p2 | 9/11 · 2.4 | 9/14 · 17.6 | 9/15 · 64.7 | 9/16 · 346 | 9/17 · 1,060 |
| smoothed | V d2 p3 α10 | 6/8 · 1.4 | 6/8 · 6.6 | 6/10 · 29.6 | 6/10 · 131 | 6/10 · 342 |
| smoothed | V d3 p3 α10 | 5/7 · 1.6 | 5/7 · 7.3 | 5/9 · 35.2 | 5/9 · 152 | 5/9 · 400 |
| smoothed | K d2 p3 α10 | 6/8 · 1.5 | 6/8 · 6.6 | 6/9 · 29.4 | 6/10 · 142 | 6/10 · 368 |
| smoothed | K d3 p3 α10 | 5/7 · 1.6 | 5/7 · 7.4 | 5/8 · 32.6 | 5/9 · 159 | 5/9 · 422 |
| smoothed | V d2 p3 α30 | 9/12 · 2.1 | 9/12 · 9.5 | 8/13 · 39.0 | 9/14 · 176 | 11/18 · 611 |
| smoothed | V d3 p3 α30 | 7/9 · 2.1 | 6/9 · 9.2 | 8/12 · 46.2 | 12/17 · 276 | 17/24 · 1,057 |
| smoothed | K d2 p3 α30 | 9/12 · 2.2 | 9/12 · 10.2 | 8/13 · 47.8 | 8/14 · 221 | 8/15 · 672 |
| smoothed | K d3 p3 α30 | 7/9 · 2.1 | 6/9 · 9.5 | 7/10 · 47.1 | 6/11 · 222 | 7/15 · 840 |

### 1M edges: the candidates against the direct solver

Setup, per-`q` hierarchy update and solve times in ms. The direct column is
the refactor + one 3-rhs solve on the same system; "/ direct" is the solve
time alone over that.

| fixture | direct | config | setup | q-update | cold it x/y/z | cold ms | cold / direct | warm it x/y/z | warm ms | warm / direct |
|---|---:|---|---:|---:|---|---:|---:|---|---:|---:|
| grid 708 | 756 | plain V d3 p2 | 163 | 7.6 | 60/62/87 | 4,530 | 6.00× | 34/37/66 | 3,626 | 4.80× |
| grid 708 | 756 | plain K d1 p2 | 164 | 8.0 | 35/34/50 | 2,668 | 3.53× | 22/22/37 | 1,753 | 2.32× |
| grid 708 | 756 | plain K d3 p2 | 161 | 7.8 | 23/23/33 | 2,806 | 3.71× | 14/14/25 | 2,206 | 2.92× |
| grid 708 | 756 | smoothed V d1 p3 α10 | 640 | 395 | 20/20/29 | 1,049 | 1.39× | 14/14/23 | 814 | 1.08× |
| grid 708 | 756 | smoothed V d2 p3 α10 | 548 | 367 | 15/15/22 | 920 | 1.22× | 10/10/17 | 743 | 0.98× |
| grid 708 | 756 | smoothed V d3 p3 α10 | 543 | 367 | 13/13/18 | 982 | 1.30× | 8/8/14 | 774 | 1.03× |
| grid 708 | 756 | smoothed K d3 p3 α10 | 543 | 368 | 12/12/17 | 1,009 | 1.34× | 7/7/13 | 774 | 1.03× |
| irregular 546 | 358 | plain V d3 p2 | 186 | 10.0 | 42/42/54 | 2,782 | 7.78× | 16/15/38 | 1,974 | 5.52× |
| irregular 546 | 358 | plain K d1 p2 | 187 | 9.7 | 22/21/30 | 1,073 | 3.00× | 12/12/23 | 810 | 2.27× |
| irregular 546 | 358 | plain K d3 p2 | 185 | 9.9 | 15/15/19 | 1,368 | 3.82× | 8/7/15 | 1,056 | 2.95× |
| irregular 546 | 358 | smoothed V d1 p3 α10 | 538 | 367 | 14/14/20 | 596 | 1.67× | 8/8/15 | 470 | 1.32× |
| irregular 546 | 358 | smoothed V d2 p3 α10 | 534 | 322 | 10/10/13 | 556 | 1.55× | 5/5/10 | 434 | 1.21× |
| irregular 546 | 358 | smoothed V d3 p3 α10 | 521 | 335 | 9/9/11 | 606 | 1.69× | 5/4/8 | 454 | 1.27× |
| irregular 546 | 358 | smoothed K d3 p3 α10 | 532 | 358 | 8/8/10 | 598 | 1.67× | 4/4/8 | 484 | 1.35× |
| dome 342 | 458 | plain V d3 p2 | 134 | 7.1 | 47/45/59 | 2,487 | 5.44× | 20/19/42 | 1,785 | 3.90× |
| dome 342 | 458 | plain K d1 p2 | 138 | 7.3 | 25/25/35 | 1,233 | 2.70× | 14/13/25 | 889 | 1.94× |
| dome 342 | 458 | plain K d3 p2 | 131 | 7.1 | 18/18/24 | 1,534 | 3.35× | 9/9/17 | 1,060 | 2.32× |
| dome 342 | 458 | smoothed V d1 p3 α10 | 536 | 382 | 15/15/20 | 506 | 1.11× | 8/8/15 | 392 | 0.86× |
| dome 342 | 458 | smoothed V d2 p3 α10 | 471 | 323 | 11/11/14 | 463 | 1.01× | 6/6/10 | 342 | 0.75× |
| dome 342 | 458 | smoothed V d3 p3 α10 | 462 | 321 | 10/10/13 | 559 | 1.22× | 5/5/9 | 400 | 0.87× |
| dome 342 | 458 | smoothed K d3 p3 α10 | 461 | 323 | 10/10/12 | 573 | 1.25× | 5/5/9 | 422 | 0.92× |

### α sweep, smoothed aggregation, grid

Cold iterations x/z · ms (warm iterations in parentheses).

| grid | config | α = 5 | α = 10 | α = 30 | α = 100 |
|---:|---|---:|---:|---:|---:|
| 224 | smoothed V d2 p3 | 13/18 · 70.1 (warm 9/13) | 14/19 · 74.0 (warm 9/14) | 19/25 · 97.0 (warm 13/19) | 30/40 · 152 (warm 21/31) |
| 224 | smoothed V d3 p3 | 12/17 · 86.2 (warm 8/12) | 12/16 · 80.6 (warm 8/12) | 14/18 · 89.3 (warm 9/14) | 20/26 · 131 (warm 14/20) |
| 708 | smoothed V d2 p3 | 14/21 · 948 (warm 9/16) | 15/22 · 920 (warm 10/17) | 21/29 · 1,302 (warm 13/23) | 33/45 · 1,932 (warm 21/36) |
| 708 | smoothed V d3 p3 | 13/19 · 1,097 (warm 8/14) | 13/18 · 982 (warm 8/14) | 15/21 · 1,192 (warm 10/16) | 22/30 · 1,716 (warm 14/24) |

### 100× force-density ratio

Same fixtures at 224 with `q` spanning 100× instead of 10× across the
domain (plus the grid at 708). Iteration counts are within ±1 of the 10×
runs for every configuration (the 708 timings here were taken under load
and are 15–25% above the table above).

| fixture | config | cold it x/z · ms | warm it x/z · ms |
|---|---|---:|---:|
| grid 224 | plain V d2 p2 | 42/55 · 208 | 27/42 · 159 |
| grid 224 | plain V d3 p2 | 33/45 · 230 | 21/34 · 174 |
| grid 224 | plain K d2 p2 | 31/42 · 232 | 19/31 · 172 |
| grid 224 | plain K d3 p2 | 20/27 · 203 | 13/20 · 148 |
| grid 224 | smoothed V d2 p3 α10 | 13/18 · 74.1 | 9/14 · 57.8 |
| grid 224 | smoothed V d3 p3 α10 | 12/15 · 78.4 | 7/12 · 63.8 |
| grid 224 | smoothed K d2 p3 α10 | 13/17 · 80.1 | 9/13 · 56.9 |
| grid 224 | smoothed K d3 p3 α10 | 11/15 · 84.1 | 7/11 · 61.4 |
| irregular 172 | plain V d2 p2 | 22/27 · 89.7 | 11/20 · 66.4 |
| irregular 172 | plain V d3 p2 | 20/24 · 108 | 10/18 · 81.5 |
| irregular 172 | plain K d2 p2 | 17/21 · 86.2 | 10/16 · 65.4 |
| irregular 172 | plain K d3 p2 | 15/18 · 98.4 | 8/14 · 76.0 |
| irregular 172 | smoothed V d2 p3 α10 | 10/12 · 44.7 | 6/9 · 33.0 |
| irregular 172 | smoothed V d3 p3 α10 | 9/10 · 48.4 | 5/8 · 40.1 |
| irregular 172 | smoothed K d2 p3 α10 | 10/12 · 47.1 | 6/9 · 34.8 |
| irregular 172 | smoothed K d3 p3 α10 | 8/10 · 53.9 | 5/7 · 38.4 |
| dome 103 | plain V d2 p2 | 23/30 · 76.8 | 12/22 · 57.2 |
| dome 103 | plain V d3 p2 | 22/27 · 93.8 | 11/19 · 66.9 |
| dome 103 | plain K d2 p2 | 17/22 · 73.6 | 10/16 · 53.5 |
| dome 103 | plain K d3 p2 | 16/20 · 87.4 | 9/15 · 64.9 |
| dome 103 | smoothed V d2 p3 α10 | 11/13 · 38.4 | 6/10 · 29.9 |
| dome 103 | smoothed V d3 p3 α10 | 10/12 · 46.7 | 5/8 · 31.7 |
| dome 103 | smoothed K d2 p3 α10 | 10/13 · 41.9 | 6/9 · 29.0 |
| dome 103 | smoothed K d3 p3 α10 | 9/11 · 45.4 | 5/8 · 32.9 |
| grid 708 | plain K d3 p2 | 22/32 · 3,203 | 14/24 · 2,521 |
| grid 708 | smoothed V d2 p3 α10 | 15/21 · 1,141 | 9/16 · 888 |
| grid 708 | smoothed V d3 p3 α10 | 13/18 · 1,261 | 8/14 · 1,022 |

### Over-correction of the coarse-grid correction, plain V d3 p2, grid

| grid | over-correction | cold it x/z · ms | warm it x/z · ms |
|---:|---:|---:|---:|
| 224 | 1.0 | 37/47 · 236 | 23/36 · 181 |
| 224 | 1.3 | 29/38 · 190 | 18/28 · 140 |
| 224 | 1.5 | 26/34 · 170 | 17/26 · 131 |
| 224 | 1.8 | 25/33 · 167 | 16/24 · 129 |
| 708 | 1.0 | 60/87 · 5,052 | 34/66 · 3,725 |
| 708 | 1.5 | 35/51 · 2,775 | 20/38 · 2,054 |

### Reading the tables

Iteration counts below are for the x and z columns (y behaves like x); the
z column (the load) always needs the most iterations and is what decides the
wall time.

**Plain aggregation, V-cycle** is not size-independent on any fixture. On
the grid the best plain V configuration (d3 p2) goes from 17/21 iterations
at 64 to 60/87 at 708 (3.5–4×); irregular 13/14 → 42/54, dome 13/15 →
47/59. The growth is with the number of levels (2 → 5), the signature of an
aggregation hierarchy whose coarse operators under-represent the fine ones
(piecewise-constant `P` on a Laplacian gives coarse operators about a
factor ~2 too stiff per level); more smoothing (d1 → d3) or one matching
pass (~2× coarsening, 8–9 levels) moves the counts but not the trend.
Over-correction of the coarse-grid correction (`--over`, the textbook
remedy for plain aggregation; the preconditioner stays symmetric) helps
more than any smoother change — a factor 1.5 takes the grid from 37/47 to
26/34 iterations at 224 and from 60/87 to 35/51 at 708 — but the growth
with size remains (+35% / +50% from 224 to 708 against +62% / +85%
without) and the 708 solve is still 2.8 s, 3.7× the direct reference.

**Plain aggregation, K-cycle** does what `AMG_PLAN.md` says it does for the
iteration counts: with two matching passes, d3, they are 17/21 → 23/33 on
the grid (+35% / +57%), 13/14 → 15/19 on irregular (+15% / +36%), 13/15 →
18/24 on the dome (+38% / +60%). That is inside the ±50% acceptance band
on the x column for irregular and dome and just outside ±30% on the grid,
and outside on the z column everywhere but irregular. The cost is the
problem: a K-cycle on a 5-level hierarchy is 8–15 fine-operator
applications per outer iteration (d2–d3, two passes; the `ktol` early exit
keeps it from doubling per level), about twice the V-cycle's 6–7, and with
one matching pass (8–9 levels) the recursion runs away (23–53 applications
per iteration). At 1M edges the plain K-cycle with two passes needs
2.7–2.8 s (grid), 1.1–1.4 s (irregular) and 1.2–1.5 s (dome) cold
(degree 1 is the cheapest in wall time, degree 3 in iterations) against
756 / 358 / 458 ms for the direct refactor + solve: 2.7–3.8× slower cold,
1.9–3.0× slower warm.

**Smoothed aggregation** removes the level dependence. `P = (I − ω D⁻¹A)
P₀` with three matching passes (aggregates of ~8, coarsening ~8× per level,
4 levels at 500k nodes), Chebyshev degree 2, α = 10, V-cycle inside plain
PCG: grid 12/15 → 15/22 (64 → 708; flat from 448 to 708 at 15–16/22),
irregular 9/10 → 10/13, dome 9/11 → 11/14; warm 8/11 → 10/17, 5/8 → 5/10,
6/8 → 6/10. x/y columns are within ±30% on the grid and ±25% elsewhere;
the z column is +47% on the grid (three levels added between 64 and 448,
then flat) and +27–30% on irregular and dome. The K-cycle buys at most one
iteration on top of smoothed aggregation at 1.3–1.5× the work per
iteration, so PCG + V-cycle is the right pairing and flexible CG is not
needed. Degree 3 saves 1–4 iterations but costs 40% more per iteration and
is slower in wall time on all three fixtures at 1M; degree 1 (1M-edge
table) needs 30–45% more iterations than degree 2 and is 7–10% slower in
wall time. α = 5 and α = 10 give the same counts, α = 30 costs +40%
iterations and α = 100 doubles them
(`[λ_max/α, λ_max]` is the range the smoother damps; the aggregates of 8
leave nothing below `λ_max/10` for the smoother to do). α also matters for
robustness: with α = 30 the degree-3 V-cycle on the dome loses its
size-independence (10/11 → 25/30 iterations) where α = 10 keeps it (8/9 →
10/13). The 100× `q` ratio changes nothing (±1 iteration on every
configuration); a 1e4 ratio was not run.

Wall time at 1M edges, smoothed V d2 p3 α10: 920 ms cold / 743 ms warm on
the grid against 756 ms direct (1.22× / 0.98×), 556 / 434 against 358 on
irregular (1.55× / 1.21×), 463 / 342 against 458 on the dome (1.01× /
0.75×). To that the smoothed hierarchy adds 320–375 ms per `q` change
(rebuilding `P` and the three Galerkin products), where plain aggregation
needs 7–10 ms. Per optimizer evaluation (one `q` update, a
forward and an adjoint solve, both counted cold) that is ~2.2 s on the
grid against the direct path's 0.83 s (661 ms factor + two 86 ms solves,
table above), 1.4 s against ~0.39 s on irregular and 1.25 s against
~0.50 s on the dome: the single-threaded prototype is 2.5–3.6× slower than
the direct solver per evaluation at 1M edges, not at parity. Where the time
goes at grid 708: the operator application on three right-hand sides costs
4.2 ms (memory bound: ~24 MB of CSR adjacency plus ~24 MB of vector
traffic per application, ~11 GB/s); a V-cycle iteration is 6.8
operator-equivalents ≈ 29 ms of operator time and 42 ms measured, the
remainder being the
Chebyshev vector updates, the smoothed prolongation/restriction (several
entries per fine row instead of one) and the PCG dot products. Setup
(aggregation, power iterations, `P`, Galerkin products, coarsest factor)
is 0.48–0.64 s at 1M and is paid once per topology.

### Recommendation

* **Plain (unsmoothed) aggregation: no-go.** V-cycle iteration counts grow
  3.5–4× from 8k to 1M edges on every fixture; the K-cycle keeps them
  roughly flat (+35% x / +57% z on the grid, +15–38% x on irregular and
  dome) but at 8–15 operator applications per iteration it is 2.7–3.8×
  slower than the direct refactor + solve at 1M edges single-threaded and
  fails the wall-time criterion of Phase 0 by that margin. Nothing in the
  parameter space measured (degree 1–3, 1–2 passes, α 10–100, ktol,
  over-correction) changes this.
* **Smoothed aggregation: fixes the growth, go for WS-C with it.**
  Recommended parameters: `P = (I − ω D⁻¹A) P₀`, `ω = 4/(3 λ_max)`; three
  pairwise matching passes per level (coarsening ~8×, 4 levels at 500k free
  nodes, coarse operators 0.27–0.45× the fine edge count, 18–28 MB at 1M);
  Chebyshev degree 2 with Jacobi scaling and α = 10; V-cycle as a fixed
  preconditioner inside standard PCG on the block of three right-hand
  sides with per-column convergence; coarsest level ≤ 2,000 nodes through
  faer's sparse Cholesky. Drop the K-cycle and flexible CG.
* Iteration flatness against the acceptance bands: grid x/y +25% (inside
  ±30%), grid z +47% (outside, but flat from 448 to 708); irregular and
  dome +11–30% on every column (inside ±50%). Warm starts after a 2% `q`
  step save 33–50% of the x/y iterations but only 23–29% of the z
  iterations that set the wall time, hence 19–26% of the solve time.
* The `hierarchy.rs` of the program document needs the general sparse
  level type of its risk register (CSR values, `A·x` as a sparse matvec on
  three right-hand sides, restriction/prolongation through a sparse `P`,
  `Pᵀ A P` per level per `q` update), not only the graph-Laplacian level.
  Keep the plain level type as the cheap `q`-update path only if WS-C finds
  a way to reuse it (below).
* Wall time at 1M edges single-threaded is 1.0–1.55× the direct refactor
  + solve for the solve alone and 2.5–3.6× per evaluation once the
  0.32–0.38 s smoothed-hierarchy update and the second (adjoint) solve are
  counted. The single-thread parity that `AMG_PLAN.md` expected at 1M is
  therefore not reached; the case for the
  iterative path at 1M rests on the parallel scaling of memory-bound
  kernels that faer's factorization does not have (WS-B), and on the
  linear scaling beyond 1M where the direct factor's `n^1.5` flops take
  over.

### Where the measurements disagree with the plans

* `AMG_PLAN.md`: "aggregation AMG on graph Laplacians converges in 10–20
  iterations independently of size". True only for the K-cycle or for
  smoothed aggregation; the plain V-cycle needs 34 → 105 iterations (d1
  p1, z column, grid 64 → 708).
* `AMG_PLAN.md` budgets the K-cycle at "≈ 10 operator-equivalents" and an
  FCG iteration at 17 ms with 4 threads at 1M. Measured: 11–15.4
  operator-equivalents with two matching passes, 23–53 with one, and 85 ms
  per iteration single-threaded (grid, d3 p2). Its Phase-0 exit criterion
  "wall time at 708 beats the direct 757 ms single-threaded" is missed by
  plain aggregation by 3–6× and met by smoothed aggregation only for warm
  solves on the grid and dome, and by nothing once the hierarchy update is
  counted.
* `AMG_PLAN.md`: coarse weights "recomputed per evaluation in one `O(nnz)`
  pass per level — no sparse matrix products". Holds for plain aggregation
  (7–10 ms at 1M) and not for the variant that converges: the smoothed
  hierarchy costs 320–380 ms per `q` change at 1M in this prototype, 40–70%
  of a cold solve.
* `AMG_PLAN.md`: warm start "cuts iterations by 2–3×" late in an
  optimisation. For a 2% `q` step the x/y columns drop 1.5–2× (15 → 10 on
  the grid, 10 → 5 on irregular, 11 → 6 on the dome) but the z column,
  which sets the wall time, only 1.3–1.4× (22 → 17, 13 → 10, 14 → 10); a
  10% step costs one or two iterations more than 2%. Smaller late-stage
  steps were not measured.
* `ITERATIVE_SOLVER_PROGRAM.md` WS-A acceptance band ±30% (grid) and ±50%
  (irregular): the smoothed V-cycle meets the x/y column on every fixture
  and the z column on irregular and dome; grid z is +47%. The plain
  K-cycle meets the irregular band and misses the grid band.
* What the plans got right: `A·x` on three right-hand sides costs
  3.0–5.1 ms at 1M (plan: 0.4 ms per 100k edges); two matching passes give
  5 levels at 500k nodes (plan: 7–8 at 5M); the hierarchy holds 0.25–0.45×
  the fine edges (plan: ~1/3); 100× `q` ratios are harmless; strength-based
  matching, the Chebyshev/Jacobi smoother and the faer coarsest solve work
  as described; the V-cycle preconditioner is symmetric to 1e-10 and the
  coarse operators equal `PᵀAP` (`--self-test`).

### Open questions

* The smoothed-hierarchy `q` update. Options to measure in WS-C: keep `P`
  (pattern and values) from the previous `q` and redo only the numeric
  triple product on the fixed pattern (any fixed `P` gives a valid SPD
  Galerkin operator, and a 2% `q` step barely moves `D⁻¹A`); recompute `P`
  only when `q` has drifted by more than a threshold; or a plain-aggregation
  update of the coarse *graph* combined with smoothing applied on the fly
  in prolongation/restriction (`P₀ᵀ(I − ωAD⁻¹) r`) with the coarse operator
  taken as the quotient graph, which changes the coarse operator and would
  have to be re-measured for convergence.
* The z column's +47% on the grid 64 → 708 comes with the levels added
  between 64 and 448 and is flat afterwards; whether a coarser stop
  (`coarsest` 5,000–10,000 with a sparse Cholesky coarsest solve) or a
  W-cycle on the two coarsest levels removes it is cheap to test in the
  same example.
* Parallel scaling. Every kernel here is memory bound (the operator streams
  ~50 MB per application at 1M); 4 threads on this VM buy at most the
  bandwidth ratio, to be measured by WS-B. The direct factorization had
  nothing to gain from faer's parallel path (table above), so the crossover
  per evaluation is expected between 1M and 3M edges on 4 threads, and
  only beyond ~3M edges single-threaded (direct `n^1.5`, iterative `n`).
* The direct reference is 2× cheaper on the irregular and dome fixtures than
  on the grid at equal edge count (fewer free nodes, less fill), so the
  crossover is topology dependent; a fixture with very few supports and
  long cables (small `anchor` mass, near-singular `A`) was not measured and
  is where aggregation AMG is known to slow down.
* Not measured: `q` ratios of 1e4 (`--qratio`), tolerances looser than 1e-8
  (the adaptive tolerance of `AMG_PLAN.md` would take the typical solve to
  1e-6, roughly 60% of the iterations), and the accuracy of the gradient
  computed from an inexact adjoint.

## WS-C: smoothed-aggregation AMG solver

Workstream WS-C of `ITERATIVE_SOLVER_PROGRAM.md`: the production
`AmgSolver<CpuBackend>` in `src/amg/` behind `LinearSolverKind::IterativeCpu`,
measured against `DirectSolver` on the same systems and against the WS-A
prototype numbers above. The configuration is the Phase-0 recommendation and
the crate defaults: smoothed aggregation, 3 matching passes per level,
`P = (I − ω D⁻¹A) P₀` with `ω = 4/(3 λ_max)`, Chebyshev degree 2 on
`[λ_max/10, λ_max]`, V-cycle (1 pre, 1 post), coarsest level ≤ 2,000 nodes
solved by faer's LLᵀ, PCG on the three columns to a fixed `1e-8` relative
residual.

```sh
# Ignored benchmark test; one markdown row per fixture.
for th in 1 4; do
  env RAYON_NUM_THREADS=$th THESEUS_AMG_SIZES=224,448,708 \
      THESEUS_AMG_FIXTURES=grid,irregular,dome THESEUS_AMG_REPS=5 \
      RUSTUP_TOOLCHAIN=1.89.0 \
      cargo test -p theseus --release --test amg_bench amg_vs_direct -- --ignored --nocapture
done
# Per-level breakdown of the numeric update:
env RAYON_NUM_THREADS=1 THESEUS_AMG_SIZES=708 THESEUS_AMG_FIXTURES=grid,irregular,dome \
    RUSTUP_TOOLCHAIN=1.89.0 \
    cargo test -p theseus --release --test amg_bench update_breakdown -- --ignored --nocapture
```

Method. Fixtures are the WS-B generators in `tests/support/fixtures`
(`grid`, `irregular`, `dome`, matched to the grid by edge count; their node
counts differ slightly from the prototype's own generators — the dome here
has 250k free nodes at 1M edges, the prototype's had 337k), with the smooth
10× force-density field `smooth_q_star` and unit downward loads. The direct
column is `DirectSolver::update` + one 3-rhs `solve` (numeric
refactorization + solve, the per-`q` cost the optimizer pays today). "Setup"
is the first `update(q)` on a fresh solver (aggregation, `P`, symbolic and
numeric Galerkin products, `λ_max`, coarsest factorization). "update(q)" is
`update` on a ±2% uniformly perturbed `q` (the pattern of `P` and of every
coarse operator is frozen; only the numeric products, the device weight
uploads, the coarsest refactorization and — when some `q_e` drifted more
than 50% since the last estimate, which a 2% step never does — `λ_max` are
redone).
"Cold" solves from `x = 0`; "warm" solves from the cold solution against
the right-hand side of the perturbed `q`. Every number is the median of 5
runs; the machine (4 vCPU Xeon, 15 GiB) was shared with one other build
during the sweep (load average 1.2–2.0 including this process; the direct
reference at grid 708 came out at 776–791 ms against the 756 ms of the
Phase-0 table taken at load ≈ 1). Iteration counts are deterministic and
identical across thread counts; the iterative solutions agree with the
direct ones to `1e-11`–`2e-9` relative (max-norm) at the `1e-8` residual
tolerance.

### update(q) — the headline number

Milliseconds for the numeric hierarchy update on a 2% `q` change, 1M edges
(grid 708 and matched fixtures), against the targets of §7.3 (≤ 100 ms on
one thread, ≤ 40 ms on four) and the prototype, which rebuilt `P` and the
coarse operators from scratch.

| fixture | edges | levels | update(q) 1 thread | update(q) 4 threads | prototype q-update (1 thread) | target |
|---|---:|---|---:|---:|---:|---|
| grid 708 | 1,001,112 | 501,260 → 62,577 → 7,772 → 929 | **34.5** | **20.1** | 367 | met |
| irregular 708 | 997,990 | 294,849 → 30,849 → 3,291 → 377 | **53.7** | **23.7** | 322 | met |
| dome 708 | 1,001,112 | 249,925 → 30,951 → 3,721 → 442 | **34.8** | **22.7** | 323 | met |

Where the single-threaded time goes (grid 708, 33 ms in this run): the
level-0 → 1 product 15.1 ms (the graph operator over `Level0Map` times
`P`, 2.4M entries in `AP`, then `Pᵀ(AP)` into 758k entries), level 1 → 2
4.1 ms, level 2 → 3 1.6 ms; the remaining ~12 ms are the level-0
weight/anchor pass over `Level0Map` and its device copy, the CSR value
uploads of the coarse operators, and the coarsest LLᵀ refactorization (929
nodes, ~1.5 ms). The irregular fixture spends 31 ms on the
first product against the grid's 15 ms for 0.6× the rows: its rows are
longer (14.7 vs 12.1 entries of `A_c` per coarse row, 5–9 neighbours per
fine node) and its neighbour accesses are not contiguous, so the product is
bound by cache misses rather than by arithmetic. The dome is the opposite:
its first product is cheaper than the grid's, its deeper levels dearer
(35–54 entries per coarse row). Four threads buy 1.5–2.3× on the update;
the products are row-parallel with a per-row deterministic accumulation
order, so the gain is bounded by memory bandwidth like everything else here.

The product is done in two stages on frozen patterns — `AP` row-parallel
over fine rows, then `Pᵀ(AP)` row-parallel over coarse rows — rather than
as a single triple loop per coarse row. The single loop visited every
`(u, v, j)` triple and was quadratic in the neighbourhood size: 529 ms on
the dome at 1M edges, 126 ms on the irregular fixture. Two stages cost one
extra matrix (`AP`, 2.4M entries / 27 MiB at level 0 of grid 708) and are
15× cheaper on the dome.

### Setup, solves and the direct reference

Single thread (`RAYON_NUM_THREADS=1`), coarsest ≤ 2,000 nodes.

| fixture | edges | free nodes | direct update + solve | AMG setup | AMG update(q) | cold it x/y/z | cold ms | cold / direct | warm it x/y/z | warm ms | warm / direct | host MiB |
|---|---:|---:|---:|---:|---:|---|---:|---:|---|---:|---:|---:|
| grid 224 | 99,904 | 50,172 | 43.3 | 31.2 | 3.7 | 15/15/20 | 99.0 | 2.29× | 10/10/15 | 74.3 | 1.72× | 36 |
| irregular 224 | 99,695 | 28,900 | 22.5 | 35.6 | 5.1 | 10/10/12 | 49.0 | 2.18× | 5/5/9 | 37.1 | 1.65× | 23 |
| dome 224 | 99,904 | 24,865 | 22.4 | 22.6 | 3.1 | 9/9/11 | 32.1 | 1.43× | 5/5/9 | 26.7 | 1.19× | 23 |
| grid 448 | 400,512 | 200,700 | 235.0 | 148.5 | 13.9 | 16/16/22 | 458.1 | 1.95× | 10/10/17 | 351.1 | 1.49× | 143 |
| irregular 448 | 400,357 | 117,649 | 113.3 | 152.6 | 21.2 | 10/10/13 | 232.4 | 2.05× | 5/5/10 | 178.6 | 1.58× | 96 |
| dome 448 | 400,512 | 99,905 | 122.6 | 96.1 | 13.3 | 9/9/12 | 152.1 | 1.24× | 5/5/9 | 113.7 | 0.93× | 95 |
| grid 708 | 1,001,112 | 501,260 | 791.4 | 407.3 | 34.5 | 17/16/23 | 1,318 | 1.67× | 10/10/18 | 1,034 | 1.31× | 358 |
| irregular 708 | 997,990 | 294,849 | 404.5 | 419.3 | 53.7 | 10/10/13 | 623.6 | 1.54× | 5/5/10 | 477.7 | 1.18× | 239 |
| dome 708 | 1,001,112 | 249,925 | 385.4 | 275.0 | 34.8 | 9/9/12 | 423.8 | 1.10× | 5/5/9 | 296.6 | 0.77× | 234 |

Four threads (`RAYON_NUM_THREADS=4`), same systems. The direct solver does
not benefit from threads (see the Phase-0 table); iteration counts are
identical to the single-threaded run, bitwise.

| fixture | direct update + solve | AMG setup | AMG update(q) | cold ms | cold / direct | warm ms | warm / direct |
|---|---:|---:|---:|---:|---:|---:|---:|
| grid 224 | 43.2 | 29.1 | 3.0 | 52.1 | 1.21× | 34.8 | 0.81× |
| irregular 224 | 22.2 | 32.0 | 3.3 | 25.3 | 1.14× | 18.5 | 0.83× |
| dome 224 | 22.5 | 20.4 | 2.7 | 18.6 | 0.83× | 15.1 | 0.67× |
| grid 448 | 232.6 | 110.9 | 8.1 | 174.6 | 0.75× | 133.9 | 0.58× |
| irregular 448 | 113.6 | 125.2 | 13.5 | 84.8 | 0.75× | 66.4 | 0.58× |
| dome 448 | 119.3 | 80.1 | 10.2 | 60.0 | 0.50× | 45.9 | 0.38× |
| grid 708 | 775.7 | 316.5 | 20.1 | 416.3 | 0.54× | 336.6 | 0.43× |
| irregular 708 | 397.5 | 373.6 | 23.7 | 202.9 | 0.51× | 161.1 | 0.41× |
| dome 708 | 392.8 | 236.3 | 22.7 | 139.6 | 0.36× | 106.3 | 0.27× |

Coarsest level ≤ 8,000 nodes instead of 2,000 (one level fewer at 1M
edges), single thread / four threads: update(q) 41.1 / 26.1 ms (grid),
55.3 / 26.3 (irregular), 37.5 / 24.6 (dome); cold solves 1,207 / 403,
609 / 197, 401 / 141 ms. The grid drops from 17/16/23 to 14/15/21
iterations, the other fixtures are unchanged; the larger LLᵀ (7.7k nodes,
31 entries per row) costs more in the update than the level it removes and
the solve time moves by ≤ 8%. The default stays at 2,000.

### Reading the tables

* The `update(q)` targets are met on every fixture with margin: 35–54 ms
  on one thread against 100, 20–24 ms on four against 40. Relative to the
  prototype's rebuild-everything update (322–367 ms) the frozen-pattern
  numeric product is 6–10× cheaper; relative to the direct update + solve
  it is 7–23× cheaper on one thread.
* Iteration counts match the prototype within one or two (grid 708:
  17/16/23 cold, 10/10/18 warm, against 15/15/22 and 10/10/17; the WS-B
  fixtures are not the prototype's generators, so the irregular and dome
  rows are not the same systems). The grid's z column still costs ~40%
  more iterations than x/y, as in Phase 0; the x/y counts at 64/128/224
  are within the +30% band checked by `tests/amg_solver.rs`.
* Per iteration the production solve is 25–35% dearer than the prototype's
  on one thread (grid 708: 57 ms per z-iteration against 42; grid 224: 5.0
  against 3.9). Part of this is load (the direct reference is 5% slower in
  the same runs), the rest is the generic kernel path: the smoothed
  restriction and prolongation are `apply_csr` products on `Pᵀ` and `P`
  (1.3M entries each at grid 708, one extra pass over the fine vectors
  compared with the prototype's fused aggregate-map form), and the three
  columns are interleaved so every kernel streams 3× the operator's
  vectors. The kernels are memory bound; four threads buy 3.0–3.2× on the
  cold solve at 1M edges.
* Against the direct solver per `q` evaluation (update + one 3-rhs solve),
  the AMG is 1.1–1.7× slower on one thread at 1M edges cold and 0.8–1.3×
  warm, and 0.36–0.54× (cold) / 0.27–0.43× (warm) on four threads; the
  crossover on one thread is the dome at 448 warm and lies beyond 1M
  edges for the grid (Phase 0 predicted 1–3M). The optimizer pays the
  adjoint solve on top of both, which favours the direct solver's second
  solve (a fraction of its refactorization) against a second PCG of the
  same cost — the accounting of §9 stands.
* Host memory (all level operators, `P`, `Pᵀ`, `AP`, the coarsest factor
  and the PCG/cycle buffers on the CPU backend): 358 MiB at grid 708, of
  which the level-0 `AP` product buffer is 27 MiB and `P` plus `Pᵀ` at
  level 0 are 30 MiB. The `f32` preconditioner option halves the device-side operators
  and vectors; it is not measured here.

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
