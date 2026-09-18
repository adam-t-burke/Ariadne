# Plan: matrix-free FDM solves with an aggregation AMG preconditioner

Goal: make Theseus form-finding feasible on networks of 1M–10M edges by
replacing the sparse Cholesky factorization of `A(q) = Cnᵀ diag(q) Cn` with a
preconditioned conjugate-gradient solve whose cost and memory are linear in
the number of edges. The direct solver stays as the default below a size
threshold and as the fallback whenever the system is not positive definite.

Numbers below are for the grid fixture in `tests/support/grid.rs` unless
stated; see `BENCHMARKS.md` for the measurements they extrapolate from.

## Why this and not something else

* At 1M edges the direct solve is 0.88 s per evaluation and 1.1 GB; 70% is
  the numeric Cholesky, which nothing around it can reduce (ordering,
  relaxation, faer's parallel path and `f32` were all measured). At 10M edges
  the factor would hold ~275M nonzeros (~2.2 GB) and cost 20–40 s per
  evaluation: fill grows like `n log n`, flops like `n^1.5`.
* The matrix-free operator costs 0.4 ms per 100k edges for three right-hand
  sides (one pass over the edges). The problem is conditioning: `A` is a
  weighted graph Laplacian, its condition number grows with the node count,
  and Jacobi-PCG needs ~1,200 iterations at 100k edges regardless of warm
  start.
* `A(q)` with `q > 0` and at least one support per component is a symmetric
  M-matrix: the class algebraic multigrid was designed for. Aggregation AMG
  on graph Laplacians converges in 10–20 iterations independently of size and
  the coarse operators are themselves weighted graph Laplacians, so every
  level uses the same edge kernel as the fine level.

## Mathematical structure to exploit

Each level `ℓ` is a weighted graph on `n_ℓ` free "nodes" (aggregates):

* edges `(u, v, w)` between free nodes, `w = Σ q` over fine edges joining the
  two aggregates;
* an anchor weight `d_u ≥ 0` per node, `d_u = Σ q` over fine edges from the
  aggregate to fixed nodes;
* `(A_ℓ x)_u = d_u x_u + Σ_{(u,v)} w (x_u − x_v)`.

With a piecewise-constant prolongation `P` (each fine node belongs to one
aggregate), the Galerkin coarse operator `Pᵀ A_ℓ P` is exactly the quotient
graph with summed weights: intra-aggregate edges vanish, inter-aggregate
edges add up, anchor weights add up. So

* the **aggregation** (which fine node goes to which aggregate) depends only
  on topology and on the *pattern* of strong connections, and is built once
  per topology (or once per solve if we want it to follow `q`);
* the **coarse weights** depend on `q` and are recomputed per evaluation in
  one `O(nnz)` pass per level — no sparse matrix products;
* every level's operator, residual and smoother is the same matrix-free edge
  kernel used at the fine level, on three right-hand sides at once.

Prolongation/restriction are gather/scatter through the aggregate index array
(`x_fine[i] += x_coarse[agg[i]]`, `r_coarse[agg[i]] += r_fine[i]`).

## Algorithm choices

* **Aggregation**: pairwise heavy-edge matching applied twice (aggregates of
  ~4 nodes, coarsening factor ~4 per level), following Notay's AGMG. Match a
  node with the neighbour of largest `w / sqrt(d_u d_v)`-type strength so
  that anisotropic force densities (bounds allow 100×) coarsen along strong
  connections. Stop when the level has fewer than ~2,000 nodes or coarsening
  stalls (`n_{ℓ+1} > 0.7 n_ℓ`). Expect 7–8 levels at 5M nodes.
* **Smoother**: Chebyshev polynomial smoothing (degree 2–3) with Jacobi
  scaling. It is made only of operator applications and vector updates, so
  it reuses the edge kernel, parallelises trivially and is deterministic;
  Gauss-Seidel would converge slightly faster per sweep but is sequential.
  The spectral bound comes from a few power iterations at setup, refreshed
  per evaluation only if `q` changed by more than a threshold (Gershgorin as
  a fallback).
* **Cycle**: plain aggregation with a V-cycle loses grid independence as
  levels are added. Use the K-cycle (two flexible-CG iterations at each
  coarse level, recursively), which is what makes plain aggregation robust,
  and wrap it in **flexible CG** (FCG, Notay) since the K-cycle is a
  nonlinear preconditioner. Keep a plain V-cycle + standard PCG as an option
  for the comparison and for small sizes.
* **Coarsest level**: faer LLT via the existing `Factorization` on a ≤2,000
  node operator (< 1 ms).
* **Right-hand sides**: x, y, z solved as a block sharing every operator
  application, with independent scalars per column and per-column
  convergence (x/y typically converge much faster than z).
* **Warm start**: forward solve from the previous free-node positions,
  adjoint from the previous multipliers. During the later part of an
  optimisation this cuts iterations by 2–3×.
* **Tolerance**: relative residual, tightened as the optimizer converges:
  `tol_k = clamp(c · ‖proj grad_k‖ / ‖proj grad_0‖, 1e-10, 1e-6)`. Line
  search parameters (`ftol = 1e-4`) tolerate the resulting inconsistency
  early on; near convergence the tolerance floor keeps loss differences of
  `1e-8` relative meaningful. The termination check keeps its unit floor.

## Integration in Theseus

* New module `crates/theseus/src/amg/` with `graph.rs` (level graphs and the
  edge kernel), `aggregate.rs`, `smoother.rs`, `cycle.rs`, `fcg.rs`, and
  `mod.rs` exposing `AmgSolver { setup(topology), update(q), solve(rhs, x0,
  tol, cancel) }`.
* `FdmCache.linear_solver: LinearSolver` enum: `Direct(Factorization)` or
  `Amg(AmgSolver)`, chosen by `SolverOptions.linear_solver: Auto | Direct |
  Amg`. `Auto` picks `Amg` above ~300k free nodes when the bounds guarantee
  `q > 0` (`FactorizationStrategy::Cholesky`), otherwise `Direct`; the LDLᵀ
  path for possibly-negative `q` stays direct (CG requires SPD).
* `fdm::factor_and_solve` and `gradients::solve_adjoint` dispatch on the
  enum. The self-weight/pressure Newton iteration already uses an operator
  form (`apply_a_xyz`, GMRES with `apply_a_preconditioner`); its
  preconditioner becomes one AMG cycle in the `Amg` case.
* The assembled `a_matrix`/`q_to_nz` are not built in the `Amg` case (saves
  ~30 MB per 1M edges and the `CnᵀCn` pattern setup).
* `inverse.rs`, `nullspace.rs`, `inverse_extra.rs` keep using the direct
  factorization; they are separate features with different systems.
* Cancellation: CG is interruptible between iterations, which is better than
  the current "one factorization is one non-interruptible boundary".
* FFI/C#: one new option field; no ABI change to results.

## The O(ne) loops must go parallel too

At 1M edges the non-solver loops (assembly, geometry, objective, gradient
accumulation) are 44 ms single-threaded; at 10M they would be ~0.5 s per
evaluation, comparable to the solve. Convert them to a node-gather form over
a CSR adjacency (each node sums its incident edges) so they can run under
rayon in fixed-size chunks with no per-split allocation and a deterministic
reduction order. The edge kernel of the AMG uses the same adjacency, so this
is one data structure. Objectives that reduce over nodes (`TargetXYZ`, etc.)
get chunked parallel sums with a fixed summation order.

## Budgets at 10M edges (5M free nodes, 4 cores)

Memory, dominated by 15M-element `f64` vectors (120 MB each):

| item | MB |
|---|---:|
| edge arrays (endpoints u32, q, lengths, forces, grad_q, bounds) | ~800 |
| CSR adjacency (20M entries u32 + offsets) | ~100 |
| CG/FCG vectors for 3 rhs (x, r, z, p, Ap) + level temporaries (~2× fine) | ~1,000 |
| hierarchy (coarse graphs, ~1/3 of fine edges total) | ~150 |
| total | ~2,000 |

For comparison the direct path at 10M would need ~2.2 GB for the factor
alone plus ~1 GB of construction transients, and 20–40 s per evaluation.

Time per evaluation, extrapolated from the measured 0.4 ms matvec per 100k
edges and assuming ~3× from 4 threads on memory-bound kernels:

| | 1M edges | 10M edges |
|---|---:|---:|
| operator application, 3 rhs | ~1.5 ms | ~15–20 ms |
| K-cycle (≈ 10 operator-equivalents) | ~15 ms | ~200 ms |
| FCG iteration | ~17 ms | ~220 ms |
| solve, cold (≈ 15 iterations) | ~0.25 s | ~3.3 s |
| solve, warm (≈ 5 iterations) | ~0.09 s | ~1.1 s |
| evaluation (forward + adjoint + loops), warm | ~0.2 s | ~2.5 s |
| evaluation today, direct, 1 thread | 0.88 s | not feasible |

So at 1M edges expect ~3–4× over the current direct path with 4 threads
(roughly parity single-threaded), and at 10M edges a 100-iteration
optimisation in the 5–10 minute range. The AMG buys feasibility and linear
scaling, not interactivity, at 10M.

## Phases, deliverables and exit criteria

**Phase 0 — de-risking prototype (in the `matrix_free_pcg` harness).**
Plain pairwise aggregation, Jacobi/Chebyshev smoother, V-cycle and K-cycle,
PCG and FCG, single-threaded, grid fixture only. Exit: iteration count to
`1e-8` relative residual on grids 64–708 is flat within ±30% for the K-cycle
(expect 10–20), and wall time at 708 beats the direct refactor + two solves
(757 ms) single-threaded. If the V-cycle is already flat, keep it and drop
FCG. If neither is flat, switch to smoothed aggregation before writing any
production code.

**Phase 1 — `amg` module with tests.** Level graphs, aggregation, coarse
weight update, Chebyshev smoother with spectral estimate, cycle, FCG, coarsest
direct solve, block right-hand sides, warm start, cancellation. Tests:
coarse operator equals explicit `Pᵀ A P` on small graphs; V-cycle is
symmetric positive definite (`⟨Mu, v⟩ = ⟨u, Mv⟩`) so PCG is valid; solution
matches the direct solve to tolerance on grids, on an irregular triangulated
mesh, on a 3D-like cable dome, and on `q` fields with 100× and 1e4× ratios;
solve of a system with two connected components each with one support.

**Phase 2 — integration.** `LinearSolver` enum, `SolverOptions.linear_solver`,
dispatch in `factor_and_solve`/`solve_adjoint`/the load Newton iteration,
adaptive tolerance in `optimizer.rs`, skip assembled-matrix construction in
the AMG case. Tests: `optimizer_contract.rs` cases pass with `Amg` forced;
gradient check against finite differences at tight tolerance; a solve with
`Amg` and one with `Direct` reach the same loss to `1e-6` relative on the
recoverable grid; `bench_scale` with `THESEUS_LINEAR_SOLVER=amg` at 224–708.

**Phase 3 — parallel O(ne) loops.** CSR adjacency shared with the AMG,
node-gather geometry/gradient/assembly under rayon with chunked deterministic
reductions, parallel edge kernel and vector operations. Exit: identical
results for 1 and 4 threads (bitwise) and ≥2.5× on the loops and operator at
1M edges on 4 cores.

**Phase 4 — 10M.** Memory audit against the budget table (peak RSS via the
benchmark), 2237×2237 grid end to end, FFI marshalling and progress
reporting at that size, tuning of aggregation depth, smoother degree and
the `Auto` threshold. Exit: `bench_scale` at 2237 completes 20 iterations
within the time and memory budgets above.

## Risks and how the plan handles them

* **Convergence on real networks** (irregular meshes, mixed member types,
  very few supports). Aggregation AMG is robust on M-matrices but iteration
  counts do move with structure; Phase 1 tests non-grid fixtures, and the
  direct solver remains one option away.
* **Extreme `q` ratios.** Strength-based matching handles the 100× the box
  bounds allow; beyond ~1e4 expect more iterations. Report iteration counts
  in the termination string so users can see it.
* **Indefinite systems** (`q` allowed negative in soft-bounds mode). Not
  supported by CG; `Auto` selects the direct LDLᵀ path, as today.
* **Inexact gradients.** Adaptive tolerance plus the loss floor in the
  termination check; Phase 2 verifies the line search does not stall and the
  final loss matches the direct path.
* **Nonlinear preconditioner.** K-cycle requires FCG rather than CG; if the
  Phase 0 V-cycle is already grid-independent, standard PCG suffices and the
  code is simpler.
* **Determinism under threads.** Node-gather kernels and fixed-order chunked
  reductions everywhere; no rayon `fold` with per-split allocation.
