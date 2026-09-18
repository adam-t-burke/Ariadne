# Iterative solver program: matrix-free AMG, CPU and GPU, 10M+ edges

Program-level design and execution plan for adding a matrix-free,
multigrid-preconditioned linear solver to Theseus, on CPU and then on GPU
(Windows and macOS), as an **explicit, opt-in** alternative to the sparse
Cholesky solve, together with the test and benchmark infrastructure needed
to know where each solver wins. The document is written so that independent
agents can each take a workstream, build against frozen interfaces, and hand
back tested, benchmarked increments.

Companion documents: `AMG_PLAN.md` (algorithmic rationale and budgets),
`BENCHMARKS.md` (measurements this program starts from).

---

## 0. Ground rules

1. **Explicit toggle only.** The iterative backends are never selected
   automatically in the Grasshopper plugin. `LinearSolverKind::Direct` is the
   default at every layer (Rust `SolverOptions`, FFI, C#, component
   persistence). An `Auto` mode is *designed for* (the crossover search
   produces the data it would need) but not exposed until this program is
   complete and a separate decision is made.
2. **Loud failure, no silent fallback.** If the user selects a GPU backend
   and no usable adapter exists, or selects an iterative backend for a
   system it cannot solve (indefinite `q`), the solve returns a typed error
   with an actionable message. The user switches the toggle; the code does
   not switch for them.
3. **The direct path is untouched in behaviour.** All existing tests keep
   passing bit-for-bit on `Direct`. Shared infrastructure changes (CSR
   adjacency, parallel loops) must reproduce current results exactly or
   within documented floating-point reassociation, verified by tests.
4. **Determinism.** For a given backend and thread count the result is
   reproducible; for CPU backends the result is identical across thread
   counts. GPU results agree with CPU to solve tolerance (bitwise equality
   across backends is not a goal).
5. **Interfaces first.** Traits and data types in §3 are landed as a stub
   commit before parallel workstreams start, so agents build against the
   same signatures. Changing an interface requires updating this document.
6. **Every increment is measured.** No workstream is done without its
   benchmark numbers recorded in `benchmarks/results/` (§6) and its tests
   running in CI on Linux, Windows and macOS.
7. **Rust 1.89 toolchain** (`RUSTUP_TOOLCHAIN=1.89.0` here; MSRV 1.87 from
   the Basin PR). Python tooling uses `uv`.

---

## 1. Baseline and target

Measured on this branch (4 vCPU Linux, single thread, grid fixture,
`Direct`):

| edges | eval (ms) | factorization share | per L-BFGS-B iteration | peak RSS |
|------:|----------:|--------------------:|-----------------------:|---------:|
| 100k | 49 | 79% | 80 ms | 110 MB |
| 1M | 877 | 75% | 1.48 s | 1.1 GB |
| 10M | ~20–40 s (extrapolated) | | | ~3 GB+ |

Targets (from `AMG_PLAN.md`, to be confirmed by measurement):

| backend | 1M edges / eval | 10M edges / eval | memory at 10M |
|---|---:|---:|---:|
| `IterativeCpu`, 4 cores | ~0.2 s | ~2.5 s | ~2 GB |
| `IterativeGpu`, mid-range discrete GPU | ~30–50 ms | ~0.3–0.5 s | ~2 GB host + ~1.5 GB device |

Acceptance for the program (§8): `IterativeGpu` completes a 20-iteration
solve on the 2237×2237 grid (10.0M edges) on a Windows/NVIDIA machine and on
an Apple-silicon Mac within 1.5× of the 10M target; the toggle works in
Grasshopper; the crossover report exists for at least three machines.

---

## 2. Architecture

### 2.1 Solver selection model

```rust
/// Which linear solver evaluates A(q) x = b and the adjoint system.
/// Serialised as i32 across FFI; values are stable.
#[repr(i32)]
pub enum LinearSolverKind {
    Direct = 0,        // sparse Cholesky/LDLT via faer (today's path, default)
    IterativeCpu = 1,  // matrix-free FCG + aggregation AMG, CPU backend
    IterativeGpu = 2,  // same algorithm, preconditioner (and where possible the
                       // outer iteration) on a wgpu device
}

pub struct IterativeSolverOptions {
    pub tolerance: TolerancePolicy,         // Fixed(f64) | Adaptive { floor, ceiling, factor }
    pub max_iterations: u32,                // per solve, default 200
    pub cycle: CycleKind,                   // V | K  (default K)
    pub smoother_degree: u8,                // Chebyshev degree, default 2
    pub aggregation_passes: u8,             // pairwise matching passes per level, default 2
    pub coarsest_size: u32,                 // stop coarsening below this many nodes, default 2000
    pub precondition_precision: Option<Precision>, // None = backend default (F64 CPU, F32 GPU);
                                                   // FFI maps -1 -> None
    pub gpu: GpuOptions,                    // adapter preference, memory cap, outer-loop placement
}

pub struct SolverOptions {
    // ... existing fields unchanged ...
    pub linear_solver: LinearSolverKind,           // default Direct
    pub iterative: IterativeSolverOptions,         // ignored when Direct
}
```

`Direct` is the default in `SolverOptions::default()`, in the FFI handle's
initial state, in `TheseusInterop.cs`, and in `OptConfigComponent`'s
persisted state (missing key ⇒ `Direct`).

### 2.2 Module layout (Rust, `crates/theseus/src`)

```
linear_solver/
  mod.rs          LinearSolver enum + dispatch, LinearSolverKind, options, SolveStats, errors
  direct.rs       thin adapter over types::Factorization (existing code moves here unchanged)
  tolerance.rs    TolerancePolicy and the per-evaluation tolerance schedule
graph/
  mod.rs          CsrAdjacency (node -> incident edges with sign), EdgeList, LevelGraph
  build.rs        construction from NetworkTopology; boundary edges; anchor weights
amg/
  mod.rs          AmgSolver { setup(&LevelGraph), update(&q), precondition(...) }
  aggregate.rs    pairwise heavy-edge matching, aggregate maps, coarsening control
  hierarchy.rs    levels, coarse weight update (Σ q over inter-aggregate edges), transfers
  smoother.rs     Chebyshev + Jacobi scaling, spectral radius estimate
  cycle.rs        V-cycle, K-cycle
  coarsest.rs     dense/faer LLT on the coarsest level
  fcg.rs          flexible CG, block of 3 right-hand sides, per-column convergence
backend/
  mod.rs          trait Backend (kernel set + buffers), enum BackendHandle { Cpu, Gpu }
  cpu.rs          rayon kernels, chunked deterministic reductions
  gpu/
    mod.rs        wgpu device/adapter probing, buffer pools, command batching
    kernels.wgsl  edge kernel (gather form), axpy/scale/dot, restrict/prolong, Chebyshev step
    f64.rs        optional SHADER_F64 outer loop where the adapter supports it
factor_solve.rs   (existing) fast triangular solves for the direct path
fdm.rs, gradients.rs   dispatch through linear_solver; parallel O(ne) loops over CsrAdjacency
```

Ownership of these paths per workstream is in §7.4; nobody edits outside
their owned paths without coordinating through the integrator.

### 2.3 Frozen interfaces

```rust
/// Block of K right-hand sides / solutions, row-major: v[node * K + k].
pub struct BlockVec<const K: usize> { pub data: Vec<f64>, pub n: usize }

/// Node-centred adjacency. For node i, incident edges are
/// edges[offsets[i]..offsets[i+1]] with sign[j] = +1 if i is the edge end, -1 if start.
/// Built once per topology; shared by AMG levels, geometry, gradients.
pub struct CsrAdjacency { pub offsets: Vec<u32>, pub edges: Vec<u32>, pub sign: Vec<i8>, pub other: Vec<u32> }

/// One multigrid level as a weighted graph on free nodes with anchor weights.
pub struct LevelGraph {
    pub n: usize,
    pub adjacency: CsrAdjacency,   // over the level's edges
    pub weight: Vec<f64>,          // per edge: Σ q of fine edges it represents
    pub anchor: Vec<f64>,          // per node: Σ q of fine edges to fixed nodes
    pub aggregate_of: Vec<u32>,    // fine node -> coarse node (empty on the coarsest)
}

pub struct SolveRequest<'a> {
    pub rhs: &'a [f64],            // n * 3
    pub x0: Option<&'a [f64]>,     // warm start, n * 3
    pub tolerance: f64,            // relative residual per column
    pub max_iterations: u32,
    pub cancel: Option<&'a AtomicBool>,
}

pub struct SolveStats {
    pub iterations: [u32; 3],
    pub relative_residual: [f64; 3],
    pub converged: bool,
    pub setup_ms: f64, pub solve_ms: f64,
    pub backend: LinearSolverKind,
}

/// `Send` so a boxed solver inside `FdmCache` can move to a worker thread.
pub trait LinearSystemSolver: Send {
    /// Called whenever q changed. Direct: refactor. Iterative: recompute level weights,
    /// spectral bounds if needed, upload to device.
    fn update(&mut self, q: &[f64]) -> Result<(), TheseusError>;
    fn solve(&mut self, req: SolveRequest<'_>, x: &mut [f64]) -> Result<SolveStats, TheseusError>;
    fn kind(&self) -> LinearSolverKind;
    fn memory_bytes(&self) -> MemoryReport;   // host, device
    // Defaulted (WS-D): preconditioner application for the load Newton GMRES
    // (AMG overrides with one V-cycle), and downcasts used by FdmCache to
    // borrow the assembled matrix on the Direct path.
    fn precondition(&mut self, req: SolveRequest<'_>, x: &mut [f64]) -> Result<SolveStats, TheseusError> { self.solve(req, x) }
    fn as_direct(&self) -> Option<&DirectSolver> { None }
    fn as_direct_mut(&mut self) -> Option<&mut DirectSolver> { None }
}

/// Factory (unit struct, not an enum): `LinearSolver::new(kind, &NetworkTopology, &Bounds,
/// &IterativeSolverOptions) -> Result<Box<dyn LinearSystemSolver>, TheseusError>`.
pub struct LinearSolver;

pub trait Backend {
    type Buf;                                   // device or host vector (block of 3 per node)
    type LevelGraphBuf;                         // one uploaded LevelGraph (CSR + weights + Jacobi diagonal)
    type AggBuf;                                // one uploaded aggregate map
    fn handle(&self) -> BackendHandle;
    fn alloc(&self, len: usize, precision: Precision) -> Self::Buf;
    fn len(&self, buf: &Self::Buf) -> usize;
    fn upload(&self, src: &[f64], dst: &mut Self::Buf);
    fn download(&self, src: &Self::Buf, dst: &mut [f64]);
    fn copy(&self, src: &Self::Buf, dst: &mut Self::Buf);
    fn zero(&self, buf: &mut Self::Buf);
    fn upload_level(&self, level: &LevelGraph, precision: Precision) -> Self::LevelGraphBuf;
    fn update_level_weights(&self, level: &LevelGraph, dst: &mut Self::LevelGraphBuf); // per update(q)
    fn upload_aggregates(&self, aggregate_of: &[u32], n_coarse: usize) -> Self::AggBuf;
    fn apply_graph(&self, level: &Self::LevelGraphBuf, x: &Self::Buf, y: &mut Self::Buf);      // y = A_l x
    fn residual(&self, level: &Self::LevelGraphBuf, b: &Self::Buf, x: &Self::Buf, r: &mut Self::Buf); // r = b - A_l x
    /// One Chebyshev update with shared scalars (same matrix for all 3 columns) using the
    /// level's Jacobi diagonal: d = alpha * D⁻¹ r + beta * d; x += d.
    fn chebyshev_step(&self, level: &Self::LevelGraphBuf, alpha: f64, beta: f64,
                      r: &Self::Buf, d: &mut Self::Buf, x: &mut Self::Buf);
    fn restrict(&self, agg: &Self::AggBuf, fine: &Self::Buf, coarse: &mut Self::Buf);
    fn prolong_add(&self, agg: &Self::AggBuf, coarse: &Self::Buf, fine: &mut Self::Buf);
    fn axpy(&self, alpha: [f64; 3], x: &Self::Buf, y: &mut Self::Buf);                   // per-column scalars
    fn scale(&self, alpha: [f64; 3], x: &mut Self::Buf);
    fn dot3(&self, a: &Self::Buf, b: &Self::Buf) -> [f64; 3];
    fn norm3(&self, a: &Self::Buf) -> [f64; 3];
    fn sync(&self);
}
```

Conventions fixed by WS-I: the `Direct` implementation reports
`SolveStats { iterations: [0; 3], relative_residual: [0.0; 3], .. }` (not
measured) and `setup_ms` = duration of the preceding `update`; the diagonal
perturbation used by `factor_and_solve` is solver state
(`DirectSolver::set_perturbation`) because `update(q)` carries no
perturbation argument — iterative solvers apply it as an anchor-weight
shift set the same way. `DirectSolver` exposes `a_matrix()` and
`factorization()` accessors; WS-D decides whether `FdmCache` keeps its own
`a_matrix` for `apply_a_xyz`/gradient loops or reads it from the solver
(recommended: the solver is the single owner, and `FdmCache` borrows it).

Errors (added to `TheseusError`):

```rust
IterativeSolverDidNotConverge { iterations: u32, relative_residual: f64, kind: LinearSolverKind },
IterativeSolverUnsupported(String),   // e.g. "q may be negative under these bounds; iterative solvers need q > 0"
GpuUnavailable(String),               // adapter probe result, with the adapters seen
GpuOutOfMemory { requested: u64, available: u64 },
```

FFI codes (`ffi.rs`, `error_code(&TheseusError)`): all pre-existing
variants keep `-1`, panics `-2`; `IterativeSolverDidNotConverge` `-3`,
`IterativeSolverUnsupported` `-4`, `GpuUnavailable` `-5`, `GpuOutOfMemory`
`-6`. `theseus_last_error` carries the message. New codes are appended, never
renumbered.

### 2.4 FDM integration points

* `FdmCache` gains `linear_solver: Box<dyn LinearSystemSolver>` and
  `adjacency: CsrAdjacency`; `a_matrix`/`q_to_nz` are built only for
  `Direct`.
* `fdm::factor_and_solve` → `linear_solver.update(q)` +
  `solve(rhs, warm = previous x)`; `gradients::solve_adjoint` → `solve(grad_x,
  warm = previous λ)`. The perturbation argument applies as an anchor-weight
  shift in the iterative case.
* The self-weight/pressure Newton iteration uses `apply_a_xyz` and GMRES with
  `apply_a_preconditioner`: the preconditioner becomes one AMG cycle when
  iterative; the operator uses the level-0 graph.
* Tolerance schedule: `optimizer.rs` computes `tol_k` per evaluation from
  the projected-gradient ratio (`AMG_PLAN.md` §"Algorithm choices") and
  passes it through `FdmProblem` to the cache; `Fixed` policy bypasses it.
* `SolveStats` accumulate into `SolverResult` (`linear_solver_iterations:
  Vec<u32>` per evaluation, and totals in `termination_reason`).
* Cancellation: checked between FCG iterations and between Newton
  iterations.

### 2.5 GPU design (wgpu)

* **Portability baseline**: wgpu with WGSL, `f32` storage buffers. Runs on
  DX12/Vulkan (Windows), Metal (macOS), Vulkan (Linux). No CUDA, no
  vendor SDK redistribution; the plugin's `cdylib` links wgpu statically.
* **Precision**: the preconditioner (all AMG levels, smoothers, transfers)
  runs in `f32`. The outer FCG runs in `f64`: on the GPU when the adapter
  reports `Features::SHADER_F64` (typically Vulkan/DX12 on discrete GPUs),
  otherwise on the CPU with the `f32` preconditioner on the device
  (`gpu.outer_loop: Auto | Device | Host`). The residual used for
  convergence is always `f64`.
* **Kernels** (`kernels.wgsl`), one thread per node or per element:
  `apply_graph` (gather over CSR incident edges), `chebyshev_step`,
  `restrict`, `prolong_add`, `axpy`, `scale`, `dot_partial` (workgroup
  reductions to a partials buffer, finished in fixed order on the host or in
  a second pass), `residual`, `coarse_weight_update` (per coarse edge, sum
  over its fine-edge list). Coarsest-level solve stays on the host (faer
  LLT on ≤ 2,000 nodes; the vector is 6 KB–48 KB).
* **Buffers**: allocated once per topology in a `BufferPool`; per-`q`
  updates upload only `weight` arrays; no allocation inside a solve. Host
  staging uses mapped buffers with double buffering. Memory report includes
  device bytes.
* **Command batching**: one command encoder per K-cycle; norms read back
  once per FCG iteration through a small staging buffer (or every 2
  iterations if the readback latency dominates at small sizes; measured in
  WS6).
* **Adapter selection**: enumerate adapters; prefer `DiscreteGpu`, then
  `IntegratedGpu`; reject `Cpu`-type (software) adapters unless
  `THESEUS_GPU_ALLOW_SOFTWARE=1` (used for CI correctness runs with
  lavapipe/SwiftShader). Report the chosen adapter in `SolveStats` and in
  the FFI probe. Require `max_storage_buffer_binding_size` ≥ the largest
  vector; split kernels over multiple bindings if the adapter limit is
  below 2 GB.
* **Failure modes**: adapter absent → `GpuUnavailable`; allocation failure →
  `GpuOutOfMemory`; device lost mid-solve → `GpuUnavailable` with the
  message, no retry.
* **macOS specifics**: unified memory, so host/device copies are cheap;
  Metal has no `f64`, so the outer loop runs on the host there. Validate on
  Apple silicon (M-series) and, if available, an Intel Mac with AMD GPU.
* **Windows specifics**: prefer DX12 for driver robustness, allow Vulkan via
  `WGPU_BACKEND`. Validate on NVIDIA (discrete, `SHADER_F64`), AMD, and an
  Intel iGPU (integrated, likely no `f64`).

### 2.6 FFI, C#, Grasshopper

* FFI (`ffi.rs`), following the `theseus_set_q_parameterization_mode`
  pattern:
  * `theseus_set_linear_solver(handle, kind: i32) -> i32`
  * `theseus_set_iterative_options(handle, tolerance_mode: i32, tolerance: f64, max_iterations: u32, cycle: i32, smoother_degree: u32, aggregation_passes: u32, coarsest_size: u32, precondition_precision: i32, gpu_outer_loop: i32) -> i32`
  * `theseus_gpu_probe(out_json: *mut c_char, cap: usize) -> i32` — adapters
    seen (name, backend, type, limits, `SHADER_F64`), chosen adapter.
  * `theseus_get_linear_solver_stats(handle, out: *mut TheseusLinearSolverStats) -> i32`
    (totals: solves, iterations, time, backend actually used).
* C# (`Theseus/TheseusInterop.cs`, `TheseusSolver.cs`): `enum
  LinearSolverKind { Direct, IterativeCpu, IterativeGpu }`, options record,
  probe wrapper, stats in the result model; default `Direct`.
* Grasshopper (`Solver/Components/OptConfigComponent.cs`): context-menu
  group "Linear solver" with three radio items, persisted under key
  `LinearSolverKind` (missing ⇒ `Direct`), nickname suffix `lin: Direct |
  Iter-CPU | Iter-GPU`; an "Advanced iterative options…" sub-menu or a
  separate `IterativeSolverOptions` component feeding the config. Solve
  component shows the backend used and total linear-solver iterations in its
  runtime messages; a GPU probe failure surfaces as a component error
  (red), not a warning.
* Docs: `README.md` section "Linear solver backends" (what each does, when
  to try it, that it is opt-in), `Theseus/NATIVE_ARTIFACTS.md` (wgpu adds
  no runtime files but note driver requirements).

---

## 3. Algorithms (normative summary; rationale in `AMG_PLAN.md`)

**Revised after Phase 0 (WS-A, `BENCHMARKS.md` "Phase 0").** Plain
(unsmoothed) aggregation is a no-go: V-cycle iterations grow 3.5–4× from
8k to 1M edges on every fixture, and the K-cycle that keeps them flat costs
8–15 operator applications per iteration (2.7–3.8× slower than the direct
solve at 1M). Smoothed aggregation with a plain V-cycle inside PCG is the
algorithm this program builds; the K-cycle and FCG are dropped.

* **Fine operator (level 0)**: `(A x)_u = anchor_u x_u + Σ_{e ∋ u} q_e (x_u −
  x_other)` in gather form over `CsrAdjacency`, matrix-free. Fixed nodes
  never appear in vectors; edges to them contribute to `anchor` and to the
  right-hand side.
* **Aggregation (tentative prolongator P₀)**: pairwise matching on strength
  `s_uv = w_uv / sqrt((anchor_u + Σ_e w_ue)(anchor_v + Σ_e w_ve))`; greedy
  in decreasing node degree order with a deterministic tie-break (lowest
  index); unmatched nodes form singletons or join the strongest neighbour's
  pair if it has < 3 members; **three passes** (aggregates of ≤ 8; 4 levels
  at 500k free nodes). Stop when `n_{l+1} < coarsest_size` or `n_{l+1} > 0.7
  n_l`. Aggregation uses the `q` current at setup and is rebuilt only when
  topology changes or on explicit drift re-setup (below).
* **Smoothed prolongator**: `P = (I − ω D⁻¹ A) P₀` with `ω = 4 / (3 λ_max)`
  (Jacobi-smoothed aggregation). Coarse operators `A_{l+1} = Pᵀ A_l P` are
  **general sparse SPD matrices**, not graph Laplacians: levels ≥ 1 are
  stored as block-of-3-compatible CSR (`LevelMatrix { row_ptr, col_idx,
  values, diag }`); level 0 stays matrix-free. Coarse-edge/column ordering
  within a row is ascending column index (the convention WS-G's GPU tests
  assume).
* **Update on `q` change** (the cost that decides CPU viability; the
  prototype's naive update took 320–380 ms at 1M edges): the symbolic
  pattern of every `P` and `A_l` is frozen at setup; `update(q)` refills
  level-0 weights/anchors (O(ne), `Level0Map::update_weights`) and reruns
  only the **numeric** triple products `Pᵀ A_l P` level by level with the
  frozen pattern (row-parallel sparse accumulation, deterministic order).
  `P`'s values are also frozen (computed from the setup `q`): the
  prolongator then no longer tracks `q`, which is admissible because `P`
  only needs to approximate the near-nullspace (constants), which is
  `q`-independent. Re-setup (new `P`, new pattern) when `max_e |Δq_e /
  q_e^{setup}| > 2` or when iterations exceed 2× the setup count twice in a
  row. Target: `update(q)` ≤ 100 ms single-threaded at 1M edges, ≤ 40 ms on
  4 threads (WS-C acceptance).
* **Smoother**: Chebyshev **degree 2** with Jacobi scaling `D⁻¹`, spectral
  interval `[λ_max/α, λ_max]` with **α = 10** (α = 30 costs +40%
  iterations; degree 3 with α = 30 loses size-independence on the dome).
  `λ_max` per level from 10 power iterations at setup; re-estimated when
  `max_e |Δq_e / q_e| > 0.5` since the last estimate, else scaled by the
  ratio of max weights (cheap upper bound).
* **Cycle**: V-cycle (one pre- and one post-smoothing sweep) as a fixed
  SPD preconditioner. Coarsest level: faer LLT via `Factorization`,
  refactored on `update(q)`; `coarsest_size` default 2,000, with 5–10k plus
  sparse Cholesky to be tried against the grid z-column growth (+47%
  between 64 and 448, flat beyond) in WS-C/WS-J.
* **Outer**: standard PCG on the block of 3 columns, each column with its
  own scalars; stop each column at `‖r_k‖ ≤ tol ‖b‖` (with `‖b‖` floored
  at `1e-300`); a column with `b = 0` returns `x = 0`. The `CycleKind::K`
  option remains in `IterativeSolverOptions` for experiments but is not the
  default and may be removed by WS-K.
* **Warm starts**: forward from the previous `x` of the same cache; adjoint
  from the previous `λ`. First solve of a session cold. Measured gain
  1.5–2× on x/y, 1.3–1.4× on the z column that sets the time.
* **Expected standing versus `Direct` (single thread, 1M edges, from
  Phase 0)**: SA-PCG solve 1.0–1.55× direct cold, 0.75–1.2× warm, before
  the update cost. A CPU win therefore depends on (i) the frozen-pattern
  numeric update above, (ii) 4-thread scaling of the memory-bound kernels
  (WS-B infrastructure), and (iii) the adjoint reusing the same hierarchy.
  The GPU backend is where the large factor is expected; the crossover
  search (§5.3) decides where each pays and the toggle stays explicit.
* **Tolerance policy**: `Adaptive { floor = 1e-10, ceiling = 1e-6, factor =
  1e-2 }`: `tol_k = clamp(factor · ‖g⁺_k‖ / ‖g⁺_0‖, floor, ceiling)` where
  `g⁺` is the projected gradient (bounded) or gradient (unbounded) from the
  previous accepted iterate; the first evaluation uses `ceiling`.
  `Fixed(tol)` overrides.
* **Determinism (CPU)**: node-gather kernels; dot products by fixed 4,096-
  element chunks summed in index order (parallel over chunks, sequential
  combine); no `rayon::fold` with per-split state.

---

## 4. Testing strategy

Levels, all mandatory unless marked optional:

1. **Unit (per module)** — `graph`: adjacency round-trips edge lists,
   boundary edges and anchors match the incidence definition; `amg`: coarse
   operator equals explicit `Pᵀ A P` on random small graphs; Chebyshev
   damps the high end of the spectrum on a 1-D chain; V-cycle symmetric and
   positive definite (`⟨M u, v⟩ = ⟨u, M v⟩`, `⟨M u, u⟩ > 0`) so PCG is
   valid; K-cycle + FCG converge on grids 16–256 with iteration counts flat
   within ±30%; `fcg`: solves an SPD dense system to 1e-12; `backend`: every
   kernel checked CPU vs reference scalar implementation.
2. **Cross-backend equivalence** — GPU kernels vs CPU on random inputs to
   `1e-6` relative (f32) / `1e-12` (f64); full solves agree to the requested
   tolerance; runs in CI under a software adapter
   (`THESEUS_GPU_ALLOW_SOFTWARE=1`, lavapipe on Linux, WARP on Windows,
   Metal on macOS runners) and on real GPUs manually.
3. **Integration** — `optimizer_contract.rs` and `optimization_diagnostic.rs`
   executed with each `LinearSolverKind` forced (parameterised); the
   recoverable grid reaches the same loss as `Direct` to `1e-6` relative;
   gradient vs central finite differences at `Fixed(1e-12)`; self-weight and
   pressure fixtures with the iterative preconditioner in the load Newton
   iteration.
4. **Regression fixtures** (`tests/support/fixtures/`): grids; triangulated
   irregular mesh (jittered points, Delaunay-like connectivity); cable dome
   / spoke-wheel (mixed cable/strut topology, many edges per node);
   two supports only (near-singular); disconnected components each with a
   support; `q` fields with 100× and 1e4× ratios; a 1,000-node fixture with
   every objective type.
5. **Determinism** — CPU results identical for `RAYON_NUM_THREADS ∈ {1, 2,
   4}` (bitwise) for both `Direct` and `IterativeCpu`.
6. **Failure modes** — soft-bounds problem that permits negative `q` with
   `IterativeCpu` ⇒ `IterativeSolverUnsupported`; GPU requested with no
   adapter (`WGPU_ADAPTER_NAME=none` hook) ⇒ `GpuUnavailable`; `max_iterations
   = 1` ⇒ `IterativeSolverDidNotConverge` carrying the residual;
   cancellation flag set inside FCG ⇒ `Cancelled` within one iteration.
7. **FFI/C#** — round-trip of every new setter; `Ariadne.Tests` covers
   `LinearSolverKind` persistence and defaulting; probe JSON parses.
8. **Scale smoke** (ignored, nightly or manual) — 1M-edge grid solve with
   each backend completes and matches the direct loss to `1e-4` relative
   after 10 iterations.

CI (`.github/workflows/ci.yml`): CPU tests on `ubuntu-latest`,
`windows-latest`, `macos-latest`; software-adapter GPU tests on all three;
an opt-in `gpu-hardware` job on self-hosted runners labelled `gpu-nvidia`,
`gpu-amd`, `gpu-apple` when available (skipped otherwise). Clippy and fmt
gates apply to new files.

---

## 5. Benchmark program

### 5.1 Harness

`crates/theseus/tests/bench_scale.rs` is extended (not replaced) to take
`THESEUS_LINEAR_SOLVER`, `THESEUS_FIXTURE`, `THESEUS_SCALE_GRIDS`,
`THESEUS_SCALE_ITERS`, `THESEUS_ITER_*` (AMG parameters) and to emit one JSON
line per run in addition to the table: fixture, size (edges, free nodes),
backend, adapter, threads, parameters, setup ms, eval ms (median of ≥ 5),
evaluations, iterations, linear-solver iterations (forward/adjoint,
mean/max), total ms, peak RSS, device bytes, final loss, git SHA, machine
id. `examples/profile_phases.rs` gains the same backend switch.

`scripts/bench_sweep.py` (uv, `pyproject` in `scripts/`) drives the harness
over a parameter grid, stores raw JSON under
`benchmarks/results/<machine-id>/<date>-<sha>/`, and renders Markdown
tables and SVG plots into `benchmarks/reports/`. Machine id = OS + CPU model
+ GPU adapter + RAM, stored in `machine.json` next to the results.

### 5.2 Protocol

* Release build, `--locked`; CPU governor/power plan noted; laptops on
  mains.
* Warm-up: one full run discarded per configuration.
* Repetitions: ≥ 5 per cell for eval timings, ≥ 3 for full solves; report
  median and interquartile range; flag cells whose IQR exceeds 10% of the
  median.
* Threads: `RAYON_NUM_THREADS ∈ {1, physical cores}`; both recorded.
* Sizes: grids 72, 160, 224, 320, 448, 708, 1000, 1414, 2237 (10k → 10M
  edges); irregular and dome fixtures at matched edge counts.
* Both backends on every size they can complete within a 10-minute wall
  cap and the machine's memory; cells that cannot complete are recorded as
  such (they are data for the crossover).
* Fixed iteration budgets (`tolerances = 0`) for comparability, plus one
  run per size to the default tolerances to record real iteration counts.

### 5.3 Crossover parameter search

Goal: for each machine class and fixture class, find the edge count above
which `IterativeCpu` (and separately `IterativeGpu`) beats `Direct` per
evaluation and per solve, and the AMG parameters that make the iterative
path fastest. Output feeds a future `Auto` mode; it is not wired into
selection in this program.

Stage A — **algorithm parameters** (CPU, three sizes: 224, 708, 1414;
grid + irregular): full factorial over `cycle ∈ {V, K}`, `smoother_degree ∈
{1, 2, 3}`, `aggregation_passes ∈ {1, 2, 3}`, `coarsest_size ∈ {500, 2000,
8000}`, spectral factor `α ∈ {10, 30, 100}` → 162 cells × 2 fixtures × 3
sizes. Metric: time to `1e-8` from cold and from a 2% warm start.
Select the configuration minimising the geometric mean of normalised times;
require iteration counts flat within ±30% across the three sizes, else
discard.

Stage B — **size sweep** with the Stage A configuration: all sizes, both
backends, both thread settings, three fixtures, cold and warm (q steps
0.5%, 2%, 10%). Fit per-machine cost models by least squares on log-log:

* `Direct`: `t = a·n^1.5 + b·n·log n + c` (factor + solves + loops);
* `Iterative`: `t = (d + e·iters(n, step))·n + f`, with `iters` fitted per
  fixture class as `i0 + i1·log n` (should be ~flat).

Crossover = smallest `n` where the fitted iterative time is below direct
for (i) a cold evaluation, (ii) a warm evaluation at 2% step, (iii) a full
40-iteration solve. Report the three numbers with bootstrap 90% intervals
from the repetitions.

Stage C — **validation**: run both backends at `0.5×`, `1×`, `2×` the
predicted crossover; accept if the measured winner matches the prediction
at `0.5×` and `2×`.

Stage D — **GPU tuning** (per GPU): `precondition_precision ∈ {F32, F64
where available}`, `outer_loop ∈ {Device, Host}`, workgroup size `∈ {64,
128, 256}`, norm readback cadence `∈ {1, 2}`; then repeat Stage B for
`IterativeGpu`.

Deliverable: `benchmarks/reports/crossover-<machine-id>.md` and
`crossover.json` (`{fixture_class: {cold, warm_2pct, solve40: {n_edges, ci90}}}`)
for each machine, plus a consolidated table in `BENCHMARKS.md`.

### 5.4 Machines

Minimum: (1) Windows, NVIDIA discrete GPU (`SHADER_F64`); (2) Windows,
integrated Intel or AMD GPU (laptop class, the common Grasshopper user);
(3) macOS Apple silicon; plus the Linux CI VM for CPU regression tracking.
Optional: Windows AMD discrete; Intel Mac with AMD GPU.

---

## 6. Deliverables tree

```
crates/theseus/src/{linear_solver,graph,amg,backend}/…      (§2.2)
crates/theseus/src/ffi.rs                                    new setters, probe, stats
crates/theseus/tests/{linear_solver_*.rs, backend_*.rs, support/fixtures/*.rs}
crates/theseus/tests/bench_scale.rs, examples/profile_phases.rs   backend switch + JSON
scripts/bench_sweep.py, scripts/pyproject.toml, scripts/crossover_fit.py
benchmarks/results/<machine>/…, benchmarks/reports/…
Theseus/TheseusInterop.cs, Theseus/TheseusSolver.cs         enum, options, probe, stats
Solver/Components/OptConfigComponent.cs (+ IterativeSolverOptionsComponent.cs)
Ariadne.Tests/LinearSolverKindTests.cs
.github/workflows/ci.yml                                    software-adapter GPU jobs, hardware job
README.md, Theseus/NATIVE_ARTIFACTS.md, BENCHMARKS.md, this document
```

---

## 7. Workstreams for agent dispatch

Each workstream lists scope, owned paths, dependencies, deliverables,
acceptance criteria and a dispatch prompt. Complexity is described in terms
of components touched, not calendar time.

### 7.1 Dependency graph and waves

```
WS-I  interfaces stub ─┬─► WS-A  Phase-0 prototype (gates algorithm choices)
                       ├─► WS-B  graph + parallel loops
                       ├─► WS-E  benchmark harness + fixtures + sweep tooling
                       ├─► WS-F  FFI + C# + Grasshopper toggle
                       └─► WS-G  GPU kernel spike (kernels are algorithm-agnostic)
WS-A + WS-B ───────────► WS-C  AMG core (CPU)
WS-C + WS-I ───────────► WS-D  FDM integration + tolerance policy + diagnostics
WS-D + WS-G ───────────► WS-H  GPU backend integration + cross-platform validation
WS-D + WS-E ───────────► WS-J  crossover search (CPU), then WS-H + WS-E ─► GPU tuning
WS-D, WS-F, WS-H ──────► WS-K  robustness, docs, CI hardware jobs, release checklist
```

Wave 1 (parallel): WS-I (short, first), then WS-A, WS-B, WS-E, WS-F, WS-G.
Wave 2: WS-C, then WS-D. Wave 3 (parallel): WS-H, WS-J. Wave 4: WS-K.

### 7.2 Branch and merge protocol

* Base branch for the program: `cursor/theseus-scale-100k-2d27` (this
  branch) until it is merged; agents branch as
  `cursor/its-<ws-letter>-<short-name>-2d27`.
* **Private worktree per agent.** Agents that share a machine never work
  in the shared checkout (`/workspace`): the first action is `git worktree
  add <private-dir> -b <branch> <base>`, and all builds run there (a
  private `CARGO_TARGET_DIR` avoids lock contention). Switching HEAD or
  leaving untracked files in the shared checkout breaks the other agents.
* One PR per workstream increment against the program base; the integrator
  (WS-I owner) merges in dependency order and rebases open branches after
  each merge. Interface changes go through the integrator as a separate
  small PR first.
* PR checklist: tests green on three OSes; clippy/fmt clean for owned
  files; benchmark JSON for the affected path attached under
  `benchmarks/results/`; `BENCHMARKS.md` updated if numbers changed;
  this document updated if interfaces or plans changed.

### 7.3 Workstream specifications

#### WS-I — Interfaces and integration (the integrator)

Scope: land §2.1–2.3 as compiling stubs (`linear_solver/`, `graph/mod.rs`
types, `backend/mod.rs` trait, error variants, FFI codes reserved),
`Direct` adapter wrapping the existing `Factorization` behind
`LinearSystemSolver` with all existing tests passing; own this document;
merge and rebase for others.
Owned: `src/linear_solver/mod.rs`, `src/linear_solver/direct.rs`,
`src/backend/mod.rs`, `src/graph/mod.rs` (types only), `TheseusError`
additions, this document.
Acceptance: `cargo test --workspace` unchanged results; `Direct` dispatch
through the trait adds < 1% to the 100k evaluation.

#### WS-A — Phase-0 convergence prototype

Scope: in `examples/matrix_free_pcg.rs` (or a sibling example), implement
pairwise aggregation (1–2 passes), Chebyshev/Jacobi smoother, V- and
K-cycle, PCG and FCG, single-threaded, on the grid, irregular and dome
fixtures. Measure iterations to `1e-8` cold and 2%-warm at grids 64, 128,
224, 448, 708 and matched irregular sizes.
Owned: `examples/matrix_free_pcg.rs`, `examples/amg_prototype.rs`,
`tests/support/fixtures/` (read-only use of WS-E's fixtures; may add
minimal generators if WS-E has not landed, to be reconciled).
Deliverable: table of iterations and wall time per configuration in
`BENCHMARKS.md` ("Phase 0"); a written recommendation: V or K cycle,
smoother degree, aggregation passes; go/no-go for plain aggregation.
Acceptance (go): K-cycle iterations flat within ±30% from 64 to 708 on the
grid and within ±50% on the irregular fixture; wall time at 708 below the
direct refactor + two solves (757 ms) single-threaded. No-go: document,
and prototype smoothed aggregation before WS-C starts.

#### WS-B — Shared graph infrastructure and parallel O(ne) loops

Scope: `CsrAdjacency` construction from `NetworkTopology`; rewrite
`compute_geometry`, `assemble_a` (direct path), `assemble_rhs`,
`accumulate_implicit_gradients`, `accumulate_explicit_gradients`' node loops
and the objective reductions to node-gather form with rayon chunking and
fixed-order reductions; `BlockVec`; deterministic `dot`/`norm`.
Owned: `src/graph/build.rs`, `src/backend/cpu.rs`, the listed functions in
`src/fdm.rs`, `src/gradients.rs`, `src/objectives.rs`.
Acceptance: bitwise-identical results between 1 and 4 threads and identical
to the pre-change sequential code where summation order is unchanged
(documented reassociation otherwise, with tests at `1e-12` relative);
≥ 2.5× speed-up of the non-solver phases at 1M edges on 4 cores; all
existing tests pass; `profile_phases` numbers recorded.

#### WS-C — AMG core (CPU)

Scope: `amg/` per §2.2 and §3 using the `Backend` trait with the CPU
backend; block right-hand sides; warm start; cancellation; `SolveStats`;
unit tests of §4.1; memory report.
Owned: `src/amg/**`, `tests/amg_*.rs`.
Depends on: WS-A recommendation, WS-B `CsrAdjacency` and CPU kernels, WS-I
traits.
Acceptance: unit tests green; solves on all regression fixtures match
`Direct` to the requested tolerance; iteration counts within the WS-A
envelope; no heap allocation inside `solve` after the first call (checked
with the counting allocator from `crates/lbfgsb/tests/allocations.rs`
pattern); single-thread time at 708 ≤ 0.6× the direct path per forward
solve.

#### WS-D — FDM integration, tolerance policy, diagnostics

Scope: `LinearSolver` dispatch in `FdmCache`, `factor_and_solve`,
`solve_adjoint`, the load Newton preconditioner; `tolerance.rs` and the
schedule in `optimizer.rs`; `SolveStats` into `SolverResult` and
termination string; typed errors; `IterativeSolverUnsupported` when bounds
permit `q ≤ 0`; parameterised contract tests; finite-difference gradient
tests; determinism tests.
Owned: `src/linear_solver/tolerance.rs`, dispatch sites in `src/fdm.rs`,
`src/gradients.rs`, `src/optimizer.rs`, `src/types.rs` (`SolverOptions`,
`SolverResult` fields), `tests/linear_solver_*.rs`.
Acceptance: all contract tests pass with each `LinearSolverKind` forced;
recoverable grid reaches the same loss as `Direct` to `1e-6` relative with
`IterativeCpu`; gradient FD check passes at `Fixed(1e-12)`; `bench_scale`
runs with `THESEUS_LINEAR_SOLVER=iterative-cpu` at 224–708 and the numbers
are recorded.

#### WS-E — Benchmark harness, fixtures, sweep tooling

Scope: fixture generators (irregular triangulated mesh, cable dome, few-
support, disconnected, anisotropic `q`) in `tests/support/fixtures/`;
`bench_scale.rs` JSON output and env switches; `profile_phases` backend
switch; `scripts/bench_sweep.py`, `scripts/crossover_fit.py` (uv project,
`numpy`/`scipy` for the fits, `matplotlib` for plots); `benchmarks/`
layout; machine-id capture; report renderer.
Owned: `tests/support/fixtures/**`, `tests/bench_scale.rs`,
`examples/profile_phases.rs`, `scripts/**`, `benchmarks/**`.
Acceptance: a sweep over `Direct` at 72–708 on the CI VM produces results,
a report and fitted cost model with residuals < 15%; fixtures have unit
tests (connectivity, support count, edge count) and are used by WS-A/WS-C.

#### WS-F — FFI, C#, Grasshopper toggle

Scope: §2.6 in full; default `Direct` everywhere; probe; stats; component
UI and persistence; `Ariadne.Tests`; README section. Until WS-D lands,
`IterativeCpu`/`IterativeGpu` return `IterativeSolverUnsupported("not yet
available")` through the same code path so the UI can be built and tested
end to end.
Owned: `src/ffi.rs` (new functions only), `Theseus/*.cs`,
`Solver/Components/OptConfigComponent.cs`,
`Solver/Components/IterativeSolverOptionsComponent.cs`, `Solver/SolverModels.cs`,
`Ariadne.Tests/LinearSolverKindTests.cs`, `README.md` section.
Acceptance: toggling in Grasshopper changes the backend reported in the
result; a saved definition reloads with its selection; a definition saved
before this change loads as `Direct`; selecting GPU on a machine without an
adapter produces a red component error with the probe message.

#### WS-G — GPU kernel spike (wgpu)

Scope: `backend/gpu/` with adapter probing, buffer pool, the WGSL kernel set
of §2.5 in `f32` (and `f64` variants behind `SHADER_F64`), CPU-vs-GPU
kernel equivalence tests runnable under a software adapter, microbenchmarks
of `apply_graph`, `dot3` and readback latency at 1M and 10M nodes.
Owned: `src/backend/gpu/**`, `tests/backend_gpu_*.rs`.
Depends on: WS-I `Backend` trait, WS-B `CsrAdjacency` layout (can start
against the frozen struct).
Acceptance: all kernels match CPU to `1e-6` (f32) / `1e-12` (f64) relative
under lavapipe on Linux CI, WARP on Windows CI and Metal on macOS CI;
`apply_graph` at 10M edges achieves ≥ 50% of the adapter's advertised
bandwidth on at least one real GPU; documented per-adapter limits
encountered.

#### WS-H — GPU backend integration and cross-platform validation

Scope: `AmgSolver` over the GPU backend; mixed precision (`f32`
preconditioner, `f64` outer on device or host); command batching; norm
readback cadence; `IterativeGpu` end to end through FFI and Grasshopper;
validation matrix on the §5.4 machines; per-platform notes.
Owned: `src/backend/gpu/**` (with WS-G), `src/amg/` GPU specialisations,
`docs` per platform.
Acceptance: `IterativeGpu` reaches the same loss as `IterativeCpu` to `1e-6`
relative on the recoverable grid at 224–2237; 10M-edge grid, 20
iterations, completes on Windows/NVIDIA and Apple silicon within 1.5× of
the §1 target; device memory within the budget in `AMG_PLAN.md`; no
allocation inside a solve.

#### WS-J — Crossover search and cost models

Scope: execute §5.3 Stages A–D on the §5.4 machines with WS-E tooling;
produce `crossover-<machine>.md`/`.json`; recommend default
`IterativeSolverOptions`; write the consolidated table and a proposal for a
future `Auto` rule (not implemented).
Owned: `benchmarks/**`, `BENCHMARKS.md` crossover section, default values
in `IterativeSolverOptions::default()` (via the integrator).
Acceptance: Stage C validation passes on every machine; fits have
residuals < 15%; defaults committed with the evidence linked.

#### WS-K — Robustness, documentation, release readiness

Scope: adversarial fixtures and failure-mode tests of §4.6; scale smoke
tests; CI hardware job wiring; `NATIVE_ARTIFACTS.md`, README and in-plugin
help text; a release checklist (toolchain, wgpu version pin, tested driver
versions, known limits: indefinite `q`, extreme `q` ratios, adapter
limits).
Owned: `tests/robustness_*.rs`, `.github/workflows/ci.yml` hardware job,
docs.
Acceptance: every failure mode in §4.6 has a test; the release checklist is
complete; program acceptance in §8 verified and recorded.

### 7.4 File ownership matrix (conflict avoidance)

| path | owner | others may |
|---|---|---|
| `src/linear_solver/mod.rs`, `direct.rs`, `src/backend/mod.rs`, `TheseusError` | WS-I | propose changes via integrator |
| `src/linear_solver/tolerance.rs`, dispatch sites in `fdm.rs`/`gradients.rs`/`optimizer.rs`, `SolverOptions`/`SolverResult` fields | WS-D | — |
| `src/graph/**`, `src/backend/cpu.rs`, loop bodies in `fdm.rs`/`gradients.rs`/`objectives.rs` | WS-B | WS-C read |
| `src/amg/**` | WS-C (WS-H for GPU specialisations) | — |
| `src/backend/gpu/**` | WS-G then WS-H | — |
| `src/ffi.rs` new functions, all C#, component files, `Ariadne.Tests` | WS-F | WS-H adds GPU probe details |
| `tests/support/fixtures/**`, `bench_scale.rs`, `profile_phases.rs`, `scripts/**`, `benchmarks/**` | WS-E | WS-A/WS-J add results |
| `examples/matrix_free_pcg.rs`, `examples/amg_prototype.rs` | WS-A | — |
| `.github/workflows/ci.yml` | WS-K (software-adapter jobs may be added by WS-G) | — |
| `BENCHMARKS.md` | shared, append-only sections per workstream | — |
| this document, `AMG_PLAN.md` | WS-I | PR comments |

### 7.5 Dispatch prompt template

Every agent receives: (1) this document and `AMG_PLAN.md`; (2) its
workstream section; (3) the branch name to create and the base branch; (4)
the instruction to run `cargo test --workspace --release`, clippy and fmt
on its owned files before each push, to record benchmark JSON under
`benchmarks/results/`, and to report back with: what was built, test
evidence, numbers, deviations from this document, and open questions for the
integrator. Agents do not edit files outside their ownership row without
noting it in the report.

---

## 8. Milestones and definition of done

| milestone | done when |
|---|---|
| M0 Interfaces | **reached** (`9a5f52d`): WS-I merged; `DirectSolver` bitwise-equal to the cache path, `factor_and_solve` shares its refactor code; full `Box<dyn LinearSystemSolver>` dispatch inside `FdmCache` deferred to WS-D |
| M1 Phase-0 verdict | **reached** (`e3d0cf6`): plain aggregation no-go; smoothed aggregation + V-cycle PCG adopted (§3); wall-time criterion not met single-threaded (2.5–3.6× direct per evaluation at 1M) — CPU viability hinges on the frozen-pattern update and 4-thread kernels, GPU on WS-H |
| M2 Infrastructure | **reached** — WS-B, WS-E, WS-F merged (*WS-B done*, `16da146`; *WS-F done*, `cb562e6`: toggle, options and probe through FFI/C#/Grasshopper, iterative kinds return code `-4`); toggle visible in Grasshopper returning "not yet available" for iterative kinds; sweep tooling produces a `Direct` cost model. *WS-E done* (`6df10ba`): fixtures, JSON harness, `scripts/bench_sweep.py`, `crossover_fit.py`, first `Direct` sweep 10k–1M on the CI VM with fits at RMS < 15% (`benchmarks/reports/`) |
| M3 CPU iterative | WS-C, WS-D merged; `IterativeCpu` passes all contract tests; `bench_scale` numbers at 224–708 recorded. *WS-D done* (`d52b2b1`): every FDM solve dispatches through `Box<dyn LinearSystemSolver>`, `Direct` byte-identical (golden dump at 1 and 4 threads), tolerance schedule, totals, typed-error propagation; `IterativeCpu` legs of the contract tests skip until WS-C lands |
| M4 GPU kernels | WS-G merged (`1c51afc`); kernels validated under lavapipe (17/17, f32 + f64); Windows/macOS software-adapter runs pending the CI job, real-GPU validation pending WS-H machines |
| M5 GPU end to end | WS-H merged; `IterativeGpu` validated on Windows/NVIDIA and Apple silicon; 10M-edge run recorded |
| M6 Crossover | WS-J reports for ≥ 3 machines; defaults committed |
| M7 Release ready | WS-K complete; §4 all green; docs and checklist done |

Program acceptance = M7 plus: `Direct` remains the default everywhere;
selecting `IterativeGpu` without an adapter fails loudly; a 10M-edge grid
solves with `IterativeGpu` on both target platforms; crossover data exists
for a future `Auto` decision.

---

## 9. Risk register

| risk | signal | mitigation / kill criterion |
|---|---|---|
| ~~Plain aggregation not grid-independent~~ **materialised** (WS-A: 3.5–4× growth) | — | switched to smoothed aggregation (§3); new risk below |
| SA hierarchy update too slow per `q` change | `update(q)` > 100 ms single-threaded at 1M after WS-C | frozen-pattern numeric triple product; coarse-level updates every k-th evaluation with the fine level exact (preconditioner staleness is harmless for PCG correctness, only iteration counts move); kill criterion: if `IterativeCpu` per-evaluation time on 4 threads exceeds `Direct` at 1M on all fixtures after WS-C+WS-D, the CPU backend is kept for correctness/reference only and the program's performance case rests on `IterativeGpu` |
| Inexact gradients stall L-BFGS-B line search | evaluations per iteration > 2 with `IterativeCpu` vs `Direct` on the same fixture | tighten `Adaptive.factor`; `Fixed(1e-10)` fallback; report in WS-D |
| `f32` preconditioner degrades convergence | GPU iterations > 1.5× CPU `f64` iterations | `f64` preconditioner where `SHADER_F64`; mixed-precision smoothing only on fine levels |
| Adapter limits (`max_storage_buffer_binding_size` < vector size) | probe on Intel iGPU / Metal | split bindings; document maximum size per adapter class |
| Readback latency dominates at small sizes | GPU slower than CPU below ~300k edges | expected; the toggle is explicit and the crossover report says where GPU pays |
| Device loss / driver bugs | validation matrix failures | error out with adapter info; document tested driver versions; no retry loops |
| Memory at 10M on 16 GB machines | peak RSS > 4 GB | audit against budget; `f32` edge scratch where accuracy permits; stream fixture construction |
| wgpu version churn | build breaks on update | pin exact version; upgrade only in WS-K with the validation matrix rerun |
| Parallel workstreams conflict in `fdm.rs`/`gradients.rs` | merge conflicts | ownership matrix; WS-B lands before WS-D touches dispatch sites; integrator rebases |
| Indefinite systems requested with iterative solver | soft bounds allowing `q < 0` | typed error at solve start, message names the toggle |

---

## 9a. Status at pause (2026-09-18) and how to resume

Work was paused with the program branch `cursor/theseus-scale-100k-2d27`
(PR #14) at the WS-C merge. Read this section first when picking up.

### What is on the program branch

| workstream | state | branch (merged) |
|---|---|---|
| WS-I interfaces, `DirectSolver` | merged | `cursor/its-i-interfaces-2d27` |
| WS-A Phase-0 prototype + verdict | merged | `cursor/its-a-phase0-2d27` |
| WS-B `graph::build`, `CpuBackend`, parallel loops | merged | `cursor/its-b-graph-loops-2d27` |
| WS-E fixtures, JSON harness, `scripts/`, first `Direct` sweep | merged | `cursor/its-e-bench-tooling-2d27` |
| WS-F FFI/C#/Grasshopper toggle (+ `max_device_bytes`) | merged | `cursor/its-f-toggle-2d27` |
| WS-G wgpu backend behind feature `gpu`, WGSL kernels, CI job | merged | `cursor/its-g-gpu-kernels-2d27` |
| WS-D `FdmCache` dispatch, tolerance schedule, diagnostics | merged | `cursor/its-d-fdm-dispatch-2d27` |
| WS-C smoothed-aggregation AMG (`IterativeCpu`), CSR kernels | merged (`f79ec1d`) — **3 test binaries fail**, see below | `cursor/its-c-amg-core-2d27` |
| WS-C/WS-D integration fixes | **in flight** at pause: branch `cursor/its-d2-amg-integration-2d27` (may or may not have been pushed; check `git ls-remote origin 'cursor/its-d2-*'`) | — |
| WS-H GPU end-to-end, WS-J crossover, WS-K robustness/docs | not started; specifications in §7.3 | — |

Milestones: M0, M1, M2 reached; M4 reached except real-hardware
validation; M3 blocked only on the three integration fixes below.

### Known failing tests at `f79ec1d` (all `IterativeCpu` legs; `Direct` is fully green)

Run `cargo test --workspace --release --no-fail-fast` to see all three at
once (`cargo test` stops at the first failing binary otherwise).

1. `optimization_diagnostic` (4 tests): soft-bounds problems drive `q < 0`;
   the AMG coarsest-level faer LLT fails and surfaces as
   `TheseusError::Linalg("CholeskyError …")`. Required: `AmgSolver::update`
   checks `q > 0` and finite up front and returns
   `IterativeSolverUnsupported` (message names the toggle); a coarsest
   Cholesky failure maps to the same typed error.
2. `linear_solver_dispatch::iterative_max_iterations_one_does_not_converge`:
   `max_iterations = 1` returned `Ok`. Check that `AmgSolver::solve` honours
   `SolveRequest.max_iterations`, reports `converged = false` on budget
   exhaustion, and that `run_recorded_solve` turns that into
   `IterativeSolverDidNotConverge`; or the tiny test problem converges in
   one iteration and the test needs a harder system.
3. `linear_solver_dispatch::every_solve_is_recorded_in_the_cache_totals`:
   `Solver("Neumann adjoint did not converge within 30 iterations
   (tolerance 1e-10)")` — `solve_modified_adjoint`'s inner iterative solves
   run at `cache.solve_tolerance` (1e-6 ceiling); they must run ≥ 100×
   tighter than the Neumann target (clamped ≥ 1e-14). `Direct` ignores
   tolerances and must stay byte-identical.

Also to verify at the same time (WS-D → WS-C hand-over, §10): `AmgSolver`
overrides `precondition` with one V-cycle; `x0` is copied, not used in
place; no 1e-12 diagonal shift in the AMG operator (matches the
matrix-free `apply_a_xyz`).

### Resume checklist

```
git fetch origin
git checkout cursor/theseus-scale-100k-2d27
git ls-remote origin 'cursor/its-d2-*'          # integration fixes pushed? merge them first
export RUSTUP_TOOLCHAIN=1.89.0
cargo test --workspace --release --no-fail-fast   # expect green after the fixes
# GPU-feature build and tests (software adapter):
sudo apt-get install -y mesa-vulkan-drivers libvulkan1
THESEUS_GPU_ALLOW_SOFTWARE=1 VK_ICD_FILENAMES=/usr/share/vulkan/icd.d/lvp_icd.json \
  cargo test -p theseus --release --features gpu
# CPU numbers (ignored benches):
RAYON_NUM_THREADS=4 cargo test -p theseus --release --test amg_bench -- --ignored --nocapture
THESEUS_LINEAR_SOLVER=iterative-cpu THESEUS_FIXTURE=grid cargo test -p theseus --release --test bench_scale -- --ignored --nocapture
```

Then dispatch, in this order (prompts derive from §7.3 + §10 decisions):

1. **WS-H** — `AmgSolver<GpuBackend>`: f32 preconditioner with f64 host
   outer loop (device outer loop where `SHADER_F64`), fallible
   `Backend::alloc`/`upload_level` (decision in §10), binding splitting for
   adapters with `max_storage_buffer_binding_size` < largest array (needed
   for 10M edges on lavapipe/iGPU), per-dispatch bind-group/param-upload
   removal, `LinearSolver::new(IterativeGpu, …)`, GPU probe cached per
   handle, `SolveStats.backend`/adapter reporting through FFI/C#. Validate
   under lavapipe here, then on Windows/NVIDIA and Apple silicon (§5.4).
2. **WS-J** — Stage A–C of §5.3 on this VM for `Direct` vs `IterativeCpu`
   (tooling: `uv run --project scripts scripts/bench_sweep.py`, then
   `crossover_fit.py`; see `benchmarks/README.md`); Stage D on GPU machines
   after WS-H. WS-C's measured AMG defaults (`amg::recommended_options()`:
   3 passes, Chebyshev 2, α = 10, V-cycle) are the Stage-A starting point;
   `IterativeSolverOptions::default()` still holds the §2.1 values (K-cycle,
   2 passes) — WS-J decides which becomes the default and the integrator
   updates §2.1 and the C# defaults together.
3. **WS-K** — §4.6 failure-mode tests, the every-objective fixture,
   `gpu-software` CI job hardening (remove `continue-on-error`), hardware
   CI job, docs, release checklist, crate-wide clippy sweep (≈70
   pre-existing findings make `-D warnings` unusable as a gate today).

### Measured standing at pause (4-vCPU Linux VM, from `BENCHMARKS.md` WS-C section)

`IterativeCpu` update+solve vs `Direct` refactor+solve at 1M edges: 4
threads 0.36–0.54× cold, 0.27–0.43× warm; 1 thread 1.10–1.67× cold,
0.77–1.31× warm; `update(q)` 35–54 ms (1 thread) / 20–24 ms (4 threads),
targets met. Break-even on 4 threads is near 100k edges; single-threaded
the direct solver still wins below ~1M edges.

## 10. Decision log

* 2026-09-18 — Direct sparse Cholesky retained as default; iterative
  backends explicit only (owner request).
* 2026-09-18 — wgpu chosen as the portable GPU layer; mixed precision
  required because Metal and WGSL lack `f64`; CUDA excluded from scope.
* 2026-09-18 — Aggregation AMG with K-cycle/FCG as the starting algorithm,
  gated by the WS-A prototype; smoothed aggregation as the documented
  fallback.
* 2026-09-18 — Crossover search produces data for a future `Auto` mode but
  does not change selection behaviour in this program.
* 2026-09-18 — WS-I landed; interface refinements adopted: `LinearSystemSolver:
  Send`, `precondition_precision: Option<Precision>`, `Backend` associated
  buffer types and fused `residual` kernel, unit-struct `LinearSolver`
  factory, perturbation as solver state, FFI codes `-3..-6`.
* 2026-09-18 — WS-E landed. Harness backend strings are fixed as
  `direct`, `iterative-cpu`, `iterative-gpu` (`THESEUS_LINEAR_SOLVER`; WS-D
  implements the parser). Under external load the sweep report keeps the
  fastest unflagged re-run per cell and records every pass in
  `config.json`; clean-machine sweeps (WS-J) use plain medians. Crossover
  JSON is keyed by the machine's max-thread setting; the single-thread
  crossover is reported alongside for the `Auto` proposal. Not yet built:
  the "every objective type" 1,000-node fixture (moves to WS-K) and
  in-harness full-solve repetitions (the sweep re-runs cells instead).
* 2026-09-18 — WS-F landed (`cb562e6`): `theseus_set_linear_solver`,
  `theseus_set_iterative_options` (flat parameters, `-1` = backend-default
  precision), `theseus_gpu_probe` (JSON, returns bytes required),
  `theseus_get_linear_solver_stats` (`TheseusLinearSolverStats`, 56 bytes);
  pre-dispatch guard `ffi::begin_linear_solver_run` returns code `-4` for
  non-`Direct` kinds until WS-D's `FdmCache` dispatch replaces it;
  `LinearSolverTotals::record` is called by WS-D from the dispatch site.
  Decisions: `max_device_bytes: u64` (0 = adapter default) is appended to
  `theseus_set_iterative_options` before the C# signature freezes (WS-F
  follow-up); the GPU pre-probe is cached per handle once WS-G lands and
  re-run only on explicit request (WS-H); preserving wires when the toggle
  rebuilds optional inputs goes to WS-K.
* 2026-09-18 — WS-B landed (`16da146`): `graph::build` (`CsrAdjacency::
  from_topology`, `Level0Map::update_weights(q)`, `LevelGraph::level0`),
  `backend::cpu::CpuBackend` with the chunked deterministic helpers
  (`CHUNK = 4096`, `PAR_MIN_LEN`, `deterministic_sum`), `FdmCache.adjacency`,
  `boundary_edges` sorted by (free row, edge). Non-solver phases 2.4–2.7×
  faster on 4 threads at 1M edges, no single-thread regression, bitwise
  identical across thread counts (tests `graph_loops_determinism`). Only the
  objective loss reductions reassociate (1e-12 tests). Decisions:
  `run_sequential`'s one-thread shortcut stays (results identical);
  `Level0Map` lives in the iterative solver state, not `FdmCache` (WS-C);
  the sequential triangular solves (~60% of non-factor time) are a WS-K
  candidate (parallel over the 3 columns or supernode-level parallelism).
* 2026-09-18 — WS-G landed (`4b82feb`, merged `1c51afc`): cargo feature
  `gpu` (`wgpu =30.0.1`, `pollster`, `bytemuck`; `ordered-float` pinned to
  5.4 for the 1.89 toolchain), `backend::gpu` with adapter policy, buffer
  pool, batched encoder, f32 kernels and f64 variants where `SHADER_F64`
  (Vulkan only in wgpu-hal 30; naga accepts `f64` without an `enable`),
  17/17 equivalence tests under lavapipe, `gpu-software` CI job
  (`continue-on-error` until WS-K). Probe types unified at merge: one
  `GpuProbe`/`GpuAdapterReport`, JSON keeps WS-F's keys and adds `driver`,
  `max_buffer_size`, `software`, `selectable`, `built_with_gpu_feature`.
  Decisions on WS-G's questions: (1) `Backend::alloc`/`upload_level`
  become fallible (`Result<_, TheseusError>`) — WS-H changes the trait and
  the CPU impl together, the integrator updates §2.3; (2) binding splitting
  for `max_storage_buffer_binding_size` < largest array is WS-H scope,
  required for 10M edges on lavapipe/iGPU-class adapters; (3) coarse-edge
  ordering convention is ascending `(min(U,V), max(U,V))` — WS-C's
  `hierarchy.rs` must emit exactly that; (4) per-dispatch bind-group and
  `write_buffer` costs are removed in WS-H once a K-cycle is measured;
  (5) Windows/macOS adapter limit rows are filled from the CI job output.
* 2026-09-18 — WS-A landed (`e3d0cf6`, `examples/amg_prototype.rs`, Phase-0
  section in `BENCHMARKS.md`). Verdict: plain aggregation no-go; smoothed
  aggregation (`P = (I − ωD⁻¹A)P₀`, 3 passes, Chebyshev 2, α = 10) with a
  V-cycle inside PCG is size-independent within +30% on irregular/dome
  and +25% (x/y) / +47% (z) on the grid; K-cycle and FCG dropped. §3
  rewritten accordingly; `AMG_PLAN.md` budgets ("≈10 operator-equivalents
  per iteration", "O(nnz) coarse update", "parity single-threaded at 1M")
  are superseded by the measured numbers. Consequences for WS-C: general
  sparse `LevelMatrix` for levels ≥ 1, frozen-pattern numeric update as
  the primary optimisation target; for WS-H: a CSR SpMV kernel replaces
  `coarse_weight_update` for levels ≥ 1.
* 2026-09-18 — WS-D landed (`d52b2b1`). `FdmCache` owns only
  `linear_solver: Box<dyn LinearSystemSolver>` (+ kind, totals, per-solve
  tolerance, warm-start vectors); `a_matrix`/`factorization`/`q_to_nz`
  moved into `DirectSolver` and are borrowed via `as_direct()`. Iterative
  kinds use a matrix-free `apply_a_xyz` over `FdmCache.adjacency` without
  the direct path's 1e-12 diagonal shift (the AMG operator must match).
  `ToleranceSchedule` lives in `optimizer.rs::run_solver` (observe after
  init and after every accepted step) and is written to
  `cache.solve_tolerance` before each evaluation. `run_recorded_solve`
  records totals, checks finiteness (iterative only) and turns
  `!converged` into `IterativeSolverDidNotConverge`; preconditioner
  applications are recorded but never fail. The `q ≤ 0` bounds check in
  `FdmCache::new` is authoritative; WS-C's is a backstop. `SolverResult`
  gains `linear_solver_iterations` (per evaluation) and
  `linear_solver_totals`; the termination suffix `linear solver: <kind>,
  N solves, M iterations` appears only for non-`Direct` kinds. Handed to
  WS-C at merge: override `precondition` with one V-cycle; `x0` is copied
  into the solution before iterating (not used in place); skip the
  dispatch-side finiteness pass if the solver guarantees finite output.
  Handed to WS-H/WS-K: C# consumes per-evaluation iteration counts as-is.
