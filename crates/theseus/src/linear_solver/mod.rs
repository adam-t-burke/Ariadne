//! Linear solver selection and the [`LinearSystemSolver`] interface.
//!
//! Every forward FDM solve and every adjoint solve is a system
//! `A(q) x = b` with `A = Cnᵀ diag(q) Cn` and a block of three right-hand
//! sides (x, y, z). This module defines
//!
//! * [`LinearSolverKind`] — the explicit, user-facing toggle between the
//!   sparse direct factorization (default) and the iterative backends;
//! * [`IterativeSolverOptions`] and its parameter enums — the knobs of the
//!   iterative solvers, ignored when the kind is `Direct`;
//! * [`SolveRequest`] / [`SolveStats`] / [`MemoryReport`] — the per-solve
//!   contract shared by all implementations;
//! * the [`LinearSystemSolver`] trait and the [`LinearSolver`] factory.
//!
//! The direct adapter lives in [`direct`]; the CPU iterative solver is
//! [`crate::amg::AmgSolver`] over [`crate::backend::CpuBackend`]. The GPU
//! kind plugs into [`LinearSolver::new`] in a later workstream; until then
//! requesting it yields [`TheseusError::IterativeSolverUnsupported`].
//!
//! Enum discriminants that cross the FFI boundary (`LinearSolverKind`,
//! `TolerancePolicy` mode, `CycleKind`, `Precision`, `GpuOuterLoop`,
//! `AdapterPreference`) are stable `i32` values and must not be renumbered.

pub mod direct;

use crate::types::{Bounds, NetworkTopology, TheseusError};
use std::fmt;
use std::sync::atomic::AtomicBool;

pub use direct::DirectSolver;

// ─────────────────────────────────────────────────────────────
//  Solver kind (the toggle)
// ─────────────────────────────────────────────────────────────

/// Which linear solver evaluates `A(q) x = b` and the adjoint system.
///
/// Serialised as `i32` across FFI; values are stable. `Direct` is the
/// default at every layer and is never changed automatically.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum LinearSolverKind {
    /// Sparse Cholesky / LDLᵀ via faer (today's path).
    #[default]
    Direct = 0,
    /// Matrix-free flexible CG with an aggregation-AMG preconditioner on the
    /// CPU backend.
    IterativeCpu = 1,
    /// Same algorithm with the preconditioner (and, where the adapter allows,
    /// the outer iteration) on a wgpu device.
    IterativeGpu = 2,
}

impl LinearSolverKind {
    /// Stable FFI discriminant.
    pub const fn as_i32(self) -> i32 {
        self as i32
    }

    /// `true` for the iterative kinds.
    pub const fn is_iterative(self) -> bool {
        !matches!(self, Self::Direct)
    }

    /// The compute backend an iterative kind runs on; `None` for `Direct`.
    pub const fn backend(self) -> Option<crate::backend::BackendHandle> {
        match self {
            Self::Direct => None,
            Self::IterativeCpu => Some(crate::backend::BackendHandle::Cpu),
            Self::IterativeGpu => Some(crate::backend::BackendHandle::Gpu),
        }
    }
}

impl TryFrom<i32> for LinearSolverKind {
    type Error = TheseusError;

    fn try_from(value: i32) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(Self::Direct),
            1 => Ok(Self::IterativeCpu),
            2 => Ok(Self::IterativeGpu),
            _ => Err(TheseusError::Shape(format!(
                "invalid linear solver kind: {value} (expected 0 = Direct, 1 = IterativeCpu, 2 = IterativeGpu)"
            ))),
        }
    }
}

impl From<LinearSolverKind> for i32 {
    fn from(kind: LinearSolverKind) -> Self {
        kind.as_i32()
    }
}

impl fmt::Display for LinearSolverKind {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Self::Direct => "Direct",
            Self::IterativeCpu => "IterativeCpu",
            Self::IterativeGpu => "IterativeGpu",
        })
    }
}

// ─────────────────────────────────────────────────────────────
//  Iterative solver options
// ─────────────────────────────────────────────────────────────

/// How the relative-residual tolerance of each iterative solve is chosen.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TolerancePolicy {
    /// The same relative residual for every solve.
    Fixed(f64),
    /// `tol_k = clamp(factor · ‖g⁺_k‖ / ‖g⁺_0‖, floor, ceiling)` from the
    /// optimizer's projected-gradient ratio; the first evaluation uses
    /// `ceiling`.
    Adaptive {
        floor: f64,
        ceiling: f64,
        factor: f64,
    },
}

impl TolerancePolicy {
    /// FFI mode discriminant of [`TolerancePolicy::Fixed`].
    pub const FIXED: i32 = 0;
    /// FFI mode discriminant of [`TolerancePolicy::Adaptive`].
    pub const ADAPTIVE: i32 = 1;

    /// Stable FFI mode discriminant.
    pub const fn mode(&self) -> i32 {
        match self {
            Self::Fixed(_) => Self::FIXED,
            Self::Adaptive { .. } => Self::ADAPTIVE,
        }
    }

    /// The tolerance used before any gradient information exists.
    pub fn initial(&self) -> f64 {
        match self {
            Self::Fixed(tol) => *tol,
            Self::Adaptive { ceiling, .. } => *ceiling,
        }
    }
}

impl Default for TolerancePolicy {
    fn default() -> Self {
        Self::Adaptive {
            floor: 1e-10,
            ceiling: 1e-6,
            factor: 1e-2,
        }
    }
}

/// Multigrid cycle used as the preconditioner.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum CycleKind {
    /// One recursive pass (one pre- and one post-smoothing sweep): a fixed
    /// SPD preconditioner for standard PCG. The configuration recommended
    /// by the Phase-0 revision of §3 (smoothed aggregation keeps its
    /// iteration counts size-independent; see `BENCHMARKS.md`, "WS-C").
    V = 0,
    /// Two flexible-CG iterations per coarse level, recursively (Notay);
    /// the outer loop then runs FCG(1). Still the `Default` because the
    /// FFI and C# layers mirror this value (§2.3); moving every layer to
    /// `V` is the integrator's / WS-J's call (§7.4).
    #[default]
    K = 1,
}

impl TryFrom<i32> for CycleKind {
    type Error = TheseusError;

    fn try_from(value: i32) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(Self::V),
            1 => Ok(Self::K),
            _ => Err(TheseusError::Shape(format!(
                "invalid multigrid cycle kind: {value} (expected 0 = V, 1 = K)"
            ))),
        }
    }
}

/// Floating-point precision of the preconditioner (the outer iteration and
/// the convergence residual are always `f64`).
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Precision {
    F64 = 0,
    F32 = 1,
}

impl Precision {
    /// Bytes per scalar.
    pub const fn size_of(self) -> usize {
        match self {
            Self::F64 => 8,
            Self::F32 => 4,
        }
    }
}

impl TryFrom<i32> for Precision {
    type Error = TheseusError;

    fn try_from(value: i32) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(Self::F64),
            1 => Ok(Self::F32),
            _ => Err(TheseusError::Shape(format!(
                "invalid precision: {value} (expected 0 = F64, 1 = F32)"
            ))),
        }
    }
}

/// Where the `f64` outer FCG iteration runs when the kind is `IterativeGpu`.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum GpuOuterLoop {
    /// On the device when it reports `SHADER_F64`, otherwise on the host.
    #[default]
    Auto = 0,
    /// Force the device; fails with `GpuUnavailable` without `SHADER_F64`.
    Device = 1,
    /// Force the host (the `f32` preconditioner stays on the device).
    Host = 2,
}

impl TryFrom<i32> for GpuOuterLoop {
    type Error = TheseusError;

    fn try_from(value: i32) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(Self::Auto),
            1 => Ok(Self::Device),
            2 => Ok(Self::Host),
            _ => Err(TheseusError::Shape(format!(
                "invalid GPU outer-loop placement: {value} (expected 0 = Auto, 1 = Device, 2 = Host)"
            ))),
        }
    }
}

/// Which adapter class to prefer when several are present. Software (CPU)
/// adapters are rejected unless `THESEUS_GPU_ALLOW_SOFTWARE=1`.
#[repr(i32)]
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Default)]
pub enum AdapterPreference {
    /// Discrete first, then integrated.
    #[default]
    Discrete = 0,
    /// Integrated first, then discrete.
    Integrated = 1,
    /// First usable adapter in enumeration order.
    Any = 2,
}

impl TryFrom<i32> for AdapterPreference {
    type Error = TheseusError;

    fn try_from(value: i32) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(Self::Discrete),
            1 => Ok(Self::Integrated),
            2 => Ok(Self::Any),
            _ => Err(TheseusError::Shape(format!(
                "invalid adapter preference: {value} (expected 0 = Discrete, 1 = Integrated, 2 = Any)"
            ))),
        }
    }
}

/// GPU-specific options; ignored unless the kind is `IterativeGpu`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct GpuOptions {
    /// Placement of the `f64` outer iteration.
    pub outer_loop: GpuOuterLoop,
    /// Adapter class to prefer.
    pub adapter_preference: AdapterPreference,
    /// Cap on device memory the solver may allocate; `None` = adapter limit.
    pub max_device_bytes: Option<u64>,
}

/// Parameters of the iterative linear solvers (§2.1 of the program plan).
///
/// Ignored when [`SolverOptions::linear_solver`](crate::types::SolverOptions)
/// is [`LinearSolverKind::Direct`].
#[derive(Debug, Clone, PartialEq)]
pub struct IterativeSolverOptions {
    /// Relative-residual tolerance schedule.
    pub tolerance: TolerancePolicy,
    /// Iteration budget per solve.
    pub max_iterations: u32,
    /// Multigrid cycle used as the preconditioner.
    pub cycle: CycleKind,
    /// Chebyshev smoother degree.
    pub smoother_degree: u8,
    /// Spectral ratio `α` of the Chebyshev smoother: the polynomial damps
    /// `[λ_max/α, λ_max]` of `D⁻¹A`. Must be `> 1`.
    pub spectral_alpha: f64,
    /// Pairwise matching passes per level (aggregates of up to
    /// `2^passes` nodes).
    pub aggregation_passes: u8,
    /// Stop coarsening once a level has fewer nodes than this.
    pub coarsest_size: u32,
    /// Precision of the preconditioner. `None` = backend default
    /// (`F64` on CPU, `F32` on GPU); see [`Self::precondition_precision_for`].
    pub precondition_precision: Option<Precision>,
    /// GPU adapter preference, memory cap and outer-loop placement.
    pub gpu: GpuOptions,
}

impl IterativeSolverOptions {
    /// Default `max_iterations`.
    pub const DEFAULT_MAX_ITERATIONS: u32 = 200;
    /// Default `smoother_degree`.
    pub const DEFAULT_SMOOTHER_DEGREE: u8 = 2;
    /// Default `spectral_alpha` (§3: α = 10; α = 30 costs +40% iterations).
    pub const DEFAULT_SPECTRAL_ALPHA: f64 = 10.0;
    /// Default `aggregation_passes`. §3 recommends three passes (aggregates
    /// of ≤ 8 nodes) with smoothed aggregation and the WS-C measurements
    /// use them; the value stays at the §2.3 interface default mirrored by
    /// the FFI and C# layers until the integrator moves every layer.
    pub const DEFAULT_AGGREGATION_PASSES: u8 = 2;
    /// Default `coarsest_size`.
    pub const DEFAULT_COARSEST_SIZE: u32 = 2000;

    /// The preconditioner precision to use with `kind`: the explicit setting
    /// if any, else `F64` on the CPU and `F32` on the GPU.
    pub fn precondition_precision_for(&self, kind: LinearSolverKind) -> Precision {
        self.precondition_precision.unwrap_or(match kind {
            LinearSolverKind::IterativeGpu => Precision::F32,
            LinearSolverKind::Direct | LinearSolverKind::IterativeCpu => Precision::F64,
        })
    }
}

impl Default for IterativeSolverOptions {
    fn default() -> Self {
        Self {
            tolerance: TolerancePolicy::default(),
            max_iterations: Self::DEFAULT_MAX_ITERATIONS,
            cycle: CycleKind::default(),
            smoother_degree: Self::DEFAULT_SMOOTHER_DEGREE,
            spectral_alpha: Self::DEFAULT_SPECTRAL_ALPHA,
            aggregation_passes: Self::DEFAULT_AGGREGATION_PASSES,
            coarsest_size: Self::DEFAULT_COARSEST_SIZE,
            precondition_precision: None,
            gpu: GpuOptions::default(),
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Per-solve contract
// ─────────────────────────────────────────────────────────────

/// One request to solve `A x = rhs` for a block of three right-hand sides.
///
/// `rhs`, `x0` and the output `x` are row-major `n × 3` (`v[node * 3 + k]`),
/// the layout of `FdmCache::rhs` / `FdmCache::x` / `FdmCache::grad_x`.
#[derive(Debug, Clone, Copy)]
pub struct SolveRequest<'a> {
    /// Right-hand sides, `n * 3`.
    pub rhs: &'a [f64],
    /// Warm start, `n * 3`. `None` starts from zero. Ignored by `Direct`.
    pub x0: Option<&'a [f64]>,
    /// Relative residual `‖r‖ ≤ tolerance · ‖b‖` per column. Ignored by `Direct`.
    pub tolerance: f64,
    /// Iteration budget for this solve. Ignored by `Direct`.
    pub max_iterations: u32,
    /// Cooperative cancellation, checked between iterations.
    pub cancel: Option<&'a AtomicBool>,
}

impl<'a> SolveRequest<'a> {
    /// A request with `x0 = None`, no cancellation, and `tolerance` /
    /// `max_iterations` taken from `options` (initial tolerance of the policy).
    pub fn new(rhs: &'a [f64], options: &IterativeSolverOptions) -> Self {
        Self {
            rhs,
            x0: None,
            tolerance: options.tolerance.initial(),
            max_iterations: options.max_iterations,
            cancel: None,
        }
    }

    /// Number of nodes (`rhs.len() / 3`).
    pub fn num_nodes(&self) -> usize {
        self.rhs.len() / 3
    }
}

/// What one [`LinearSystemSolver::solve`] did.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SolveStats {
    /// Iterations per column (`0` for the direct solver).
    pub iterations: [u32; 3],
    /// Final relative residual per column (`0.0` for the direct solver, whose
    /// residual is not measured).
    pub relative_residual: [f64; 3],
    /// Every column reached its tolerance (always `true` for a successful
    /// direct solve).
    pub converged: bool,
    /// Time spent in `update` (factorization / hierarchy update) attributed
    /// to this solve, in milliseconds.
    pub setup_ms: f64,
    /// Wall time of `solve` itself, in milliseconds.
    pub solve_ms: f64,
    /// Solver kind that actually ran.
    pub backend: LinearSolverKind,
}

impl SolveStats {
    /// Largest per-column iteration count.
    pub fn max_iterations(&self) -> u32 {
        self.iterations.iter().copied().max().unwrap_or(0)
    }
}

/// Running totals over the [`SolveStats`] of one optimisation or forward
/// solve, as exported through the FFI (`theseus_get_linear_solver_stats`).
///
/// `converged_all` starts `true` and is cleared by the first recorded solve
/// that did not converge, so a run with no recorded solves reports
/// `converged_all == true` and `solves == 0`.
///
/// TODO(WS-D): nothing calls [`LinearSolverTotals::record`] yet. The FFI
/// resets the totals at the start of every run, but `factor_and_solve`
/// (`fdm.rs`) does not produce a [`SolveStats`], so `Direct` runs report
/// `solves == 0` and zero times. Once the dispatch inside `FdmCache` returns
/// `SolveStats`, `record` each of them into the handle's totals.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LinearSolverTotals {
    /// Solver kind the run was configured with (and, once the dispatch
    /// records stats, the kind that actually ran).
    pub backend: LinearSolverKind,
    /// Number of `solve` calls recorded.
    pub solves: u64,
    /// Sum over solves of the largest per-column iteration count.
    pub iterations_total: u64,
    /// Largest per-column iteration count of any single solve.
    pub iterations_max: u32,
    /// Sum of `solve_ms`.
    pub solve_ms_total: f64,
    /// Sum of `setup_ms`.
    pub setup_ms_total: f64,
    /// Every recorded solve converged.
    pub converged_all: bool,
}

impl LinearSolverTotals {
    /// Empty totals for a run on `backend`.
    pub fn new(backend: LinearSolverKind) -> Self {
        Self {
            backend,
            solves: 0,
            iterations_total: 0,
            iterations_max: 0,
            solve_ms_total: 0.0,
            setup_ms_total: 0.0,
            converged_all: true,
        }
    }

    /// Fold one solve into the totals.
    pub fn record(&mut self, stats: &SolveStats) {
        let iters = stats.max_iterations();
        self.backend = stats.backend;
        self.solves += 1;
        self.iterations_total += u64::from(iters);
        self.iterations_max = self.iterations_max.max(iters);
        self.solve_ms_total += stats.solve_ms;
        self.setup_ms_total += stats.setup_ms;
        self.converged_all &= stats.converged;
    }
}

impl Default for LinearSolverTotals {
    fn default() -> Self {
        Self::new(LinearSolverKind::Direct)
    }
}

/// Memory held by a solver, split by where it lives.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct MemoryReport {
    pub host_bytes: u64,
    pub device_bytes: u64,
}

impl MemoryReport {
    pub fn total_bytes(&self) -> u64 {
        self.host_bytes + self.device_bytes
    }
}

// ─────────────────────────────────────────────────────────────
//  The solver interface
// ─────────────────────────────────────────────────────────────

/// A solver for `A(q) x = b`, `A = Cnᵀ diag(q) Cn`, on a fixed topology.
///
/// Implementations are constructed for one topology (see [`LinearSolver::new`])
/// and then driven by the FDM code as `update(q)` once per new `q`, followed
/// by any number of `solve` calls against that `q` (forward and adjoint
/// systems share `A` because it is symmetric).
///
/// The trait requires `Send` so that a boxed solver can live inside
/// `FdmCache` when the cache is moved to a worker thread.
pub trait LinearSystemSolver: Send {
    /// Called whenever `q` changed. Direct: reassemble and refactor.
    /// Iterative: recompute level weights, spectral bounds if needed, upload
    /// to the device.
    fn update(&mut self, q: &[f64]) -> Result<(), TheseusError>;

    /// Solve `A x = req.rhs` for the three columns, writing `x` (`n * 3`).
    fn solve(&mut self, req: SolveRequest<'_>, x: &mut [f64]) -> Result<SolveStats, TheseusError>;

    /// The kind this solver implements.
    fn kind(&self) -> LinearSolverKind;

    /// Host and device memory currently held.
    fn memory_bytes(&self) -> MemoryReport;
}

// ─────────────────────────────────────────────────────────────
//  Factory
// ─────────────────────────────────────────────────────────────

/// Factory for [`LinearSystemSolver`] implementations, dispatching on
/// [`LinearSolverKind`].
///
/// This is the single place where a kind is turned into a solver; later
/// workstreams register the iterative implementations here.
#[derive(Debug, Clone, Copy)]
pub struct LinearSolver;

impl LinearSolver {
    /// Build the solver for `kind` on `topology`.
    ///
    /// `bounds` decides the direct factorization strategy (Cholesky when every
    /// lower bound is positive, LDLᵀ otherwise) and, for the iterative kinds,
    /// whether `q > 0` is guaranteed. `options` is ignored by `Direct`.
    // The factory deliberately returns the trait object, not `Self`: callers
    // hold a `Box<dyn LinearSystemSolver>` and never name the concrete type.
    #[allow(clippy::new_ret_no_self)]
    pub fn new(
        kind: LinearSolverKind,
        topology: &NetworkTopology,
        bounds: &Bounds,
        options: &IterativeSolverOptions,
    ) -> Result<Box<dyn LinearSystemSolver>, TheseusError> {
        match kind {
            LinearSolverKind::Direct => Ok(Box::new(DirectSolver::from_bounds(topology, bounds)?)),
            LinearSolverKind::IterativeCpu => Ok(Box::new(crate::amg::AmgSolver::cpu(
                topology, bounds, options,
            )?)),
            LinearSolverKind::IterativeGpu => {
                Err(TheseusError::IterativeSolverUnsupported(format!(
                    "the GPU iterative solver is not yet available (requested '{kind}'); \
                 use IterativeCpu or Direct"
                )))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kind_round_trips_through_i32() {
        for kind in [
            LinearSolverKind::Direct,
            LinearSolverKind::IterativeCpu,
            LinearSolverKind::IterativeGpu,
        ] {
            let code: i32 = kind.into();
            assert_eq!(LinearSolverKind::try_from(code).unwrap(), kind);
            assert_eq!(kind.as_i32(), code);
        }
        assert_eq!(LinearSolverKind::Direct.as_i32(), 0);
        assert_eq!(LinearSolverKind::IterativeCpu.as_i32(), 1);
        assert_eq!(LinearSolverKind::IterativeGpu.as_i32(), 2);
        assert!(LinearSolverKind::try_from(3).is_err());
        assert!(LinearSolverKind::try_from(-1).is_err());
    }

    #[test]
    fn direct_is_the_default_everywhere() {
        assert_eq!(LinearSolverKind::default(), LinearSolverKind::Direct);
        let options = crate::types::SolverOptions::default();
        assert_eq!(options.linear_solver, LinearSolverKind::Direct);
        assert_eq!(options.iterative, IterativeSolverOptions::default());
        assert!(!LinearSolverKind::Direct.is_iterative());
        assert!(LinearSolverKind::IterativeCpu.is_iterative());
        assert_eq!(LinearSolverKind::Direct.backend(), None);
    }

    #[test]
    fn iterative_option_defaults_match_the_program_plan() {
        let o = IterativeSolverOptions::default();
        assert_eq!(
            o.tolerance,
            TolerancePolicy::Adaptive {
                floor: 1e-10,
                ceiling: 1e-6,
                factor: 1e-2
            }
        );
        assert_eq!(o.tolerance.initial(), 1e-6);
        assert_eq!(o.max_iterations, 200);
        assert_eq!(o.cycle, CycleKind::K);
        assert_eq!(o.smoother_degree, 2);
        assert_eq!(o.spectral_alpha, 10.0);
        assert_eq!(o.aggregation_passes, 2);
        assert_eq!(o.coarsest_size, 2000);
        assert_eq!(
            o.precondition_precision_for(LinearSolverKind::IterativeCpu),
            Precision::F64
        );
        assert_eq!(
            o.precondition_precision_for(LinearSolverKind::IterativeGpu),
            Precision::F32
        );
        assert_eq!(o.gpu.outer_loop, GpuOuterLoop::Auto);
        assert_eq!(o.gpu.adapter_preference, AdapterPreference::Discrete);
        assert_eq!(o.gpu.max_device_bytes, None);
    }

    #[test]
    fn parameter_enums_round_trip_through_i32() {
        for v in [CycleKind::V, CycleKind::K] {
            assert_eq!(CycleKind::try_from(v as i32).unwrap(), v);
        }
        for v in [Precision::F64, Precision::F32] {
            assert_eq!(Precision::try_from(v as i32).unwrap(), v);
        }
        for v in [GpuOuterLoop::Auto, GpuOuterLoop::Device, GpuOuterLoop::Host] {
            assert_eq!(GpuOuterLoop::try_from(v as i32).unwrap(), v);
        }
        for v in [
            AdapterPreference::Discrete,
            AdapterPreference::Integrated,
            AdapterPreference::Any,
        ] {
            assert_eq!(AdapterPreference::try_from(v as i32).unwrap(), v);
        }
        assert!(CycleKind::try_from(2).is_err());
        assert!(Precision::try_from(2).is_err());
        assert!(GpuOuterLoop::try_from(3).is_err());
        assert!(AdapterPreference::try_from(3).is_err());
        assert_eq!(TolerancePolicy::Fixed(1e-8).mode(), TolerancePolicy::FIXED);
        assert_eq!(TolerancePolicy::default().mode(), TolerancePolicy::ADAPTIVE);
    }

    #[test]
    fn totals_accumulate_solve_stats() {
        let mut totals = LinearSolverTotals::new(LinearSolverKind::IterativeCpu);
        assert_eq!(
            totals,
            LinearSolverTotals {
                backend: LinearSolverKind::IterativeCpu,
                solves: 0,
                iterations_total: 0,
                iterations_max: 0,
                solve_ms_total: 0.0,
                setup_ms_total: 0.0,
                converged_all: true,
            }
        );
        totals.record(&SolveStats {
            iterations: [3, 7, 5],
            relative_residual: [1e-9; 3],
            converged: true,
            setup_ms: 1.5,
            solve_ms: 4.0,
            backend: LinearSolverKind::IterativeCpu,
        });
        totals.record(&SolveStats {
            iterations: [2, 2, 2],
            relative_residual: [1e-3; 3],
            converged: false,
            setup_ms: 0.5,
            solve_ms: 1.0,
            backend: LinearSolverKind::IterativeCpu,
        });
        assert_eq!(totals.solves, 2);
        assert_eq!(totals.iterations_total, 9);
        assert_eq!(totals.iterations_max, 7);
        assert_eq!(totals.solve_ms_total, 5.0);
        assert_eq!(totals.setup_ms_total, 2.0);
        assert!(!totals.converged_all);
        assert_eq!(
            LinearSolverTotals::default().backend,
            LinearSolverKind::Direct
        );
    }

    #[test]
    fn error_messages_name_the_toggle() {
        let e = TheseusError::IterativeSolverDidNotConverge {
            iterations: 200,
            relative_residual: 1e-3,
            kind: LinearSolverKind::IterativeCpu,
        }
        .to_string();
        assert!(e.contains("IterativeCpu") && e.contains("Direct"), "{e}");
        let e = TheseusError::IterativeSolverUnsupported("q may be negative".into()).to_string();
        assert!(
            e.contains("q may be negative") && e.contains("Direct"),
            "{e}"
        );
        let e = TheseusError::GpuUnavailable("no adapters".into()).to_string();
        assert!(e.contains("IterativeGpu") && e.contains("Direct"), "{e}");
        let e = TheseusError::GpuOutOfMemory {
            requested: 10,
            available: 5,
        }
        .to_string();
        assert!(
            e.contains("IterativeGpu") && e.contains("10") && e.contains("5"),
            "{e}"
        );
    }
}
