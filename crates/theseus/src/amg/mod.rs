//! Smoothed-aggregation algebraic multigrid (program plan §3) behind the
//! [`LinearSystemSolver`] interface: [`AmgSolver`].
//!
//! Module layout:
//!
//! * [`aggregate`] — pairwise heavy-edge matching → `aggregate_of` (`P₀`);
//! * [`hierarchy`] — [`LevelMatrix`](hierarchy::LevelMatrix) (CSR), the
//!   smoothed prolongator `P = (I − ω D⁻¹ A) P₀`, the symbolic pattern of
//!   `Pᵀ A P` and its numeric-only refill;
//! * [`smoother`] — Chebyshev with Jacobi scaling over a
//!   [`DeviceLevel`](smoother::DeviceLevel) (graph or CSR), `λ_max` power
//!   iteration;
//! * [`cycle`] — V-cycle (default) and K-cycle over the device hierarchy;
//! * [`coarsest`] — faer LLᵀ of the coarsest level;
//! * [`pcg`] — block-of-3 PCG with per-column stopping.
//!
//! # Life cycle
//!
//! `AmgSolver::new` builds the topology-only data (level-0 adjacency,
//! `Level0Map`, outer PCG buffers). The first `update(q)` runs **setup**:
//! aggregation, `λ_max`, `P`, the frozen coarse patterns, the device
//! uploads and the coarsest factorization. Later `update(q)` calls refill
//! level-0 weights and the coarse operators numerically on the frozen
//! patterns (`P` values frozen from the setup `q`) and refactor the
//! coarsest level; they run setup again only when `max_e |Δq_e /
//! q_e^{setup}| > 2` or when a solve needed more than twice the iterations
//! of the first post-setup solve two times in a row.
//!
//! After the first `solve`, `solve` allocates nothing on the heap.

pub mod aggregate;
pub mod coarsest;
pub mod cycle;
pub mod hierarchy;
pub mod pcg;
pub mod smoother;

use crate::backend::{Backend, BackendHandle, CpuBackend};
use crate::graph::build::free_node_map;
use crate::graph::{CsrAdjacency, Level0Map, LevelGraph};
use crate::linear_solver::{
    CycleKind, IterativeSolverOptions, LinearSolverKind, LinearSystemSolver, MemoryReport,
    Precision, SolveRequest, SolveStats,
};
use crate::types::{Bounds, NetworkTopology, TheseusError};
use aggregate::{aggregate, StrengthGraph};
use coarsest::CoarsestSolver;
use cycle::{Cycle, LevelWork, Transfer};
use hierarchy::{
    galerkin_numeric, galerkin_pattern, smoothed_prolongator, LevelMatrix, LevelRef, ScratchPool,
};
use pcg::{PcgBuffers, PcgOutcome};
use smoother::{
    estimate_lambda_max, max_relative_change, Chebyshev, DeviceLevel, LAMBDA_REESTIMATE_DRIFT,
};
use std::time::Instant;

/// Power iterations per level for `λ_max` at setup.
pub const POWER_ITERATIONS: u32 = 10;
/// Relative `q` drift since setup that forces a new setup (§3).
pub const RESETUP_DRIFT: f64 = 2.0;
/// A solve needing more than this multiple of the post-setup iteration
/// count is "slow"; two slow solves in a row force a new setup.
pub const RESETUP_ITERATION_FACTOR: u32 = 2;
/// Coarsening stops when a level would keep more than this fraction of the
/// nodes of the level above.
pub const MIN_COARSENING: f64 = 0.7;
/// Hard cap on the number of levels.
pub const MAX_LEVELS: usize = 60;
/// K-cycle inner tolerance (prototype default).
pub const KCYCLE_TOL: f64 = 0.25;

/// Counters for tests and benchmarks.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct AmgCounters {
    /// Full setups run (aggregation, `P`, patterns).
    pub setups: u64,
    /// `update` calls that refilled the frozen patterns only.
    pub numeric_updates: u64,
    /// `λ_max` re-estimations outside setup.
    pub lambda_reestimates: u64,
    /// `solve` calls.
    pub solves: u64,
}

/// Smoothed-aggregation AMG-preconditioned CG on backend `B`.
pub struct AmgSolver<B: Backend> {
    backend: B,
    options: IterativeSolverOptions,
    precision: Precision,
    smoother: Chebyshev,
    n: usize,
    num_edges: usize,
    level0: LevelGraph,
    map: Level0Map,
    perturbation: f64,

    // Host hierarchy (empty before the first update).
    matrices: Vec<LevelMatrix>,
    p: Vec<LevelMatrix>,
    pt: Vec<LevelMatrix>,
    /// `A_l P_l` on its frozen pattern (first stage of the Galerkin refill).
    ap: Vec<LevelMatrix>,
    aggregate_of: Vec<Vec<u32>>,
    lambda_max: Vec<f64>,
    /// Explicit level-0 matrix when the hierarchy has a single level.
    level0_matrix: Option<LevelMatrix>,
    coarsest: Option<CoarsestSolver>,
    scratch: ScratchPool,

    // Device hierarchy.
    levels: Vec<DeviceLevel<B>>,
    /// Level 0 in `f64` for the outer PCG when the preconditioner runs in
    /// another precision (`None` when `precision == F64`: `levels[0]` is
    /// used directly).
    outer_level0: Option<DeviceLevel<B>>,
    transfers: Vec<Transfer<B>>,
    work: Vec<LevelWork<B::Buf>>,
    outer: PcgBuffers<B::Buf>,
    x_dev: B::Buf,
    b_dev: B::Buf,
    /// Host scratch (`n * 3`): `λ_max` start vectors and, when the
    /// preconditioner precision differs from `f64`, the conversion staging.
    host_scratch: Vec<f64>,

    // Drift bookkeeping.
    q_setup: Vec<f64>,
    q_lambda: Vec<f64>,
    setup_iterations: Option<u32>,
    slow_solves: u32,
    resetup_requested: bool,
    last_update_ms: f64,
    counters: AmgCounters,
}

impl<B: Backend> std::fmt::Debug for AmgSolver<B> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("AmgSolver")
            .field("backend", &self.backend.handle())
            .field("n", &self.n)
            .field("num_edges", &self.num_edges)
            .field("precision", &self.precision)
            .field("levels", &self.level_sizes())
            .field("lambda_max", &self.lambda_max)
            .field("counters", &self.counters)
            .finish()
    }
}

/// Mutable pieces of one preconditioner application.
struct CycleScratch<'a, B: Backend> {
    precision: Precision,
    work: &'a mut [LevelWork<B::Buf>],
    outer: &'a mut PcgBuffers<B::Buf>,
    x_dev: &'a mut B::Buf,
    b_dev: &'a mut B::Buf,
    host: &'a mut [f64],
}

/// `z = M r` on `f64` device buffers, converting through `host` when the
/// preconditioner runs in another precision.
fn precondition_with<B: Backend>(
    cycle: &Cycle<'_, B>,
    precision: Precision,
    work: &mut [LevelWork<B::Buf>],
    host: &mut [f64],
    r: &B::Buf,
    z: &mut B::Buf,
) {
    let be = cycle.backend;
    if precision == Precision::F64 {
        be.copy(r, &mut work[0].b);
    } else {
        be.download(r, host);
        be.upload(host, &mut work[0].b);
    }
    cycle.apply(work);
    if precision == Precision::F64 {
        be.copy(&work[0].x, z);
    } else {
        be.download(&work[0].x, host);
        be.upload(host, z);
    }
}

/// Error for bounds that admit `q ≤ 0`.
fn check_bounds(bounds: &Bounds, kind: LinearSolverKind) -> Result<(), TheseusError> {
    if let Some((e, &lb)) = bounds
        .lower
        .iter()
        .enumerate()
        .find(|(_, &lb)| lb.is_nan() || lb <= 0.0)
    {
        return Err(TheseusError::IterativeSolverUnsupported(format!(
            "{kind} needs q > 0 on every edge (A(q) must be SPD), but the bounds allow \
             q ≤ 0 (lower bound {lb} on edge {e}); raise the lower bounds or select \
             LinearSolverKind::Direct through SolverOptions::linear_solver"
        )));
    }
    Ok(())
}

impl AmgSolver<CpuBackend> {
    /// CPU solver (`LinearSolverKind::IterativeCpu`).
    pub fn cpu(
        topology: &NetworkTopology,
        bounds: &Bounds,
        options: &IterativeSolverOptions,
    ) -> Result<Self, TheseusError> {
        Self::new(CpuBackend::new(), topology, bounds, options)
    }
}

impl<B: Backend> AmgSolver<B> {
    /// Topology-only construction on `backend`; the hierarchy is built by
    /// the first [`LinearSystemSolver::update`].
    ///
    /// Fails with [`TheseusError::IterativeSolverUnsupported`] when `bounds`
    /// admit `q ≤ 0` (the operator would not be SPD).
    pub fn new(
        backend: B,
        topology: &NetworkTopology,
        bounds: &Bounds,
        options: &IterativeSolverOptions,
    ) -> Result<Self, TheseusError> {
        let kind = match backend.handle() {
            BackendHandle::Cpu => LinearSolverKind::IterativeCpu,
            BackendHandle::Gpu => LinearSolverKind::IterativeGpu,
        };
        if bounds.lower.len() != topology.num_edges {
            return Err(TheseusError::Shape(format!(
                "{kind}: bounds have {} entries, topology has {} edges",
                bounds.lower.len(),
                topology.num_edges
            )));
        }
        check_bounds(bounds, kind)?;
        if options.smoother_degree == 0 {
            return Err(TheseusError::Shape(
                "iterative solver: smoother_degree must be ≥ 1".into(),
            ));
        }
        if options.spectral_alpha.is_nan() || options.spectral_alpha <= 1.0 {
            return Err(TheseusError::Shape(format!(
                "iterative solver: spectral_alpha must be > 1, got {}",
                options.spectral_alpha
            )));
        }
        let precision = options.precondition_precision_for(kind);
        let adjacency = CsrAdjacency::from_topology(topology);
        let node_to_free = free_node_map(topology);
        let (level0, map) = Level0Map::new(&adjacency, &node_to_free, &topology.free_node_indices);
        let n = level0.n;
        let outer = PcgBuffers::new(&backend, n, Precision::F64);
        let x_dev = backend.alloc(n * 3, Precision::F64);
        let b_dev = backend.alloc(n * 3, Precision::F64);
        Ok(Self {
            smoother: Chebyshev::new(u32::from(options.smoother_degree), options.spectral_alpha),
            backend,
            options: options.clone(),
            precision,
            n,
            num_edges: topology.num_edges,
            level0,
            map,
            perturbation: 0.0,
            matrices: Vec::new(),
            p: Vec::new(),
            pt: Vec::new(),
            ap: Vec::new(),
            aggregate_of: Vec::new(),
            lambda_max: Vec::new(),
            level0_matrix: None,
            coarsest: None,
            scratch: ScratchPool::new(),
            levels: Vec::new(),
            outer_level0: None,
            transfers: Vec::new(),
            work: Vec::new(),
            outer,
            x_dev,
            b_dev,
            host_scratch: vec![0.0; n * 3],
            q_setup: Vec::new(),
            q_lambda: Vec::new(),
            setup_iterations: None,
            slow_solves: 0,
            resetup_requested: false,
            last_update_ms: 0.0,
            counters: AmgCounters::default(),
        })
    }

    /// Diagonal (anchor) shift applied to `A` on the next `update`, like
    /// `DirectSolver::set_perturbation`.
    pub fn set_perturbation(&mut self, perturbation: f64) {
        self.perturbation = perturbation;
    }

    /// Current diagonal shift.
    pub fn perturbation(&self) -> f64 {
        self.perturbation
    }

    /// The backend.
    pub fn backend(&self) -> &B {
        &self.backend
    }

    /// Number of free nodes.
    pub fn num_nodes(&self) -> usize {
        self.n
    }

    /// Preconditioner precision.
    pub fn precision(&self) -> Precision {
        self.precision
    }

    /// Nodes per level (empty before the first `update`).
    pub fn level_sizes(&self) -> Vec<usize> {
        if self.levels.is_empty() {
            return Vec::new();
        }
        std::iter::once(self.n)
            .chain(self.matrices.iter().map(|m| m.n))
            .collect()
    }

    /// `λ_max(D⁻¹ A_l)` per level.
    pub fn lambda_max(&self) -> &[f64] {
        &self.lambda_max
    }

    /// Setup / update / solve counters.
    pub fn counters(&self) -> AmgCounters {
        self.counters
    }

    /// The level-0 graph with the weights of the last `update`.
    pub fn level0(&self) -> &LevelGraph {
        &self.level0
    }

    /// Coarse operator `A_l` for `l ≥ 1`.
    pub fn coarse_matrix(&self, l: usize) -> &LevelMatrix {
        &self.matrices[l - 1]
    }

    /// Smoothed prolongator between levels `l` and `l + 1`.
    pub fn prolongator(&self, l: usize) -> &LevelMatrix {
        &self.p[l]
    }

    /// Fine → coarse aggregate map between levels `l` and `l + 1`.
    pub fn aggregate_of(&self, l: usize) -> &[u32] {
        &self.aggregate_of[l]
    }

    /// Iterations of the first solve after the last setup, once known.
    pub fn setup_iterations(&self) -> Option<u32> {
        self.setup_iterations
    }

    /// `true` when the next `update` will run a full setup.
    pub fn resetup_pending(&self) -> bool {
        self.resetup_requested || self.levels.is_empty()
    }

    fn level_ref(&self, l: usize) -> LevelRef<'_> {
        if l == 0 {
            LevelRef::Graph(&self.level0)
        } else {
            LevelRef::Csr(&self.matrices[l - 1])
        }
    }

    fn fill_level0(&mut self, q: &[f64]) {
        self.map.update_weights(q, &mut self.level0);
        if self.perturbation != 0.0 {
            let s = self.perturbation;
            for a in self.level0.anchor.iter_mut() {
                *a += s;
            }
        }
    }

    fn estimate_lambda(&mut self, l: usize) -> f64 {
        let n = self.levels_n(l);
        let w = &mut self.work[l];
        let host = &mut self.host_scratch[..n * 3];
        // `v` = b, `w` = r, `d` = d, zero = x.
        estimate_lambda_max(
            &self.backend,
            &self.levels[l],
            POWER_ITERATIONS,
            0x5EED_0000 + l as u64,
            host,
            &mut w.b,
            &mut w.r,
            &mut w.d,
            &mut w.x,
        )
    }

    fn levels_n(&self, l: usize) -> usize {
        if l == 0 {
            self.n
        } else {
            self.matrices[l - 1].n
        }
    }

    fn push_device_level(&mut self, level: DeviceLevel<B>, n: usize) {
        let kcycle = self.options.cycle == CycleKind::K;
        self.levels.push(level);
        self.work
            .push(LevelWork::new(&self.backend, n, self.precision, kcycle));
    }

    /// Full setup for `q`: aggregation, `λ_max`, smoothed `P`, frozen
    /// coarse patterns, device uploads, coarsest factorization.
    fn setup(&mut self, q: &[f64]) -> Result<(), TheseusError> {
        self.matrices.clear();
        self.p.clear();
        self.pt.clear();
        self.ap.clear();
        self.aggregate_of.clear();
        self.lambda_max.clear();
        self.levels.clear();
        self.transfers.clear();
        self.work.clear();
        self.level0_matrix = None;
        self.coarsest = None;

        self.fill_level0(q);
        let level0 = self.backend.upload_level(&self.level0, self.precision);
        self.push_device_level(DeviceLevel::Graph(level0), self.n);
        self.outer_level0 = (self.precision != Precision::F64)
            .then(|| DeviceLevel::Graph(self.backend.upload_level(&self.level0, Precision::F64)));
        let lambda0 = self.estimate_lambda(0);
        self.lambda_max.push(lambda0);

        let passes = usize::from(self.options.aggregation_passes.max(1));
        let coarsest_size = self.options.coarsest_size as usize;
        loop {
            let l = self.levels.len() - 1;
            let n_l = self.levels_n(l);
            if n_l < coarsest_size || self.levels.len() >= MAX_LEVELS {
                break;
            }
            let strength = match self.level_ref(l) {
                LevelRef::Graph(g) => StrengthGraph::from_level_graph(g),
                LevelRef::Csr(m) => StrengthGraph::from_matrix(m),
            };
            let (agg, nc) = aggregate(&strength, passes);
            if nc as f64 > MIN_COARSENING * n_l as f64 || nc == 0 {
                break;
            }
            let omega = 4.0 / (3.0 * self.lambda_max[l]);
            let a = self.level_ref(l);
            let p = smoothed_prolongator(a, &agg, nc, omega);
            let pt = p.transpose();
            let (mut ap, mut coarse) = galerkin_pattern(&pt, a, &p);
            galerkin_numeric(&pt, a, &p, &mut ap, &mut coarse, &self.scratch);
            let dev_p = self.backend.upload_csr(&p, self.precision);
            let dev_pt = self.backend.upload_csr(&pt, self.precision);
            let dev_coarse = self.backend.upload_csr(&coarse, self.precision);
            self.transfers.push(Transfer {
                p: dev_p,
                pt: dev_pt,
            });
            self.p.push(p);
            self.pt.push(pt);
            self.ap.push(ap);
            self.aggregate_of.push(agg);
            self.matrices.push(coarse);
            self.push_device_level(DeviceLevel::Csr(dev_coarse), nc);
            let lambda = self.estimate_lambda(l + 1);
            self.lambda_max.push(lambda);
        }

        // Coarsest factorization.
        let coarsest_matrix = match self.matrices.last() {
            Some(m) => m,
            None => {
                self.level0_matrix = Some(LevelMatrix::from_graph(&self.level0));
                self.level0_matrix.as_ref().unwrap()
            }
        };
        let nc = coarsest_matrix.n;
        self.coarsest = Some(CoarsestSolver::new(coarsest_matrix)?);
        let last = self.work.pop().expect("coarsest work");
        self.work.push(last.with_coarsest_scratch(nc));

        self.q_setup.clear();
        self.q_setup.extend_from_slice(q);
        self.q_lambda.clear();
        self.q_lambda.extend_from_slice(q);
        self.setup_iterations = None;
        self.slow_solves = 0;
        self.resetup_requested = false;
        self.counters.setups += 1;
        self.backend.sync();
        Ok(())
    }

    /// Numeric update on the frozen patterns.
    fn numeric_update(&mut self, q: &[f64]) -> Result<(), TheseusError> {
        self.fill_level0(q);
        if let DeviceLevel::Graph(l0) = &mut self.levels[0] {
            self.backend
                .update_level_weights(l0, &self.level0.weight, &self.level0.anchor);
        }
        if let Some(DeviceLevel::Graph(l0)) = &mut self.outer_level0 {
            self.backend
                .update_level_weights(l0, &self.level0.weight, &self.level0.anchor);
        }
        for l in 0..self.matrices.len() {
            let (done, rest) = self.matrices.split_at_mut(l);
            let a = if l == 0 {
                LevelRef::Graph(&self.level0)
            } else {
                LevelRef::Csr(&done[l - 1])
            };
            galerkin_numeric(
                &self.pt[l],
                a,
                &self.p[l],
                &mut self.ap[l],
                &mut rest[0],
                &self.scratch,
            );
            if let DeviceLevel::Csr(dev) = &mut self.levels[l + 1] {
                self.backend.update_csr_values(&rest[0], dev);
            }
        }
        let coarsest = self.coarsest.as_mut().expect("coarsest present");
        match self.matrices.last() {
            Some(m) => coarsest.update(m)?,
            None => {
                let m = self.level0_matrix.as_mut().expect("level-0 matrix");
                m.set_from_graph(&self.level0);
                coarsest.update(m)?;
            }
        }
        if max_relative_change(q, &self.q_lambda) > LAMBDA_REESTIMATE_DRIFT {
            for l in 0..self.levels.len() {
                let lambda = self.estimate_lambda(l);
                self.lambda_max[l] = lambda;
            }
            self.q_lambda.clear();
            self.q_lambda.extend_from_slice(q);
            self.counters.lambda_reestimates += 1;
        }
        self.counters.numeric_updates += 1;
        self.backend.sync();
        Ok(())
    }

    /// Split `self` into the immutable cycle view and the mutable pieces a
    /// preconditioner application needs (disjoint field borrows).
    fn split_for_cycle(&mut self) -> (Cycle<'_, B>, &DeviceLevel<B>, CycleScratch<'_, B>) {
        let Self {
            backend,
            options,
            precision,
            smoother,
            lambda_max,
            coarsest,
            levels,
            outer_level0,
            transfers,
            work,
            outer,
            x_dev,
            b_dev,
            host_scratch,
            ..
        } = self;
        let cycle = Cycle {
            backend,
            levels,
            transfers,
            lambda_max,
            smoother: *smoother,
            kind: options.cycle,
            ktol: KCYCLE_TOL,
            coarsest: coarsest.as_ref().expect("coarsest present"),
        };
        let scratch = CycleScratch {
            precision: *precision,
            work,
            outer,
            x_dev,
            b_dev,
            host: host_scratch,
        };
        let outer_a = outer_level0.as_ref().unwrap_or(&cycle.levels[0]);
        (cycle, outer_a, scratch)
    }

    /// One V-cycle (or K-cycle) on host vectors: `z = M r`. For tests of
    /// the preconditioner's symmetry and positivity.
    pub fn precondition(&mut self, r: &[f64], z: &mut [f64]) -> Result<(), TheseusError> {
        if self.levels.is_empty() {
            return Err(TheseusError::MissingFactorization);
        }
        assert_eq!(r.len(), self.n * 3);
        assert_eq!(z.len(), self.n * 3);
        let (cycle, _, scratch) = self.split_for_cycle();
        let be = cycle.backend;
        be.upload(r, scratch.x_dev);
        precondition_with(
            &cycle,
            scratch.precision,
            scratch.work,
            scratch.host,
            scratch.x_dev,
            scratch.b_dev,
        );
        be.download(scratch.b_dev, z);
        Ok(())
    }

    fn run_pcg(&mut self, req: &SolveRequest<'_>) -> Result<PcgOutcome, TheseusError> {
        let flexible = self.options.cycle == CycleKind::K;
        let (cycle, outer_a, scratch) = self.split_for_cycle();
        let CycleScratch {
            precision,
            work,
            outer,
            x_dev,
            b_dev,
            host,
        } = scratch;
        pcg::pcg(
            cycle.backend,
            outer_a,
            b_dev,
            x_dev,
            outer,
            req.tolerance,
            req.max_iterations,
            flexible,
            req.cancel,
            |r, z| precondition_with(&cycle, precision, work, host, r, z),
        )
    }

    /// Bytes of the host-side hierarchy and scratch.
    fn host_bytes(&self) -> u64 {
        let g = &self.level0;
        let mut total = g.adjacency.offsets.len() * 4
            + g.adjacency.edges.len() * 4
            + g.adjacency.sign.len()
            + g.adjacency.other.len() * 4
            + g.weight.len() * 8
            + g.anchor.len() * 8
            + self.map.fine_edge.len() * 4
            + self.map.anchor_offsets.len() * 4
            + self.map.anchor_edges.len() * 4;
        total += self.matrices.iter().map(LevelMatrix::bytes).sum::<usize>();
        total += self.p.iter().map(LevelMatrix::bytes).sum::<usize>();
        total += self.pt.iter().map(LevelMatrix::bytes).sum::<usize>();
        total += self.ap.iter().map(LevelMatrix::bytes).sum::<usize>();
        total += self.aggregate_of.iter().map(|a| a.len() * 4).sum::<usize>();
        total += self.level0_matrix.as_ref().map_or(0, LevelMatrix::bytes);
        total += self.coarsest.as_ref().map_or(0, CoarsestSolver::bytes);
        total += self.scratch.bytes();
        total += (self.host_scratch.len() + self.q_setup.len() + self.q_lambda.len()) * 8;
        for w in &self.work {
            total += (w.host_b.len() + w.host_x.len() + w.solve_work.len()) * 8;
        }
        total as u64
    }

    /// Bytes of the backend buffers (vectors and uploaded operators).
    fn buffer_bytes(&self) -> u64 {
        let be = &self.backend;
        let sz = |buf: &B::Buf, precision: Precision| (be.len(buf) * precision.size_of()) as u64;
        let mut total = sz(&self.outer.r, Precision::F64)
            + sz(&self.outer.z, Precision::F64)
            + sz(&self.outer.p, Precision::F64)
            + sz(&self.outer.ap, Precision::F64)
            + sz(&self.x_dev, Precision::F64)
            + sz(&self.b_dev, Precision::F64);
        for w in &self.work {
            for buf in [&w.b, &w.x, &w.r, &w.d, &w.c1, &w.v1, &w.v2, &w.rt] {
                total += sz(buf, self.precision);
            }
        }
        // Uploaded operators: pattern (u32) + values in the buffer precision.
        let s = self.precision.size_of();
        let g = &self.level0;
        total += (g.adjacency.offsets.len() * 4
            + g.adjacency.edges.len() * 4
            + g.adjacency.other.len() * 4
            + (g.weight.len() + g.anchor.len() + g.n) * s) as u64;
        if self.outer_level0.is_some() {
            total += (g.adjacency.offsets.len() * 4
                + g.adjacency.edges.len() * 4
                + g.adjacency.other.len() * 4
                + (g.weight.len() + g.anchor.len() + g.n) * 8) as u64;
        }
        for m in self.matrices.iter().chain(&self.p).chain(&self.pt) {
            total += (m.row_ptr.len() * 4
                + m.col_idx.len() * 4
                + (m.values.len() + m.diag.len()) * s) as u64;
        }
        total
    }
}

impl<B> LinearSystemSolver for AmgSolver<B>
where
    B: Backend + Send,
    B::Buf: Send,
    B::LevelGraphBuf: Send,
    B::AggBuf: Send,
    B::CsrBuf: Send,
{
    fn update(&mut self, q: &[f64]) -> Result<(), TheseusError> {
        if q.len() != self.num_edges {
            return Err(TheseusError::Shape(format!(
                "AmgSolver::update: q has {} entries, topology has {} edges",
                q.len(),
                self.num_edges
            )));
        }
        let start = Instant::now();
        let needs_setup = self.levels.is_empty()
            || self.resetup_requested
            || max_relative_change(q, &self.q_setup) > RESETUP_DRIFT;
        if needs_setup {
            self.setup(q)?;
        } else {
            self.numeric_update(q)?;
        }
        self.last_update_ms += start.elapsed().as_secs_f64() * 1e3;
        Ok(())
    }

    fn solve(&mut self, req: SolveRequest<'_>, x: &mut [f64]) -> Result<SolveStats, TheseusError> {
        let n3 = self.n * 3;
        if req.rhs.len() != n3 || x.len() != n3 {
            return Err(TheseusError::Shape(format!(
                "AmgSolver::solve: expected rhs and x of length {n3} (n = {} × 3), got {} and {}",
                self.n,
                req.rhs.len(),
                x.len()
            )));
        }
        if let Some(x0) = req.x0 {
            if x0.len() != n3 {
                return Err(TheseusError::Shape(format!(
                    "AmgSolver::solve: x0 has length {}, expected {n3}",
                    x0.len()
                )));
            }
        }
        if self.levels.is_empty() {
            return Err(TheseusError::MissingFactorization);
        }
        let start = Instant::now();
        self.backend.upload(req.rhs, &mut self.b_dev);
        match req.x0 {
            Some(x0) => self.backend.upload(x0, &mut self.x_dev),
            None => self.backend.zero(&mut self.x_dev),
        }
        let outcome = self.run_pcg(&req)?;
        self.backend.download(&self.x_dev, x);
        self.counters.solves += 1;
        let solve_ms = start.elapsed().as_secs_f64() * 1e3;

        // Iteration-based re-setup rule.
        let iters = outcome.iterations.iter().copied().max().unwrap_or(0);
        match self.setup_iterations {
            None => self.setup_iterations = Some(iters),
            Some(base) => {
                if iters > RESETUP_ITERATION_FACTOR * base.max(1) {
                    self.slow_solves += 1;
                    if self.slow_solves >= 2 {
                        self.resetup_requested = true;
                    }
                } else {
                    self.slow_solves = 0;
                }
            }
        }

        let stats = SolveStats {
            iterations: outcome.iterations,
            relative_residual: outcome.relative_residual,
            converged: outcome.converged,
            setup_ms: std::mem::take(&mut self.last_update_ms),
            solve_ms,
            backend: self.kind(),
        };
        if !outcome.converged {
            return Err(TheseusError::IterativeSolverDidNotConverge {
                iterations: iters,
                relative_residual: outcome
                    .relative_residual
                    .iter()
                    .copied()
                    .fold(0.0, f64::max),
                kind: self.kind(),
            });
        }
        if !x.iter().all(|v| v.is_finite()) {
            return Err(TheseusError::Solver(
                crate::linear_solver::direct::NON_FINITE_SOLUTION_MSG.into(),
            ));
        }
        Ok(stats)
    }

    fn kind(&self) -> LinearSolverKind {
        match self.backend.handle() {
            BackendHandle::Cpu => LinearSolverKind::IterativeCpu,
            BackendHandle::Gpu => LinearSolverKind::IterativeGpu,
        }
    }

    fn memory_bytes(&self) -> MemoryReport {
        let buffers = self.buffer_bytes();
        match self.backend.handle() {
            BackendHandle::Cpu => MemoryReport {
                host_bytes: self.host_bytes() + buffers,
                device_bytes: 0,
            },
            BackendHandle::Gpu => MemoryReport {
                host_bytes: self.host_bytes(),
                device_bytes: buffers,
            },
        }
    }
}
