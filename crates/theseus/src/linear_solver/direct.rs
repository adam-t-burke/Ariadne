//! Direct sparse solver: the existing faer Cholesky / LDLᵀ path behind
//! [`LinearSystemSolver`].
//!
//! Two layers live here:
//!
//! * free functions ([`build_system_pattern`], [`refactor`]) that hold the
//!   direct-path logic once, used both by `FdmCache` / `fdm::factor_and_solve`
//!   (the production path) and by [`DirectSolver`];
//! * [`DirectSolver`], a self-contained owner of `A(q)`, its factorization
//!   and the solve workspace, implementing the trait so that the FDM code can
//!   be dispatched through `Box<dyn LinearSystemSolver>`.
//!
//! Nothing here changes the numerical behaviour of the direct path: the
//! assembly gather, the Cholesky → LDLᵀ fallback and the specialised
//! three-column triangular solve are the same code the cache has always run.

use super::{LinearSolverKind, LinearSystemSolver, MemoryReport, SolveRequest, SolveStats};
use crate::sparse::SparseColMatOwned;
use crate::types::{
    find_nz_index, Bounds, Factorization, FactorizationStrategy, NetworkTopology, QToNz,
    TheseusError,
};
use std::sync::atomic::Ordering;
use std::time::Instant;

/// Message of the non-finite-solution error, shared with `fdm::factor_and_solve`.
pub(crate) const NON_FINITE_SOLUTION_MSG: &str =
    "FDM linear solve produced non-finite solution (singular or ill-conditioned equilibrium matrix). \
     Check network connectivity, supports, and initial force densities.";

/// Build the sparsity pattern of `A = Cnᵀ Cn` over the free nodes and the
/// gather map from force densities to its nonzeros.
///
/// The returned matrix has zero values; call [`assemble_values`] to fill it
/// for a given `q`. Fails with [`TheseusError::SparsityMismatch`] if the
/// product pattern is inconsistent with the incidence (which cannot happen
/// for a well-formed topology).
pub fn build_system_pattern(
    topo: &NetworkTopology,
) -> Result<(SparseColMatOwned, QToNz), TheseusError> {
    let ne = topo.num_edges;
    let nn_free = topo.free_node_indices.len();

    // A's sparsity pattern from Cn^T * Cn.
    let cn = &topo.free_incidence; // ne × nn_free
    let cn_t = cn.transpose();
    let a_matrix =
        SparseColMatOwned::sparse_times_sparse(&cn_t, cn).map_err(TheseusError::Solver)?;

    // For each edge k, find which free nodes it touches in Cn, then map those
    // (n1, n2) pairs to indices in a_matrix.values.
    let mut edge_to_free_nodes: Vec<Vec<(usize, f64)>> = vec![Vec::new(); ne];
    for col in 0..nn_free {
        let start = cn.col_ptrs[col] as usize;
        let end_ = cn.col_ptrs[col + 1] as usize;
        for idx in start..end_ {
            let row = cn.row_indices[idx] as usize;
            let val = cn.values[idx];
            edge_to_free_nodes[row].push((col, val));
        }
    }

    let mut q_to_nz_entries: Vec<Vec<(usize, f64)>> = vec![Vec::new(); ne];
    for k in 0..ne {
        let nodes = &edge_to_free_nodes[k];
        for &(n1, v1) in nodes {
            for &(n2, v2) in nodes {
                let nz_idx = find_nz_index(&a_matrix.col_ptrs, &a_matrix.row_indices, n1, n2)
                    .ok_or(TheseusError::SparsityMismatch {
                        edge: k,
                        row: n1,
                        col: n2,
                    })?;
                q_to_nz_entries[k].push((nz_idx, v1 * v2));
            }
        }
    }

    let q_to_nz = QToNz::from_edge_entries(&q_to_nz_entries, a_matrix.values.len());
    Ok((a_matrix, q_to_nz))
}

/// Fill `values` (the `nzval` array of `A`) from `q` through the gather map:
/// `A = Cnᵀ diag(q) Cn`. Same summation order as `fdm::assemble_a`.
pub fn assemble_values(map: &QToNz, q: &[f64], values: &mut [f64]) {
    for (nz, value) in values.iter_mut().enumerate() {
        let range = map.nz_offsets[nz]..map.nz_offsets[nz + 1];
        *value = map.edge[range.clone()]
            .iter()
            .zip(&map.coeff[range])
            .map(|(&k, &c)| q[k as usize] * c)
            .sum();
    }
}

/// Factor `a` into `factorization`, reusing the symbolic analysis when one
/// exists, with the Cholesky → LDLᵀ fallback of the direct path.
///
/// On the first call `factorization` is `None` and a fresh symbolic + numeric
/// factorization with `strategy` is created. Later calls refactor numerically
/// in place. If the Cholesky path fails (the matrix is no longer SPD because
/// `q` drifted), `strategy` is switched to `LDL` permanently and the
/// factorization is rebuilt from scratch.
pub fn refactor(
    factorization: &mut Option<Factorization>,
    strategy: &mut FactorizationStrategy,
    a: &SparseColMatOwned,
    stack: &mut dyn_stack::GlobalPodBuffer,
) -> Result<(), TheseusError> {
    let mut need_ldl_fallback = false;

    match factorization {
        Some(fac) => {
            if let Err(e) = fac.update(a, stack) {
                if fac.strategy() == FactorizationStrategy::Cholesky {
                    need_ldl_fallback = true;
                } else {
                    return Err(e);
                }
            }
        }
        None => match Factorization::new(a, *strategy, stack) {
            Ok(fac) => {
                *factorization = Some(fac);
            }
            Err(_e) if *strategy == FactorizationStrategy::Cholesky => {
                need_ldl_fallback = true;
            }
            Err(e) => return Err(e),
        },
    }

    if need_ldl_fallback {
        *strategy = FactorizationStrategy::LDL;
        *factorization = None;
        *factorization = Some(Factorization::new(a, FactorizationStrategy::LDL, stack)?);
    }

    Ok(())
}

/// `Err(Solver)` with the direct path's message if any entry is not finite.
pub(crate) fn check_finite(x: &[f64]) -> Result<(), TheseusError> {
    if x.iter().all(|v| v.is_finite()) {
        Ok(())
    } else {
        Err(TheseusError::Solver(NON_FINITE_SOLUTION_MSG.into()))
    }
}

// ─────────────────────────────────────────────────────────────
//  DirectSolver
// ─────────────────────────────────────────────────────────────

/// The sparse direct path (`faer` Cholesky / LDLᵀ) as a [`LinearSystemSolver`].
///
/// Owns `A(q)`, its `q → nzval` gather map, the factorization and the
/// triangular-solve workspace for one topology. `update(q)` reassembles and
/// refactors (falling back from Cholesky to LDLᵀ if `A` stops being SPD, as
/// the cache path does); `solve` runs the specialised three-column
/// triangular solve of [`crate::factor_solve`].
pub struct DirectSolver {
    a_matrix: SparseColMatOwned,
    q_to_nz: QToNz,
    factorization: Option<Factorization>,
    strategy: FactorizationStrategy,
    /// Diagonal shift added to `A` before factoring (the `perturbation`
    /// argument of `fdm::factor_and_solve`); `0.0` = none.
    perturbation: f64,
    num_edges: usize,
    factor_stack: dyn_stack::GlobalPodBuffer,
    solve_workspace: Vec<f64>,
    last_update_ms: f64,
}

impl std::fmt::Debug for DirectSolver {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("DirectSolver")
            .field("n", &self.a_matrix.nrows)
            .field("nnz", &self.a_matrix.nnz())
            .field("num_edges", &self.num_edges)
            .field("strategy", &self.strategy)
            .field("factorization", &self.factorization)
            .field("perturbation", &self.perturbation)
            .finish()
    }
}

impl DirectSolver {
    /// Build the pattern of `A` for `topology` and choose `strategy`
    /// (Cholesky requires every `q > 0`; use LDLᵀ otherwise). No numeric
    /// work happens until [`LinearSystemSolver::update`].
    pub fn new(
        topology: &NetworkTopology,
        strategy: FactorizationStrategy,
    ) -> Result<Self, TheseusError> {
        let (a_matrix, q_to_nz) = build_system_pattern(topology)?;
        let n = a_matrix.nrows;
        Ok(Self {
            a_matrix,
            q_to_nz,
            factorization: None,
            strategy,
            perturbation: 0.0,
            num_edges: topology.num_edges,
            factor_stack: dyn_stack::GlobalPodBuffer::new(dyn_stack::StackReq::empty()),
            solve_workspace: vec![0.0; n * 6],
            last_update_ms: 0.0,
        })
    }

    /// [`Self::new`] with the strategy chosen from `bounds` as the cache does
    /// ([`FactorizationStrategy::from_bounds`]).
    pub fn from_bounds(topology: &NetworkTopology, bounds: &Bounds) -> Result<Self, TheseusError> {
        Self::new(topology, FactorizationStrategy::from_bounds(bounds))
    }

    /// Diagonal shift applied to `A` on the next `update`.
    pub fn set_perturbation(&mut self, perturbation: f64) {
        self.perturbation = perturbation;
    }

    /// Diagonal shift currently applied on `update`.
    pub fn perturbation(&self) -> f64 {
        self.perturbation
    }

    /// Current factorization strategy (may have fallen back to LDLᵀ).
    pub fn strategy(&self) -> FactorizationStrategy {
        self.strategy
    }

    /// Number of free nodes.
    pub fn num_nodes(&self) -> usize {
        self.a_matrix.nrows
    }

    /// The assembled `A(q)` from the last `update` (zero before the first).
    pub fn a_matrix(&self) -> &SparseColMatOwned {
        &self.a_matrix
    }

    /// The factorization from the last `update`, if any.
    pub fn factorization(&self) -> Option<&Factorization> {
        self.factorization.as_ref()
    }
}

impl LinearSystemSolver for DirectSolver {
    fn update(&mut self, q: &[f64]) -> Result<(), TheseusError> {
        if q.len() != self.num_edges {
            return Err(TheseusError::Shape(format!(
                "DirectSolver::update: q has {} entries, topology has {} edges",
                q.len(),
                self.num_edges
            )));
        }
        let start = Instant::now();
        assemble_values(&self.q_to_nz, q, &mut self.a_matrix.values);
        if self.perturbation > 0.0 {
            self.a_matrix.add_diagonal(self.perturbation);
        }
        refactor(
            &mut self.factorization,
            &mut self.strategy,
            &self.a_matrix,
            &mut self.factor_stack,
        )?;
        self.last_update_ms = start.elapsed().as_secs_f64() * 1e3;
        Ok(())
    }

    fn solve(&mut self, req: SolveRequest<'_>, x: &mut [f64]) -> Result<SolveStats, TheseusError> {
        let n = self.a_matrix.nrows;
        if req.rhs.len() != n * 3 || x.len() != n * 3 {
            return Err(TheseusError::Shape(format!(
                "DirectSolver::solve: expected rhs and x of length {} (n = {n} × 3), got {} and {}",
                n * 3,
                req.rhs.len(),
                x.len()
            )));
        }
        if req.cancel.is_some_and(|flag| flag.load(Ordering::Acquire)) {
            return Err(TheseusError::Cancelled);
        }
        let fac = self
            .factorization
            .as_ref()
            .ok_or(TheseusError::MissingFactorization)?;

        let start = Instant::now();
        fac.solve_slices::<3>(req.rhs, x, &mut self.solve_workspace);
        let solve_ms = start.elapsed().as_secs_f64() * 1e3;
        check_finite(x)?;

        Ok(SolveStats {
            iterations: [0; 3],
            relative_residual: [0.0; 3],
            converged: true,
            setup_ms: std::mem::take(&mut self.last_update_ms),
            solve_ms,
            backend: LinearSolverKind::Direct,
        })
    }

    fn kind(&self) -> LinearSolverKind {
        LinearSolverKind::Direct
    }

    fn memory_bytes(&self) -> MemoryReport {
        let f64s = |n: usize| (n * std::mem::size_of::<f64>()) as u64;
        let u32s = |n: usize| (n * std::mem::size_of::<u32>()) as u64;
        let a = &self.a_matrix;
        let mut host_bytes =
            u32s(a.col_ptrs.len()) + u32s(a.row_indices.len()) + f64s(a.values.len());
        host_bytes += (self.q_to_nz.nz_offsets.len() * std::mem::size_of::<usize>()) as u64
            + u32s(self.q_to_nz.edge.len())
            + f64s(self.q_to_nz.coeff.len());
        host_bytes += f64s(self.solve_workspace.len()) + self.factor_stack.len() as u64;
        if let Some(fac) = &self.factorization {
            host_bytes += f64s(fac.len_values());
        }
        MemoryReport {
            host_bytes,
            device_bytes: 0,
        }
    }
}
