use ndarray::Array2;
use sprs::{CsMat, FillInReduction, SymmetryCheck};
use sprs_ldl::{Ldl, LdlNumeric};
use std::fmt;
use std::fmt::Debug;

// ─────────────────────────────────────────────────────────────
//  Error type
// ─────────────────────────────────────────────────────────────

/// Unified error type for all fallible operations in the crate.
///
/// Every function in the public Rust API returns `Result<T, TheseusError>`
/// instead of panicking.  The FFI layer translates these into integer
/// return codes + a thread-local error message.
#[derive(Debug)]
pub enum TheseusError {
    /// Linear algebra failure (singular / not-SPD matrix, etc.).
    Linalg(sprs::errors::LinalgError),
    /// Sparsity pattern is inconsistent (should never happen after
    /// a correct `FdmCache::new`).
    SparsityMismatch { edge: usize, row: usize, col: usize },
    /// The factorization has not been computed yet.
    MissingFactorization,
    /// Argmin solver returned an error.
    Solver(String),
    /// Shape mismatch in input data.
    Shape(String),
    /// Optimization was cancelled by the caller via the progress callback.
    Cancelled,
}

impl fmt::Display for TheseusError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Linalg(e) => write!(f, "linear algebra error: {e}"),
            Self::SparsityMismatch { edge, row, col } =>
                write!(f, "sparsity pattern mismatch: edge {edge}, ({row},{col}) not in A"),
            Self::MissingFactorization =>
                write!(f, "factorization not computed (call solve_fdm first)"),
            Self::Solver(msg) => write!(f, "solver error: {msg}"),
            Self::Shape(msg) => write!(f, "shape error: {msg}"),
            Self::Cancelled => write!(f, "optimization cancelled by user"),
        }
    }
}

impl std::error::Error for TheseusError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Linalg(e) => Some(e),
            _ => None,
        }
    }
}

impl From<sprs::errors::LinalgError> for TheseusError {
    fn from(e: sprs::errors::LinalgError) -> Self {
        Self::Linalg(e)
    }
}

impl From<argmin::core::Error> for TheseusError {
    fn from(e: argmin::core::Error) -> Self {
        Self::Solver(e.to_string())
    }
}

// ─────────────────────────────────────────────────────────────
//  Constants
// ─────────────────────────────────────────────────────────────

pub const DEFAULT_BARRIER_SHARPNESS: f64 = 10.0;

// ─────────────────────────────────────────────────────────────
//  Objective trait  (extensible — implement for custom objectives)
// ─────────────────────────────────────────────────────────────

/// Trait for form-finding objectives.
///
/// Implement `loss` and `accumulate_gradient` to add custom objectives.
/// The gradient method must accumulate into `cache.grad_x` (for implicit
/// adjoint contributions) and/or `cache.grad_q` (for explicit q gradients).
pub trait ObjectiveTrait: Debug + Send + Sync {
    /// Scalar loss contribution from this objective.
    fn loss(&self, snap: &GeometrySnapshot) -> f64;

    /// Accumulate dJ/dx̂ into `cache.grad_x` and explicit dJ/dq into
    /// `cache.grad_q`.  Called before the adjoint solve.
    fn accumulate_gradient(
        &self,
        cache: &mut FdmCache,
        problem: &Problem,
    );

    /// Weight of this objective (used for display/debugging).
    fn weight(&self) -> f64;
}

// ─────────────────────────────────────────────────────────────
//  Built-in objective structs  (13 types from the Julia code)
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct TargetXYZ {
    pub weight: f64,
    pub node_indices: Vec<usize>,
    pub target: Array2<f64>, // n × 3
}

#[derive(Debug, Clone)]
pub struct TargetXY {
    pub weight: f64,
    pub node_indices: Vec<usize>,
    pub target: Array2<f64>,
}

/// Target positions on an arbitrary plane. Origin and axes are in world coordinates;
/// axes should be unit and orthogonal (e.g. Rhino plane Origin, XAxis, YAxis).
#[derive(Debug, Clone)]
pub struct TargetPlane {
    pub weight: f64,
    pub node_indices: Vec<usize>,
    pub target: Array2<f64>, // n × 3 world positions
    pub origin: [f64; 3],
    pub x_axis: [f64; 3],
    pub y_axis: [f64; 3],
}

#[derive(Debug, Clone)]
pub struct TargetLength {
    pub weight: f64,
    pub edge_indices: Vec<usize>,
    pub target: Vec<f64>,
}

#[derive(Debug, Clone)]
pub struct LengthVariation {
    pub weight: f64,
    pub edge_indices: Vec<usize>,
    pub sharpness: f64,
}

#[derive(Debug, Clone)]
pub struct ForceVariation {
    pub weight: f64,
    pub edge_indices: Vec<usize>,
    pub sharpness: f64,
}

#[derive(Debug, Clone)]
pub struct SumForceLength {
    pub weight: f64,
    pub edge_indices: Vec<usize>,
}

#[derive(Debug, Clone)]
pub struct MinLength {
    pub weight: f64,
    pub edge_indices: Vec<usize>,
    pub threshold: Vec<f64>,
    pub sharpness: f64,
}

#[derive(Debug, Clone)]
pub struct MaxLength {
    pub weight: f64,
    pub edge_indices: Vec<usize>,
    pub threshold: Vec<f64>,
    pub sharpness: f64,
}

#[derive(Debug, Clone)]
pub struct MinForce {
    pub weight: f64,
    pub edge_indices: Vec<usize>,
    pub threshold: Vec<f64>,
    pub sharpness: f64,
}

#[derive(Debug, Clone)]
pub struct MaxForce {
    pub weight: f64,
    pub edge_indices: Vec<usize>,
    pub threshold: Vec<f64>,
    pub sharpness: f64,
}

#[derive(Debug, Clone)]
pub struct RigidSetCompare {
    pub weight: f64,
    pub node_indices: Vec<usize>,
    pub target: Array2<f64>,
}

#[derive(Debug, Clone)]
pub struct ReactionDirection {
    pub weight: f64,
    pub anchor_indices: Vec<usize>,
    pub target_directions: Array2<f64>, // n × 3, unit rows
}

#[derive(Debug, Clone)]
pub struct ReactionDirectionMagnitude {
    pub weight: f64,
    pub anchor_indices: Vec<usize>,
    pub target_directions: Array2<f64>,
    pub target_magnitudes: Vec<f64>,
}

// ─────────────────────────────────────────────────────────────
//  Bounds
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct Bounds {
    pub lower: Vec<f64>,
    pub upper: Vec<f64>,
}

impl Bounds {
    pub fn default_for(num_edges: usize) -> Self {
        Self {
            lower: vec![1e-8; num_edges],
            upper: vec![f64::INFINITY; num_edges],
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Solver / Tracing options
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct SolverOptions {
    pub absolute_tolerance: f64,
    pub relative_tolerance: f64,
    pub max_iterations: usize,
    pub report_frequency: usize,
    pub barrier_weight: f64,
    pub barrier_sharpness: f64,
}

impl Default for SolverOptions {
    fn default() -> Self {
        Self {
            absolute_tolerance: 1e-6,
            relative_tolerance: 1e-6,
            max_iterations: 500,
            report_frequency: 1,
            barrier_weight: 10.0,
            barrier_sharpness: DEFAULT_BARRIER_SHARPNESS,
        }
    }
}


// ─────────────────────────────────────────────────────────────
//  Network topology
// ─────────────────────────────────────────────────────────────

/// Compressed connectivity information built once from the incidence matrix.
#[derive(Debug, Clone)]
pub struct NetworkTopology {
    /// Full incidence matrix  (ne × nn)  with ±1 entries.
    pub incidence: CsMat<f64>,
    /// Free-node incidence    (ne × nn_free)
    pub free_incidence: CsMat<f64>,
    /// Fixed-node incidence   (ne × nn_fixed)
    pub fixed_incidence: CsMat<f64>,
    pub num_edges: usize,
    pub num_nodes: usize,
    pub free_node_indices: Vec<usize>,
    pub fixed_node_indices: Vec<usize>,
}

// ─────────────────────────────────────────────────────────────
//  Anchor info  (variable / fixed supports)
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct AnchorInfo {
    pub variable_indices: Vec<usize>,
    pub fixed_indices: Vec<usize>,
    pub reference_positions: Array2<f64>,       // n_fixed × 3
    pub initial_variable_positions: Array2<f64>, // n_var × 3
}

impl AnchorInfo {
    /// All anchors fixed – no movable supports.
    pub fn all_fixed(reference_positions: Array2<f64>) -> Self {
        let n = reference_positions.nrows();
        Self {
            variable_indices: Vec::new(),
            fixed_indices: (0..n).collect(),
            reference_positions,
            initial_variable_positions: Array2::zeros((0, 3)),
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Problem definition  (immutable after construction)
// ─────────────────────────────────────────────────────────────

#[derive(Debug)]
pub struct Problem {
    pub topology: NetworkTopology,
    pub free_node_loads: Array2<f64>,  // nn_free × 3
    pub fixed_node_positions: Array2<f64>, // n_fixed × 3  (reference)
    pub anchors: AnchorInfo,
    pub objectives: Vec<Box<dyn ObjectiveTrait>>,
    pub bounds: Bounds,
    pub solver: SolverOptions,
}

// ─────────────────────────────────────────────────────────────
//  Sparsity mapping  q_k  →  A.data[] indices
// ─────────────────────────────────────────────────────────────

/// Pre-computed contribution of edge `k` to the CSC `nzval` array of A.
#[derive(Debug, Clone)]
pub struct QToNz {
    /// For each edge k: list of (nz_index_in_A_data, coefficient)
    pub entries: Vec<Vec<(usize, f64)>>,
}

// ─────────────────────────────────────────────────────────────
//  Factorization strategy
// ─────────────────────────────────────────────────────────────

/// Adaptive factorization for A(q) = Cn^T diag(q) Cn.
/// Cholesky when bounds guarantee sign-definiteness; LDL for mixed sign.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FactorizationStrategy {
    /// All q_k > 0  (or all < 0):  A is SPD → Cholesky.
    /// Uses AMD fill-in reduction for better sparsity in L.
    Cholesky,
    /// Mixed sign q allowed:  A is symmetric indefinite → LDL.
    LDL,
}

impl FactorizationStrategy {
    /// Choose strategy from the bounds on q.
    pub fn from_bounds(bounds: &Bounds) -> Self {
        let all_positive = bounds.lower.iter().all(|&lb| lb > 0.0);
        let all_negative = bounds.upper.iter().all(|&ub| ub < 0.0);
        if all_positive || all_negative {
            Self::Cholesky
        } else {
            Self::LDL
        }
    }
}

/// Holds a numeric LDL^T (or Cholesky) factorization.
///
/// Both variants use `sprs-ldl`'s `LdlNumeric` internally.
/// The Cholesky path uses AMD fill-in reduction and validates D > 0.
/// The LDL path allows indefinite D.
pub enum Factorization {
    /// SPD path: AMD-ordered, D > 0 validated
    Cholesky(LdlNumeric<f64, usize>),
    /// Indefinite path: no sign constraint on D
    Ldl(LdlNumeric<f64, usize>),
}

impl std::fmt::Debug for Factorization {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Cholesky(_) => write!(f, "Factorization::Cholesky(...)"),
            Self::Ldl(_) => write!(f, "Factorization::Ldl(...)"),
        }
    }
}

impl Factorization {
    /// Create an initial factorization from A and the chosen strategy.
    pub fn new(a: sprs::CsMatView<f64>, strategy: FactorizationStrategy) -> Result<Self, sprs::errors::LinalgError> {
        match strategy {
            FactorizationStrategy::Cholesky => {
                let ldl = Ldl::new()
                    .fill_in_reduction(FillInReduction::ReverseCuthillMcKee)
                    .check_symmetry(SymmetryCheck::DontCheckSymmetry)
                    .numeric(a)?;
                // Validate positive-definiteness: all diagonal D entries > 0
                for (i, &di) in ldl.d().iter().enumerate() {
                    if di <= 0.0 {
                        return Err(sprs::errors::LinalgError::SingularMatrix(
                            sprs::errors::SingularMatrixInfo {
                                index: i,
                                reason: "D <= 0 in Cholesky factorization (not SPD)",
                            },
                        ));
                    }
                }
                Ok(Self::Cholesky(ldl))
            }
            FactorizationStrategy::LDL => {
                let ldl = Ldl::new()
                    .fill_in_reduction(FillInReduction::ReverseCuthillMcKee)
                    .check_symmetry(SymmetryCheck::DontCheckSymmetry)
                    .numeric(a)?;
                Ok(Self::Ldl(ldl))
            }
        }
    }

    /// Re-factor with updated numeric values (same sparsity pattern).
    pub fn update(&mut self, a: sprs::CsMatView<f64>) -> Result<(), sprs::errors::LinalgError> {
        match self {
            Self::Cholesky(ldl) => {
                ldl.update(a)?;
                for (i, &di) in ldl.d().iter().enumerate() {
                    if di <= 0.0 {
                        return Err(sprs::errors::LinalgError::SingularMatrix(
                            sprs::errors::SingularMatrixInfo {
                                index: i,
                                reason: "D <= 0 in Cholesky re-factor (not SPD)",
                            },
                        ));
                    }
                }
                Ok(())
            }
            Self::Ldl(ldl) => {
                ldl.update(a)?;
                Ok(())
            }
        }
    }

    /// Solve A x = rhs using the stored factorization.
    pub fn solve(&self, rhs: &[f64]) -> Vec<f64> {
        match self {
            Self::Cholesky(ldl) | Self::Ldl(ldl) => ldl.solve(rhs),
        }
    }

    /// The strategy this factorization was built with.
    pub fn strategy(&self) -> FactorizationStrategy {
        match self {
            Self::Cholesky(_) => FactorizationStrategy::Cholesky,
            Self::Ldl(_) => FactorizationStrategy::LDL,
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Pre-allocated solver cache
// ─────────────────────────────────────────────────────────────

/// All mutable workspace for the forward solve, adjoint, and gradient
/// accumulation.  Built once from a [`Problem`], reused across iterations.
#[derive(Debug)]
pub struct FdmCache {
    // ── Sparse system ──────────────────────────────────────
    /// System matrix A = Cn^T diag(q) Cn  (CSC, nn_free × nn_free).
    /// Sparsity pattern is fixed; values are updated in-place each iteration.
    pub a_matrix: CsMat<f64>,

    /// Numeric factorization — Cholesky (SPD) or LDL (indefinite).
    /// Created on first factor, reused via `.update()` thereafter.
    pub factorization: Option<Factorization>,

    pub q_to_nz: QToNz,

    /// Start / end node of each edge (global node indices, 0-based)
    pub edge_starts: Vec<usize>,
    pub edge_ends: Vec<usize>,
    /// Global-node → free-index mapping  (`None` if fixed)
    pub node_to_free_idx: Vec<Option<usize>>,

    /// Cn  (ne × nn_free)  and  Cf  (ne × nn_fixed)  stored as CSC
    pub cn: CsMat<f64>,
    pub cf: CsMat<f64>,

    // ── Primal buffers ─────────────────────────────────────
    /// Free-node positions         (nn_free × 3, column-major)
    pub x: Array2<f64>,
    /// Adjoint variables           (nn_free × 3)
    pub lambda: Array2<f64>,
    /// dJ / d(free-node positions) (nn_free × 3)
    pub grad_x: Array2<f64>,
    /// Force densities
    pub q: Vec<f64>,
    /// dJ / dq
    pub grad_q: Vec<f64>,
    /// dJ / dNf  (all nodes × 3, only fixed rows used)
    pub grad_nf: Array2<f64>,

    // ── Derived geometry ───────────────────────────────────
    pub member_lengths: Vec<f64>,
    pub member_forces: Vec<f64>,
    pub reactions: Array2<f64>, // nn × 3

    // ── Intermediate RHS buffers ───────────────────────────
    pub cf_nf: Array2<f64>,    // ne × 3
    pub q_cf_nf: Array2<f64>,  // ne × 3
    pub pn: Array2<f64>,       // nn_free × 3  (copy of free-node loads)
    pub nf: Array2<f64>,       // nn × 3       (full node positions)
    pub nf_fixed: Array2<f64>, // nn_fixed × 3

    // ── RHS buffer (reusable for linear solve input) ───────
    pub rhs: Array2<f64>,      // nn_free × 3

    // ── Factorization ──────────────────────────────────────
    pub strategy: FactorizationStrategy,
}

impl FdmCache {
    /// Build a fully pre-allocated cache from a [`Problem`].
    ///
    /// Returns `Err` if the incidence sparsity pattern is inconsistent.
    pub fn new(problem: &Problem) -> Result<Self, TheseusError> {
        let topo = &problem.topology;
        let ne = topo.num_edges;
        let nn = topo.num_nodes;
        let nn_free = topo.free_node_indices.len();
        let nn_fixed = topo.fixed_node_indices.len();

        // ── 1. Build A's sparsity pattern from Cn^T * Cn ──
        let cn = &topo.free_incidence; // ne × nn_free
        let cn_t = cn.transpose_view().to_csc();
        // Symbolic Cn^T * Cn to get the pattern
        let a_template = &cn_t * cn;
        let a_matrix = a_template.to_csc();

        // ── 2. Build q_to_nz mapping ──────────────────────
        // For each edge k, find which free nodes it touches in Cn,
        // then map those (n1, n2) pairs to indices in a_matrix.data().
        let mut edge_to_free_nodes: Vec<Vec<(usize, f64)>> = vec![Vec::new(); ne];
        // Iterate columns of Cn (CSC: each column = a free node)
        let cn_csc = cn.to_csc();
        for col in 0..nn_free {
            let start = cn_csc.indptr().raw_storage()[col];
            let end_ = cn_csc.indptr().raw_storage()[col + 1];
            for idx in start..end_ {
                let row = cn_csc.indices()[idx]; // edge index
                let val = cn_csc.data()[idx];
                edge_to_free_nodes[row].push((col, val));
            }
        }

        let mut q_to_nz_entries: Vec<Vec<(usize, f64)>> = vec![Vec::new(); ne];
        for k in 0..ne {
            let nodes = &edge_to_free_nodes[k];
            for &(n1, v1) in nodes {
                for &(n2, v2) in nodes {
                    // Find nz index of (n1, n2) in CSC  [row=n1, col=n2]
                    let indptr = a_matrix.indptr();
                    let nz_idx = find_nz_index(indptr.raw_storage(), a_matrix.indices(), n1, n2)
                        .ok_or(TheseusError::SparsityMismatch { edge: k, row: n1, col: n2 })?;
                    q_to_nz_entries[k].push((nz_idx, v1 * v2));
                }
            }
        }

        // ── 3. Edge start / end from incidence ────────────
        let mut edge_starts = vec![0usize; ne];
        let mut edge_ends = vec![0usize; ne];
        let inc = &topo.incidence;
        let inc_csc = inc.to_csc();
        for col in 0..nn {
            let start = inc_csc.indptr().raw_storage()[col];
            let end_ = inc_csc.indptr().raw_storage()[col + 1];
            for idx in start..end_ {
                let row = inc_csc.indices()[idx];
                let val = inc_csc.data()[idx];
                if val == -1.0 {
                    edge_starts[row] = col;
                } else if val == 1.0 {
                    edge_ends[row] = col;
                }
            }
        }

        // ── 4. node_to_free_idx ───────────────────────────
        let mut node_to_free_idx = vec![None; nn];
        for (i, &node) in topo.free_node_indices.iter().enumerate() {
            node_to_free_idx[node] = Some(i);
        }

        // ── 5. Factorization strategy ─────────────────────
        let strategy = FactorizationStrategy::from_bounds(&problem.bounds);

        // ── 6. Pre-allocate all buffers ───────────────────
        let cf = topo.fixed_incidence.clone();
        let cn_owned = cn.clone();

        Ok(FdmCache {
            a_matrix,
            factorization: None,
            q_to_nz: QToNz { entries: q_to_nz_entries },
            edge_starts,
            edge_ends,
            node_to_free_idx,
            cn: cn_owned,
            cf,
            x: Array2::zeros((nn_free, 3)),
            lambda: Array2::zeros((nn_free, 3)),
            grad_x: Array2::zeros((nn_free, 3)),
            q: vec![0.0; ne],
            grad_q: vec![0.0; ne],
            grad_nf: Array2::zeros((nn, 3)),
            member_lengths: vec![0.0; ne],
            member_forces: vec![0.0; ne],
            reactions: Array2::zeros((nn, 3)),
            cf_nf: Array2::zeros((ne, 3)),
            q_cf_nf: Array2::zeros((ne, 3)),
            pn: problem.free_node_loads.clone(),
            nf: Array2::zeros((nn, 3)),
            nf_fixed: Array2::zeros((nn_fixed, 3)),
            rhs: Array2::zeros((nn_free, 3)),
            strategy,
        })
    }
}

// ─────────────────────────────────────────────────────────────
//  Geometry snapshot  (read-only view after forward solve)
// ─────────────────────────────────────────────────────────────

/// Immutable snapshot of the geometry after a forward FDM solve.
/// Views borrow from `FdmCache` buffers.
pub struct GeometrySnapshot<'a> {
    pub xyz_full: &'a Array2<f64>,     // nn × 3
    pub member_lengths: &'a [f64],
    pub member_forces: &'a [f64],
    pub reactions: &'a Array2<f64>,     // nn × 3
}

// ─────────────────────────────────────────────────────────────
//  Optimisation state  (mutable across iterations)
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct OptimizationState {
    pub force_densities: Vec<f64>,
    pub variable_anchor_positions: Array2<f64>, // n_var × 3
    pub loss_trace: Vec<f64>,
    pub iterations: usize,
}

impl OptimizationState {
    pub fn new(q: Vec<f64>, anchors: Array2<f64>) -> Self {
        Self {
            force_densities: q,
            variable_anchor_positions: anchors,
            loss_trace: Vec::new(),
            iterations: 0,
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Solver result  (returned from optimize)
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct SolverResult {
    pub q: Vec<f64>,
    pub anchor_positions: Array2<f64>,
    pub xyz: Array2<f64>,        // nn × 3
    pub member_lengths: Vec<f64>,
    pub member_forces: Vec<f64>,
    pub reactions: Array2<f64>,  // nn × 3
    pub loss_trace: Vec<f64>,
    pub iterations: usize,
    pub converged: bool,
    pub termination_reason: String,
}

// ─────────────────────────────────────────────────────────────
//  Helper: find nz index in CSC
// ─────────────────────────────────────────────────────────────

/// Given CSC indptr and indices arrays, find the position of element (row, col)
/// in the data array.  Returns `None` if the entry is not in the sparsity pattern.
pub fn find_nz_index(
    indptr: &[usize],
    indices: &[usize],
    row: usize,
    col: usize,
) -> Option<usize> {
    let start = indptr[col];
    let end_ = indptr[col + 1];
    for nz in start..end_ {
        if indices[nz] == row {
            return Some(nz);
        }
    }
    None
}
