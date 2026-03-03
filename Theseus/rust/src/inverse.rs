//! Inverse FDM solvers: find force densities q from a target geometry.
//!
//! Given target free-node positions, the FDM equilibrium per dimension d is:
//!
//!   Cn^T  diag(u_d)  q  =  Pn_d
//!
//! where u_d = C · N_target are the target member coordinate differences.
//! Stacking all 3 dimensions gives the sparse system  M q = p.
//!
//! Three solvers are provided:
//!
//!   - **Pseudoinverse L2** (Tikhonov-regularised normal equations):
//!       min ‖Mq − p‖²  (+λ‖q‖²)   →  single sparse LDL via `sprs-ldl`
//!
//!   - **Pseudoinverse L1** (iteratively reweighted least squares):
//!       min ‖Mq − p‖₁  (+λ‖q‖²)   →  IRLS: repeated weighted L2 solves
//!       More robust than L2 when a few equilibrium equations are outliers.
//!
//!   - **NNLS** (spectral projected gradient):
//!       min ‖Mq − p‖²  s.t. q ≥ 0   →  sparse matvecs only

use crate::types::{Problem, TheseusError};
use ndarray::Array2;
use sprs::{CsMat, TriMat};
use sprs_ldl::Ldl;

// ─────────────────────────────────────────────────────────────
//  Target member vectors
// ─────────────────────────────────────────────────────────────

/// Compute target member coordinate differences from full node positions.
///
/// Returns `u` (ne × 3) where `u[k,:] = N_target[end_k] − N_target[start_k]`.
fn compute_target_member_vectors(
    problem: &Problem,
    target_free_xyz: &Array2<f64>,
) -> Array2<f64> {
    let topo = &problem.topology;
    let ne = topo.num_edges;
    let nn = topo.num_nodes;

    // Build full node position array: free nodes from target, fixed from problem
    let mut nf = Array2::<f64>::zeros((nn, 3));

    for (i, &node) in topo.free_node_indices.iter().enumerate() {
        for d in 0..3 {
            nf[[node, d]] = target_free_xyz[[i, d]];
        }
    }
    for (i, &node) in topo.fixed_node_indices.iter().enumerate() {
        for d in 0..3 {
            nf[[node, d]] = problem.fixed_node_positions[[i, d]];
        }
    }

    // u = C * N  via incidence (CSC)
    let inc = topo.incidence.to_csc();
    let mut u = Array2::<f64>::zeros((ne, 3));
    for col in 0..nn {
        let start = inc.indptr().raw_storage()[col];
        let end_ = inc.indptr().raw_storage()[col + 1];
        for nz in start..end_ {
            let row = inc.indices()[nz];
            let val = inc.data()[nz];
            for d in 0..3 {
                u[[row, d]] += val * nf[[col, d]];
            }
        }
    }
    u
}

// ─────────────────────────────────────────────────────────────
//  Build sparse equilibrium matrix  M  (3·nn_free × ne)
// ─────────────────────────────────────────────────────────────

/// Build the sparse matrix M by vertically stacking Cn^T · diag(u_d) for d=x,y,z.
///
/// M has the same column-sparsity as Cn^T (at most 2 non-zeros per column per
/// dimension block → ≤ 6 per column total).  Uses COO → CSC for construction.
fn build_equilibrium_matrix(
    problem: &Problem,
    u: &Array2<f64>,
) -> (CsMat<f64>, Vec<f64>) {
    let topo = &problem.topology;
    let ne = topo.num_edges;
    let nn_free = topo.free_node_indices.len();
    let m_rows = 3 * nn_free;

    // Cn is (ne × nn_free).  Cn^T is (nn_free × ne).
    let cn_t = topo.free_incidence.transpose_view().to_csc();

    let mut tri = TriMat::new((m_rows, ne));

    // For each dimension d, M_d = Cn^T · diag(u_d): column j of Cn^T scaled by u[j,d]
    for d in 0..3usize {
        let row_offset = d * nn_free;
        for col in 0..ne {
            let scale = u[[col, d]];
            if scale == 0.0 {
                continue;
            }
            let start = cn_t.indptr().raw_storage()[col];
            let end_ = cn_t.indptr().raw_storage()[col + 1];
            for nz in start..end_ {
                let row = cn_t.indices()[nz];
                let val = cn_t.data()[nz];
                tri.add_triplet(row_offset + row, col, val * scale);
            }
        }
    }

    let m_mat = tri.to_csc();

    // p = [Pn_x; Pn_y; Pn_z]
    let loads = &problem.free_node_loads;
    let mut p = vec![0.0; m_rows];
    for d in 0..3 {
        let off = d * nn_free;
        for i in 0..nn_free {
            p[off + i] = loads[[i, d]];
        }
    }

    (m_mat, p)
}

// ─────────────────────────────────────────────────────────────
//  Augmented (saddle-point) system builder
// ─────────────────────────────────────────────────────────────

/// Build the augmented (saddle-point) system that avoids forming M^T M.
///
/// Instead of the normal equations `(M^T M + λI) q = M^T p`, assembles:
///
///   [ I    M  ] [r]   [rhs_top]
///   [ M^T  -λI] [q] = [   0   ]
///
/// The augmented matrix is symmetric indefinite with size (m+n)×(m+n)
/// where m = M.rows() and n = M.cols().  Non-zeros are O(nnz(M))
/// rather than O(nnz(M^T M)), which avoids the fill-in explosion
/// from the Gram product at large scales (>50k edges).
///
/// Returns `(K, rhs)` where K is CSC and `rhs = [rhs_top; 0]`.
fn build_augmented_system(
    m_mat: &CsMat<f64>,
    rhs_top: &[f64],
    regularization: f64,
) -> (CsMat<f64>, Vec<f64>) {
    let m = m_mat.rows();
    let n = m_mat.cols();
    let total = m + n;

    let m_csc = m_mat.to_csc();
    let nnz_m = m_csc.nnz();
    let estimated_nnz = m + 2 * nnz_m + n;

    let mut tri = TriMat::with_capacity((total, total), estimated_nnz);

    // Top-left: I  (m × m)
    for i in 0..m {
        tri.add_triplet(i, i, 1.0);
    }

    // Upper-right: M  and  lower-left: M^T  (symmetric pair)
    let indptr = m_csc.indptr();
    let raw = indptr.raw_storage();
    for col in 0..n {
        for nz in raw[col]..raw[col + 1] {
            let row = m_csc.indices()[nz];
            let val = m_csc.data()[nz];
            tri.add_triplet(row, m + col, val);
            tri.add_triplet(m + col, row, val);
        }
    }

    // Bottom-right: -λI  (n × n)
    for j in 0..n {
        tri.add_triplet(m + j, m + j, -regularization);
    }

    let k_mat = tri.to_csc();

    // RHS = [rhs_top; 0]
    let mut rhs = vec![0.0; total];
    rhs[..m].copy_from_slice(rhs_top);

    (k_mat, rhs)
}

// ─────────────────────────────────────────────────────────────
//  Implicit matvecs  (for NNLS — avoids forming M)
// ─────────────────────────────────────────────────────────────

/// Compute r = M·q  without forming M.
///
/// r_d = Cn^T · (u_d ⊙ q)  for each dimension d, concatenated.
fn apply_m(cn_t_csc: &CsMat<f64>, u: &Array2<f64>, q: &[f64], nn_free: usize) -> Vec<f64> {
    let ne = q.len();
    let mut r = vec![0.0; 3 * nn_free];

    for d in 0..3usize {
        let off = d * nn_free;
        // w = u_d ⊙ q
        // Cn^T · w  (Cn^T is nn_free × ne, CSC)
        for col in 0..ne {
            let w = u[[col, d]] * q[col];
            if w == 0.0 {
                continue;
            }
            let start = cn_t_csc.indptr().raw_storage()[col];
            let end_ = cn_t_csc.indptr().raw_storage()[col + 1];
            for nz in start..end_ {
                let row = cn_t_csc.indices()[nz];
                let val = cn_t_csc.data()[nz];
                r[off + row] += val * w;
            }
        }
    }
    r
}

/// Compute g = M^T · r  without forming M.
///
/// g_k = Σ_d  u[k,d] · (Cn · r_d)[k]
fn apply_mt(cn_csc: &CsMat<f64>, u: &Array2<f64>, r: &[f64], nn_free: usize) -> Vec<f64> {
    let ne = u.nrows();
    let mut g = vec![0.0; ne];

    for d in 0..3usize {
        let off = d * nn_free;
        let r_d = &r[off..off + nn_free];

        // v = Cn · r_d  (Cn is ne × nn_free, CSC)
        // Then g_k += u[k,d] * v_k
        for col in 0..nn_free {
            let rd_val = r_d[col];
            if rd_val == 0.0 {
                continue;
            }
            let start = cn_csc.indptr().raw_storage()[col];
            let end_ = cn_csc.indptr().raw_storage()[col + 1];
            for nz in start..end_ {
                let row = cn_csc.indices()[nz];
                let val = cn_csc.data()[nz];
                g[row] += u[[row, d]] * val * rd_val;
            }
        }
    }
    g
}

// ─────────────────────────────────────────────────────────────
//  Pseudoinverse  (Tikhonov-regularised sparse normal equations)
// ─────────────────────────────────────────────────────────────

/// Find force densities via pseudoinverse of the equilibrium system.
///
/// Solves  `(M^T M + λI) q = M^T p`  using sparse LDL factorisation.
pub fn solve_pseudoinverse(
    problem: &Problem,
    target_free_xyz: &Array2<f64>,
    regularization: f64,
) -> Result<Vec<f64>, TheseusError> {
    let ne = problem.topology.num_edges;

    let u = compute_target_member_vectors(problem, target_free_xyz);
    let (m_mat, p) = build_equilibrium_matrix(problem, &u);

    // G = M^T M  (ne × ne, sparse)
    let m_t = m_mat.transpose_view().to_csc();
    let m_csc = m_mat.to_csc();
    let mut g: CsMat<f64> = (&m_t * &m_csc).to_csc();

    // Add regularisation: G += λ I
    if regularization > 0.0 {
        add_diagonal(&mut g, regularization);
    }

    // h = M^T p  (sparse-dense matvec)
    let h = sparse_matvec(&m_t, &p);

    // Factorise G and solve
    let ldl = Ldl::new()
        .fill_in_reduction(sprs::FillInReduction::ReverseCuthillMcKee)
        .check_symmetry(sprs::SymmetryCheck::DontCheckSymmetry)
        .numeric(g.view())
        .map_err(|e| TheseusError::Linalg(e))?;

    let q = ldl.solve(&h);

    // Validate solution
    for (i, &v) in q.iter().enumerate() {
        if !v.is_finite() {
            return Err(TheseusError::Solver(format!(
                "pseudoinverse produced non-finite q at edge {i}; \
                 increase regularisation or check target geometry",
            )));
        }
    }

    // Ensure correct length (should always match, but be safe)
    if q.len() != ne {
        return Err(TheseusError::Shape(format!(
            "pseudoinverse returned {} q values, expected {ne}",
            q.len()
        )));
    }

    Ok(q)
}

// ─────────────────────────────────────────────────────────────
//  Pseudoinverse L2  (augmented saddle-point system)
// ─────────────────────────────────────────────────────────────

/// Find force densities via pseudoinverse using the augmented saddle-point system.
///
/// Mathematically equivalent to `solve_pseudoinverse` but avoids forming M^T M.
/// Instead factorises the larger but much sparser augmented system:
///
///   [ I    M  ] [r]   [p]
///   [ M^T  -λI] [q] = [0]
///
/// Asymptotically faster for large meshes (>50k edges) where the M^T M
/// fill-in explosion dominates runtime.  Requires `λ > 0`.
pub fn solve_pseudoinverse_augmented(
    problem: &Problem,
    target_free_xyz: &Array2<f64>,
    regularization: f64,
) -> Result<Vec<f64>, TheseusError> {
    if regularization <= 0.0 {
        return Err(TheseusError::Solver(
            "augmented system requires regularization > 0".into(),
        ));
    }

    let ne = problem.topology.num_edges;
    let m_rows = 3 * problem.topology.free_node_indices.len();

    let u = compute_target_member_vectors(problem, target_free_xyz);
    let (m_mat, p) = build_equilibrium_matrix(problem, &u);

    let (k_mat, rhs) = build_augmented_system(&m_mat, &p, regularization);

    let ldl = Ldl::new()
        .fill_in_reduction(sprs::FillInReduction::ReverseCuthillMcKee)
        .check_symmetry(sprs::SymmetryCheck::DontCheckSymmetry)
        .numeric(k_mat.view())
        .map_err(|e| TheseusError::Linalg(e))?;

    let sol = ldl.solve(&rhs);

    // q occupies the last ne entries of the solution vector
    let q: Vec<f64> = sol[m_rows..m_rows + ne].to_vec();

    for (i, &v) in q.iter().enumerate() {
        if !v.is_finite() {
            return Err(TheseusError::Solver(format!(
                "pseudoinverse (augmented) produced non-finite q at edge {i}; \
                 increase regularisation or check target geometry",
            )));
        }
    }

    if q.len() != ne {
        return Err(TheseusError::Shape(format!(
            "pseudoinverse (augmented) returned {} q values, expected {ne}",
            q.len()
        )));
    }

    Ok(q)
}

// ─────────────────────────────────────────────────────────────
//  Pseudoinverse L1  (IRLS — iteratively reweighted least squares)
// ─────────────────────────────────────────────────────────────

/// Find force densities via L1-minimisation of the equilibrium residual.
///
/// Minimises `‖Mq − p‖₁` (sum of absolute residuals) using IRLS:
/// each iteration solves a weighted least-squares problem
/// `(M^T W M + λI) q = M^T W p` where `W = diag(1/max(|r_i|, ε))`.
///
/// Warm-starts from the L2 pseudoinverse solution for fast convergence.
pub fn solve_pseudoinverse_l1(
    problem: &Problem,
    target_free_xyz: &Array2<f64>,
    regularization: f64,
    max_iter: usize,
) -> Result<Vec<f64>, TheseusError> {
    let ne = problem.topology.num_edges;

    let u = compute_target_member_vectors(problem, target_free_xyz);
    let (m_mat, p) = build_equilibrium_matrix(problem, &u);

    let m_t = m_mat.transpose_view().to_csc();
    let m_csc = m_mat.to_csc();
    let m_rows = p.len();

    // Warm-start: L2 solution  (M^T M + λI) q = M^T p
    let mut g_l2: CsMat<f64> = (&m_t * &m_csc).to_csc();
    if regularization > 0.0 {
        add_diagonal(&mut g_l2, regularization);
    }
    let h_l2 = sparse_matvec(&m_t, &p);
    let ldl_l2 = Ldl::new()
        .fill_in_reduction(sprs::FillInReduction::ReverseCuthillMcKee)
        .check_symmetry(sprs::SymmetryCheck::DontCheckSymmetry)
        .numeric(g_l2.view())
        .map_err(|e| TheseusError::Linalg(e))?;
    let mut q = ldl_l2.solve(&h_l2);

    const ABS_EPS: f64 = 1e-12;
    let mut prev_l1 = f64::MAX;

    for _ in 0..max_iter {
        // r = M*q − p
        let mut r = sparse_matvec(&m_csc, &q);
        for (ri, &pi) in r.iter_mut().zip(p.iter()) {
            *ri -= pi;
        }

        // Adaptive epsilon: proportional to the largest residual so that
        // the weight ratio stays bounded (~1e4:1), preventing the weighted
        // Gram matrix from becoming ill-conditioned.
        let max_abs_r = r.iter().map(|ri| ri.abs()).fold(0.0_f64, f64::max);
        let eps_iter = (1e-4 * max_abs_r).max(ABS_EPS);

        // L1 objective and IRLS weights  (w_i = 1/max(|r_i|, eps))
        let mut l1_obj = 0.0;
        let mut sqrt_w = vec![0.0; m_rows];
        let mut wp = vec![0.0; m_rows];
        for i in 0..m_rows {
            let abs_r = r[i].abs();
            l1_obj += abs_r;
            let w_i = 1.0 / abs_r.max(eps_iter);
            sqrt_w[i] = w_i.sqrt();
            wp[i] = w_i * p[i];
        }

        // Normalize weights so max(w) = 1, keeping M^T W M entries O(1).
        // Divide regularisation by the same factor to preserve the solution.
        let w_max = sqrt_w.iter().map(|s| s * s).fold(0.0_f64, f64::max);
        let effective_reg = if w_max > 0.0 {
            let inv_sqrt_wmax = 1.0 / w_max.sqrt();
            for i in 0..m_rows {
                sqrt_w[i] *= inv_sqrt_wmax;
                wp[i] /= w_max;
            }
            regularization / w_max
        } else {
            regularization
        };

        // Convergence: relative change in L1 objective
        if prev_l1 < f64::MAX {
            let rel_change = (prev_l1 - l1_obj).abs() / (prev_l1 + ABS_EPS);
            if rel_change < 1e-8 {
                break;
            }
        }
        prev_l1 = l1_obj;

        // Weighted normal equations: G = M_w^T M_w + λ'I,  h = M^T W' p
        let m_w = row_scaled_copy(&m_csc, &sqrt_w);
        let m_w_t = m_w.transpose_view().to_csc();
        let mut g: CsMat<f64> = (&m_w_t * &m_w).to_csc();
        if effective_reg > 0.0 {
            add_diagonal(&mut g, effective_reg);
        }
        let h = sparse_matvec(&m_t, &wp);

        let ldl = Ldl::new()
            .fill_in_reduction(sprs::FillInReduction::ReverseCuthillMcKee)
            .check_symmetry(sprs::SymmetryCheck::DontCheckSymmetry)
            .numeric(g.view())
            .map_err(|e| TheseusError::Linalg(e))?;

        q = ldl.solve(&h);
    }

    // Validate solution
    for (i, &v) in q.iter().enumerate() {
        if !v.is_finite() {
            return Err(TheseusError::Solver(format!(
                "pseudoinverse L1 produced non-finite q at edge {i}; \
                 increase regularisation or check target geometry",
            )));
        }
    }
    if q.len() != ne {
        return Err(TheseusError::Shape(format!(
            "pseudoinverse L1 returned {} q values, expected {ne}",
            q.len()
        )));
    }

    Ok(q)
}

// ─────────────────────────────────────────────────────────────
//  Pseudoinverse L1  (augmented saddle-point IRLS)
// ─────────────────────────────────────────────────────────────

/// Find force densities via L1-minimisation using the augmented saddle-point system.
///
/// Equivalent to `solve_pseudoinverse_l1` but each IRLS iteration factorises
/// the augmented system instead of forming M_w^T M_w.  Avoids fill-in explosion
/// at large scales.  Requires `λ > 0`.
pub fn solve_pseudoinverse_l1_augmented(
    problem: &Problem,
    target_free_xyz: &Array2<f64>,
    regularization: f64,
    max_iter: usize,
) -> Result<Vec<f64>, TheseusError> {
    if regularization <= 0.0 {
        return Err(TheseusError::Solver(
            "augmented system requires regularization > 0".into(),
        ));
    }

    let ne = problem.topology.num_edges;

    let u = compute_target_member_vectors(problem, target_free_xyz);
    let (m_mat, p) = build_equilibrium_matrix(problem, &u);

    let m_csc = m_mat.to_csc();
    let m_rows = p.len();

    // Warm-start: L2 solution via augmented system
    let (k_l2, rhs_l2) = build_augmented_system(&m_mat, &p, regularization);
    let ldl_l2 = Ldl::new()
        .fill_in_reduction(sprs::FillInReduction::ReverseCuthillMcKee)
        .check_symmetry(sprs::SymmetryCheck::DontCheckSymmetry)
        .numeric(k_l2.view())
        .map_err(|e| TheseusError::Linalg(e))?;
    let sol_l2 = ldl_l2.solve(&rhs_l2);
    let mut q: Vec<f64> = sol_l2[m_rows..m_rows + ne].to_vec();

    const ABS_EPS: f64 = 1e-12;
    let mut prev_l1 = f64::MAX;

    for _ in 0..max_iter {
        // r = M*q − p
        let mut r = sparse_matvec(&m_csc, &q);
        for (ri, &pi) in r.iter_mut().zip(p.iter()) {
            *ri -= pi;
        }

        let max_abs_r = r.iter().map(|ri| ri.abs()).fold(0.0_f64, f64::max);
        let eps_iter = (1e-4 * max_abs_r).max(ABS_EPS);

        // L1 objective and IRLS weights
        let mut l1_obj = 0.0;
        let mut sqrt_w = vec![0.0; m_rows];
        for i in 0..m_rows {
            let abs_r = r[i].abs();
            l1_obj += abs_r;
            let w_i = 1.0 / abs_r.max(eps_iter);
            sqrt_w[i] = w_i.sqrt();
        }

        // Normalize weights so max(w) = 1
        let w_max = sqrt_w.iter().map(|s| s * s).fold(0.0_f64, f64::max);
        let effective_reg = if w_max > 0.0 {
            let inv_sqrt_wmax = 1.0 / w_max.sqrt();
            for sw in sqrt_w.iter_mut() {
                *sw *= inv_sqrt_wmax;
            }
            regularization / w_max
        } else {
            regularization
        };

        // Convergence: relative change in L1 objective
        if prev_l1 < f64::MAX {
            let rel_change = (prev_l1 - l1_obj).abs() / (prev_l1 + ABS_EPS);
            if rel_change < 1e-8 {
                break;
            }
        }
        prev_l1 = l1_obj;

        // Build M_w = diag(sqrt_w) * M, then augmented system with
        // RHS top = sqrt_w ⊙ p  (so that M_w^T * rhs_top = M^T W p).
        let m_w = row_scaled_copy(&m_csc, &sqrt_w);
        let rhs_top: Vec<f64> = (0..m_rows).map(|i| sqrt_w[i] * p[i]).collect();

        let (k_mat, rhs) = build_augmented_system(&m_w, &rhs_top, effective_reg);

        let ldl = Ldl::new()
            .fill_in_reduction(sprs::FillInReduction::ReverseCuthillMcKee)
            .check_symmetry(sprs::SymmetryCheck::DontCheckSymmetry)
            .numeric(k_mat.view())
            .map_err(|e| TheseusError::Linalg(e))?;

        let sol = ldl.solve(&rhs);
        q = sol[m_rows..m_rows + ne].to_vec();
    }

    for (i, &v) in q.iter().enumerate() {
        if !v.is_finite() {
            return Err(TheseusError::Solver(format!(
                "pseudoinverse L1 (augmented) produced non-finite q at edge {i}; \
                 increase regularisation or check target geometry",
            )));
        }
    }
    if q.len() != ne {
        return Err(TheseusError::Shape(format!(
            "pseudoinverse L1 (augmented) returned {} q values, expected {ne}",
            q.len()
        )));
    }

    Ok(q)
}

// ─────────────────────────────────────────────────────────────
//  Pseudoinverse dispatcher
// ─────────────────────────────────────────────────────────────

/// Dispatch to L2 or L1 pseudoinverse, with normal-equations or augmented system.
///
/// The augmented path avoids forming M^T M and is faster for large meshes
/// (>50k edges).  It requires `regularization > 0`.
pub fn solve_pseudoinverse_dispatch(
    problem: &Problem,
    target_free_xyz: &Array2<f64>,
    regularization: f64,
    use_l2: bool,
    max_l1_iter: usize,
    use_augmented: bool,
) -> Result<Vec<f64>, TheseusError> {
    match (use_l2, use_augmented) {
        (true, false) => solve_pseudoinverse(problem, target_free_xyz, regularization),
        (true, true) => solve_pseudoinverse_augmented(problem, target_free_xyz, regularization),
        (false, false) => solve_pseudoinverse_l1(problem, target_free_xyz, regularization, max_l1_iter),
        (false, true) => solve_pseudoinverse_l1_augmented(problem, target_free_xyz, regularization, max_l1_iter),
    }
}

// ─────────────────────────────────────────────────────────────
//  NNLS  (spectral projected gradient)
// ─────────────────────────────────────────────────────────────

/// Find non-negative force densities via spectral projected gradient.
///
/// Solves  `min ‖Mq − p‖²   s.t.  q ≥ 0`  using only sparse matvecs.
pub fn solve_nnls(
    problem: &Problem,
    target_free_xyz: &Array2<f64>,
    max_iter: usize,
    tol: f64,
) -> Result<Vec<f64>, TheseusError> {
    let ne = problem.topology.num_edges;
    let nn_free = problem.topology.free_node_indices.len();

    let u = compute_target_member_vectors(problem, target_free_xyz);

    let cn_t_csc = problem.topology.free_incidence.transpose_view().to_csc();
    let cn_csc = problem.topology.free_incidence.to_csc();

    // p = [Pn_x; Pn_y; Pn_z]
    let loads = &problem.free_node_loads;
    let mut p = vec![0.0; 3 * nn_free];
    for d in 0..3 {
        let off = d * nn_free;
        for i in 0..nn_free {
            p[off + i] = loads[[i, d]];
        }
    }

    // Initialise q = 1.0 (feasible starting point)
    let mut q = vec![1.0; ne];

    let mut prev_g = vec![0.0; ne];
    let mut prev_q = vec![0.0; ne];
    let mut alpha = 1.0;

    for iter in 0..max_iter {
        // r = M·q − p
        let mut r = apply_m(&cn_t_csc, &u, &q, nn_free);
        for (ri, &pi) in r.iter_mut().zip(p.iter()) {
            *ri -= pi;
        }

        // g = M^T · r  (gradient of ½‖Mq−p‖²)
        let g = apply_mt(&cn_csc, &u, &r, nn_free);

        // Barzilai-Borwein step size (after first iteration)
        if iter > 0 {
            let mut dq_dot_dg = 0.0;
            let mut dg_dot_dg = 0.0;
            for k in 0..ne {
                let dq = q[k] - prev_q[k];
                let dg = g[k] - prev_g[k];
                dq_dot_dg += dq * dg;
                dg_dot_dg += dg * dg;
            }
            if dg_dot_dg > 0.0 && dq_dot_dg > 0.0 {
                alpha = dq_dot_dg / dg_dot_dg;
            }
        }

        // Save for BB computation
        prev_q.copy_from_slice(&q);
        prev_g.copy_from_slice(&g);

        // Projected gradient step: q = max(0, q − α·g)
        let mut proj_grad_norm_sq = 0.0;
        for k in 0..ne {
            let new_q = (q[k] - alpha * g[k]).max(0.0);
            let pg = q[k] - new_q;
            proj_grad_norm_sq += pg * pg;
            q[k] = new_q;
        }

        // Convergence check: ‖projected gradient‖ < tol
        if proj_grad_norm_sq.sqrt() < tol {
            break;
        }
    }

    // Validate solution
    for (i, &v) in q.iter().enumerate() {
        if !v.is_finite() {
            return Err(TheseusError::Solver(format!(
                "NNLS produced non-finite q at edge {i}",
            )));
        }
    }

    Ok(q)
}

// ─────────────────────────────────────────────────────────────
//  Sparse helpers
// ─────────────────────────────────────────────────────────────

/// Create a copy of a CSC matrix with each row `i` scaled by `row_scales[i]`.
fn row_scaled_copy(mat: &CsMat<f64>, row_scales: &[f64]) -> CsMat<f64> {
    let mut scaled = mat.clone();
    let indptr = scaled.indptr().to_owned();
    let raw = indptr.raw_storage();
    for col in 0..scaled.cols() {
        for nz in raw[col]..raw[col + 1] {
            let row = scaled.indices()[nz];
            scaled.data_mut()[nz] *= row_scales[row];
        }
    }
    scaled
}

/// Add a scalar to every diagonal entry of a CSC matrix (in-place).
fn add_diagonal(mat: &mut CsMat<f64>, value: f64) {
    let n = mat.cols().min(mat.rows());
    let indptr = mat.indptr().to_owned();
    let raw = indptr.raw_storage();
    for col in 0..n {
        let start = raw[col];
        let end_ = raw[col + 1];
        for nz in start..end_ {
            if mat.indices()[nz] == col {
                mat.data_mut()[nz] += value;
                break;
            }
        }
    }
}

/// Sparse matrix × dense vector  (CSC layout).
fn sparse_matvec(mat: &CsMat<f64>, x: &[f64]) -> Vec<f64> {
    let nrows = mat.rows();
    let ncols = mat.cols();
    let mut y = vec![0.0; nrows];
    let indptr = mat.indptr();
    let raw = indptr.raw_storage();
    for col in 0..ncols {
        let xc = x[col];
        if xc == 0.0 {
            continue;
        }
        let start = raw[col];
        let end_ = raw[col + 1];
        for nz in start..end_ {
            let row = mat.indices()[nz];
            let val = mat.data()[nz];
            y[row] += val * xc;
        }
    }
    y
}
