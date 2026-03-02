//! Inverse FDM solvers: find force densities q from a target geometry.
//!
//! Given target free-node positions, the FDM equilibrium per dimension d is:
//!
//!   Cn^T  diag(u_d)  q  =  Pn_d
//!
//! where u_d = C · N_target are the target member coordinate differences.
//! Stacking all 3 dimensions gives the sparse system  M q = p.
//!
//! Two solvers are provided:
//!
//!   - **Pseudoinverse** (Tikhonov-regularised normal equations):
//!       (M^T M + λI) q = M^T p    →  sparse LDL via `sprs-ldl`
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
