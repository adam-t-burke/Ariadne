//! Inverse FDM solvers: find force densities q from a target geometry.
//!
//! Given target free-node positions, the FDM equilibrium per dimension d is:
//!
//!   Cn^T  diag(u_d)  q  =  Pn_d
//!
//! where u_d = C · N_target are the target member coordinate differences.
//! Stacking all 3 dimensions gives the sparse system  M q = p.
//!
//! When `solve_for_q` is false, the system is reformulated to solve for axial
//! forces F directly by normalising u to unit vectors:
//!
//!   Cn^T  diag(e_d)  F  =  Pn_d       where  e_d = u_d / ‖u‖
//!
//! Force densities are then recovered as  q = F / L.
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

use crate::sparse::SparseColMatOwned;
use crate::types::{Factorization, FactorizationStrategy, Problem, TheseusError};
use ndarray::Array2;

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
    let inc = &topo.incidence;
    let mut u = Array2::<f64>::zeros((ne, 3));
    for col in 0..nn {
        let start = inc.col_ptrs[col] as usize;
        let end_ = inc.col_ptrs[col + 1] as usize;
        for nz in start..end_ {
            let row = inc.row_indices[nz] as usize;
            let val = inc.values[nz];
            for d in 0..3 {
                u[[row, d]] += val * nf[[col, d]];
            }
        }
    }
    u
}

// ─────────────────────────────────────────────────────────────
//  Normalise u → unit vectors (for solve-for-F mode)
// ─────────────────────────────────────────────────────────────

const ZERO_LENGTH_TOL: f64 = 1e-14;

/// Normalise member vectors in-place to unit vectors and return the edge lengths.
///
/// Each row `u[k,:]` is divided by its Euclidean norm `L_k`.  Returns `L`
/// (length `ne`).  Errors if any edge has near-zero length.
fn normalise_member_vectors(u: &mut Array2<f64>) -> Result<Vec<f64>, TheseusError> {
    let ne = u.nrows();
    let mut lengths = vec![0.0; ne];
    for k in 0..ne {
        let lsq = u[[k, 0]] * u[[k, 0]] + u[[k, 1]] * u[[k, 1]] + u[[k, 2]] * u[[k, 2]];
        let l = lsq.sqrt();
        if l < ZERO_LENGTH_TOL {
            return Err(TheseusError::Solver(format!(
                "target edge {k} has near-zero length ({l:.2e}); \
                 cannot solve for F with degenerate edges",
            )));
        }
        let inv_l = 1.0 / l;
        u[[k, 0]] *= inv_l;
        u[[k, 1]] *= inv_l;
        u[[k, 2]] *= inv_l;
        lengths[k] = l;
    }
    Ok(lengths)
}

/// Convert axial forces F back to force densities: `q[k] = F[k] / L[k]`.
fn forces_to_q(f: &[f64], lengths: &[f64]) -> Vec<f64> {
    f.iter()
        .zip(lengths.iter())
        .map(|(&fi, &li)| fi / li)
        .collect()
}

// ─────────────────────────────────────────────────────────────
//  Build sparse equilibrium matrix  M  (3·nn_free × ne)
// ─────────────────────────────────────────────────────────────

/// Build the sparse matrix M by vertically stacking Cn^T · diag(u_d) for d=x,y,z.
///
/// When `enforce_zero_rx` is true and there are fixed nodes, appends Cf^T · diag(u_x)
/// and zeros to the RHS to strictly enforce R_x = 0 at supports.
///
/// M has the same column-sparsity as Cn^T (at most 2 non-zeros per column per
/// dimension block → ≤ 6 per column total).  Uses COO → CSC for construction.
fn build_equilibrium_matrix(
    problem: &Problem,
    u: &Array2<f64>,
    enforce_zero_rx: bool,
) -> (SparseColMatOwned, Vec<f64>) {
    let topo = &problem.topology;
    let ne = topo.num_edges;
    let nn_free = topo.free_node_indices.len();
    let nn_fixed = topo.fixed_node_indices.len();
    let rx_rows = if enforce_zero_rx && nn_fixed > 0 {
        nn_fixed
    } else {
        0
    };
    let m_rows = 3 * nn_free + rx_rows;

    // Cn is (ne × nn_free).  Cn^T is (nn_free × ne).
    let cn_t = topo.free_incidence.transpose();

    let mut triplets: Vec<(u32, u32, f64)> = Vec::new();
    for d in 0..3usize {
        let row_offset = d * nn_free;
        for col in 0..ne {
            let scale = u[[col, d]];
            if scale == 0.0 {
                continue;
            }
            let start = cn_t.col_ptrs[col] as usize;
            let end_ = cn_t.col_ptrs[col + 1] as usize;
            for nz in start..end_ {
                let row = cn_t.row_indices[nz] as usize;
                let val = cn_t.values[nz];
                triplets.push(((row_offset + row) as u32, col as u32, val * scale));
            }
        }
    }

    // Cf^T · diag(u_x): append rows for R_x = 0 at fixed nodes
    if rx_rows > 0 {
        let cf = &topo.fixed_incidence; // ne × nn_fixed
        let row_offset = 3 * nn_free;
        for j in 0..nn_fixed {
            let start = cf.col_ptrs[j] as usize;
            let end_ = cf.col_ptrs[j + 1] as usize;
            for nz in start..end_ {
                let i = cf.row_indices[nz] as usize;
                let val = cf.values[nz];
                let scale = u[[i, 0]];
                if scale != 0.0 {
                    triplets.push(((row_offset + j) as u32, i as u32, val * scale));
                }
            }
        }
    }

    let m_mat = SparseColMatOwned::from_triplets(m_rows, ne, &triplets)
        .expect("build_equilibrium_matrix");

    // p = [Pn_x; Pn_y; Pn_z; 0; ...; 0]  (zeros for Rx rows when enforce_zero_rx)
    let loads = &problem.free_node_loads;
    let mut p = vec![0.0; m_rows];
    for d in 0..3 {
        let off = d * nn_free;
        for i in 0..nn_free {
            p[off + i] = loads[[i, d]];
        }
    }
    // rx_rows zeros are already in p via vec![0.0; m_rows]

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
    m_mat: &SparseColMatOwned,
    rhs_top: &[f64],
    regularization: f64,
) -> (SparseColMatOwned, Vec<f64>) {
    let m = m_mat.nrows;
    let n = m_mat.ncols;
    let total = m + n;

    let mut triplets: Vec<(u32, u32, f64)> = Vec::new();

    // Top-left: I  (m × m)
    for i in 0..m {
        triplets.push((i as u32, i as u32, 1.0));
    }

    // Upper-right: M  and  lower-left: M^T  (symmetric pair)
    for col in 0..n {
        let start = m_mat.col_ptrs[col] as usize;
        let end_ = m_mat.col_ptrs[col + 1] as usize;
        for nz in start..end_ {
            let row = m_mat.row_indices[nz] as usize;
            let val = m_mat.values[nz];
            triplets.push((row as u32, (m + col) as u32, val));
            triplets.push(((m + col) as u32, row as u32, val));
        }
    }

    // Bottom-right: -λI  (n × n)
    for j in 0..n {
        triplets.push(((m + j) as u32, (m + j) as u32, -regularization));
    }

    let k_mat = SparseColMatOwned::from_triplets(total, total, &triplets)
        .expect("build_augmented_system");

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
fn apply_m(cn_t: &SparseColMatOwned, u: &Array2<f64>, q: &[f64], nn_free: usize) -> Vec<f64> {
    let ne = q.len();
    let mut r = vec![0.0; 3 * nn_free];

    for d in 0..3usize {
        let off = d * nn_free;
        for col in 0..ne {
            let w = u[[col, d]] * q[col];
            if w == 0.0 {
                continue;
            }
            let start = cn_t.col_ptrs[col] as usize;
            let end_ = cn_t.col_ptrs[col + 1] as usize;
            for nz in start..end_ {
                let row = cn_t.row_indices[nz] as usize;
                let val = cn_t.values[nz];
                r[off + row] += val * w;
            }
        }
    }
    r
}

/// Compute g = M^T · r  without forming M.
///
/// g_k = Σ_d  u[k,d] · (Cn · r_d)[k]
fn apply_mt(cn: &SparseColMatOwned, u: &Array2<f64>, r: &[f64], nn_free: usize) -> Vec<f64> {
    let ne = u.nrows();
    let mut g = vec![0.0; ne];

    for d in 0..3usize {
        let off = d * nn_free;
        let r_d = &r[off..off + nn_free];

        for col in 0..nn_free {
            let rd_val = r_d[col];
            if rd_val == 0.0 {
                continue;
            }
            let start = cn.col_ptrs[col] as usize;
            let end_ = cn.col_ptrs[col + 1] as usize;
            for nz in start..end_ {
                let row = cn.row_indices[nz] as usize;
                let val = cn.values[nz];
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
/// When `solve_for_q` is false, solves for axial forces F instead and
/// recovers q = F / L from the target edge lengths.
pub fn solve_pseudoinverse(
    problem: &Problem,
    target_free_xyz: &Array2<f64>,
    regularization: f64,
    enforce_zero_rx: bool,
    solve_for_q: bool,
) -> Result<Vec<f64>, TheseusError> {
    let ne = problem.topology.num_edges;

    let mut u = compute_target_member_vectors(problem, target_free_xyz);
    let target_lengths = if solve_for_q {
        None
    } else {
        Some(normalise_member_vectors(&mut u)?)
    };
    let (m_mat, p) = build_equilibrium_matrix(problem, &u, enforce_zero_rx);

    // G = M^T M  (ne × ne, sparse)
    let m_t = m_mat.transpose();
    let mut g = SparseColMatOwned::sparse_times_sparse(&m_t, &m_mat)
        .map_err(|e| TheseusError::Solver(e))?;

    // Add regularisation: G += λ I
    if regularization > 0.0 {
        g.add_diagonal(regularization);
    }

    // h = M^T p  (sparse-dense matvec)
    let h = m_t.matvec(&p);

    // Factorise G and solve (LDL for normal equations)
    let fac = Factorization::new(&g, FactorizationStrategy::LDL)?;
    let sol = fac.solve(&h);

    let q = match target_lengths {
        Some(ref lengths) => forces_to_q(&sol, lengths),
        None => sol,
    };

    // Validate solution
    for (i, &v) in q.iter().enumerate() {
        if !v.is_finite() {
            return Err(TheseusError::Solver(format!(
                "pseudoinverse produced non-finite q at edge {i}; \
                 increase regularisation or check target geometry",
            )));
        }
    }

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
/// When `solve_for_q` is false, solves for axial forces F instead and
/// recovers q = F / L from the target edge lengths.
pub fn solve_pseudoinverse_augmented(
    problem: &Problem,
    target_free_xyz: &Array2<f64>,
    regularization: f64,
    enforce_zero_rx: bool,
    solve_for_q: bool,
) -> Result<Vec<f64>, TheseusError> {
    if regularization <= 0.0 {
        return Err(TheseusError::Solver(
            "augmented system requires regularization > 0".into(),
        ));
    }

    let ne = problem.topology.num_edges;

    let mut u = compute_target_member_vectors(problem, target_free_xyz);
    let target_lengths = if solve_for_q {
        None
    } else {
        Some(normalise_member_vectors(&mut u)?)
    };
    let (m_mat, p) = build_equilibrium_matrix(problem, &u, enforce_zero_rx);
    let m_rows = p.len();

    let (k_mat, rhs) = build_augmented_system(&m_mat, &p, regularization);

    let fac = Factorization::new(&k_mat, FactorizationStrategy::LDL)?;
    let sol = fac.solve(&rhs);

    let raw: Vec<f64> = sol[m_rows..m_rows + ne].to_vec();
    let q = match target_lengths {
        Some(ref lengths) => forces_to_q(&raw, lengths),
        None => raw,
    };

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
/// When `solve_for_q` is false, solves for axial forces F instead and
/// recovers q = F / L from the target edge lengths.
pub fn solve_pseudoinverse_l1(
    problem: &Problem,
    target_free_xyz: &Array2<f64>,
    regularization: f64,
    max_iter: usize,
    enforce_zero_rx: bool,
    solve_for_q: bool,
) -> Result<Vec<f64>, TheseusError> {
    let ne = problem.topology.num_edges;

    let mut u = compute_target_member_vectors(problem, target_free_xyz);
    let target_lengths = if solve_for_q {
        None
    } else {
        Some(normalise_member_vectors(&mut u)?)
    };
    let (m_mat, p) = build_equilibrium_matrix(problem, &u, enforce_zero_rx);

    let m_t = m_mat.transpose();
    let m_rows = p.len();

    // Warm-start: L2 solution  (M^T M + λI) q = M^T p
    let mut g_l2 = SparseColMatOwned::sparse_times_sparse(&m_t, &m_mat)
        .map_err(|e| TheseusError::Solver(e))?;
    if regularization > 0.0 {
        g_l2.add_diagonal(regularization);
    }
    let h_l2 = m_t.matvec(&p);
    let fac_l2 = Factorization::new(&g_l2, FactorizationStrategy::LDL)?;
    let mut sol = fac_l2.solve(&h_l2);

    const ABS_EPS: f64 = 1e-12;
    let mut prev_l1 = f64::MAX;

    for _ in 0..max_iter {
        // r = M*sol − p
        let mut r = m_mat.matvec(&sol);
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
        let m_w = row_scaled_copy(&m_mat, &sqrt_w);
        let m_w_t = m_w.transpose();
        let mut g = SparseColMatOwned::sparse_times_sparse(&m_w_t, &m_w)
            .map_err(|e| TheseusError::Solver(e))?;
        if effective_reg > 0.0 {
            g.add_diagonal(effective_reg);
        }
        let h = m_t.matvec(&wp);

        let fac = Factorization::new(&g, FactorizationStrategy::LDL)?;
        sol = fac.solve(&h);
    }

    let q = match target_lengths {
        Some(ref lengths) => forces_to_q(&sol, lengths),
        None => sol,
    };

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
/// When `solve_for_q` is false, solves for axial forces F instead and
/// recovers q = F / L from the target edge lengths.
pub fn solve_pseudoinverse_l1_augmented(
    problem: &Problem,
    target_free_xyz: &Array2<f64>,
    regularization: f64,
    max_iter: usize,
    enforce_zero_rx: bool,
    solve_for_q: bool,
) -> Result<Vec<f64>, TheseusError> {
    if regularization <= 0.0 {
        return Err(TheseusError::Solver(
            "augmented system requires regularization > 0".into(),
        ));
    }

    let ne = problem.topology.num_edges;

    let mut u = compute_target_member_vectors(problem, target_free_xyz);
    let target_lengths = if solve_for_q {
        None
    } else {
        Some(normalise_member_vectors(&mut u)?)
    };
    let (m_mat, p) = build_equilibrium_matrix(problem, &u, enforce_zero_rx);

    let m_rows = p.len();

    // Warm-start: L2 solution via augmented system
    let (k_l2, rhs_l2) = build_augmented_system(&m_mat, &p, regularization);
    let fac_l2 = Factorization::new(&k_l2, FactorizationStrategy::LDL)?;
    let sol_l2 = fac_l2.solve(&rhs_l2);
    let mut sol: Vec<f64> = sol_l2[m_rows..m_rows + ne].to_vec();

    const ABS_EPS: f64 = 1e-12;
    let mut prev_l1 = f64::MAX;

    for _ in 0..max_iter {
        // r = M*sol − p
        let mut r = m_mat.matvec(&sol);
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
        let m_w = row_scaled_copy(&m_mat, &sqrt_w);
        let rhs_top: Vec<f64> = (0..m_rows).map(|i| sqrt_w[i] * p[i]).collect();

        let (k_mat, rhs) = build_augmented_system(&m_w, &rhs_top, effective_reg);

        let fac = Factorization::new(&k_mat, FactorizationStrategy::LDL)?;
        let iter_sol = fac.solve(&rhs);
        sol = iter_sol[m_rows..m_rows + ne].to_vec();
    }

    let q = match target_lengths {
        Some(ref lengths) => forces_to_q(&sol, lengths),
        None => sol,
    };

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
///
/// When `enforce_zero_rx` is true, augments the system to strictly enforce
/// R_x = 0 at fixed supports (channels horizontal thrust into tie members).
///
/// When `solve_for_q` is false, solves for axial forces F directly (using
/// unit direction vectors instead of coordinate differences) and recovers
/// q = F / L from the target edge lengths.
pub fn solve_pseudoinverse_dispatch(
    problem: &Problem,
    target_free_xyz: &Array2<f64>,
    regularization: f64,
    use_l2: bool,
    max_l1_iter: usize,
    use_augmented: bool,
    enforce_zero_rx: bool,
    solve_for_q: bool,
) -> Result<Vec<f64>, TheseusError> {
    match (use_l2, use_augmented) {
        (true, false) => solve_pseudoinverse(problem, target_free_xyz, regularization, enforce_zero_rx, solve_for_q),
        (true, true) => solve_pseudoinverse_augmented(problem, target_free_xyz, regularization, enforce_zero_rx, solve_for_q),
        (false, false) => solve_pseudoinverse_l1(problem, target_free_xyz, regularization, max_l1_iter, enforce_zero_rx, solve_for_q),
        (false, true) => solve_pseudoinverse_l1_augmented(problem, target_free_xyz, regularization, max_l1_iter, enforce_zero_rx, solve_for_q),
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

    let cn_t = problem.topology.free_incidence.transpose();
    let cn = &problem.topology.free_incidence;

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
        let mut r = apply_m(&cn_t, &u, &q, nn_free);
        for (ri, &pi) in r.iter_mut().zip(p.iter()) {
            *ri -= pi;
        }

        // g = M^T · r  (gradient of ½‖Mq−p‖²)
        let g = apply_mt(cn, &u, &r, nn_free);

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
fn row_scaled_copy(mat: &SparseColMatOwned, row_scales: &[f64]) -> SparseColMatOwned {
    let mut scaled = mat.clone();
    for col in 0..scaled.ncols {
        let start = scaled.col_ptrs[col] as usize;
        let end_ = scaled.col_ptrs[col + 1] as usize;
        for nz in start..end_ {
            let row = scaled.row_indices[nz] as usize;
            scaled.values[nz] *= row_scales[row];
        }
    }
    scaled
}
