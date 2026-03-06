//! Forward FDM solver: assemble A(q), build RHS, factorise, triangular solve.
//!
//! Mirrors `src/FDM.jl` from the Julia code.

use crate::types::{FdmCache, Factorization, FactorizationStrategy, Problem, TheseusError};
use ndarray::Array2;
use rayon::prelude::*;

// ─────────────────────────────────────────────────────────────
//  Fixed-node position assembly
// ─────────────────────────────────────────────────────────────

/// Write current fixed-node positions into `cache.nf`, overlaying reference
/// positions and (optionally) variable anchor positions.
pub fn update_fixed_positions(cache: &mut FdmCache, problem: &Problem, anchor_positions: &Array2<f64>) {
    let fixed = &problem.topology.fixed_node_indices;
    let nn_fixed = fixed.len();

    // 1. Copy reference positions for fixed nodes
    let ref_pos = &problem.anchors.reference_positions;
    if ref_pos.nrows() == nn_fixed {
        for (i, &node) in fixed.iter().enumerate() {
            for d in 0..3 {
                cache.nf[[node, d]] = ref_pos[[i, d]];
            }
        }
    } else if ref_pos.nrows() == problem.topology.num_nodes {
        for &node in fixed {
            for d in 0..3 {
                cache.nf[[node, d]] = ref_pos[[node, d]];
            }
        }
    }

    // 2. Overlay variable anchors
    for (i, &node) in problem.anchors.variable_indices.iter().enumerate() {
        for d in 0..3 {
            cache.nf[[node, d]] = anchor_positions[[i, d]];
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  System matrix assembly  A(q)
// ─────────────────────────────────────────────────────────────

/// Zero-allocation in-place update of A's values from current q.
/// A = Cn^T diag(q) Cn  via the precomputed `q_to_nz` mapping.
pub fn assemble_a(cache: &mut FdmCache) {
    for v in cache.a_matrix.values.iter_mut() {
        *v = 0.0;
    }
    for (k, entries) in cache.q_to_nz.entries.iter().enumerate() {
        let qk = cache.q[k];
        for &(nz_idx, coeff) in entries {
            cache.a_matrix.values[nz_idx] += qk * coeff;
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  RHS assembly:  b = Pn − Cn^T diag(q) Cf Nf_fixed
// ─────────────────────────────────────────────────────────────

/// Build the right-hand side for A x = b.
///
/// Steps:
///   1. Copy fixed-node positions into dense buffer `nf_fixed`
///   2. cf_nf  = Cf * nf_fixed        (ne × 3)
///   3. q_cf_nf = diag(q) * cf_nf     (ne × 3)
///   4. rhs    = Pn − Cn^T * q_cf_nf  (nn_free × 3)
pub fn assemble_rhs(cache: &mut FdmCache, problem: &Problem) {
    let fixed = &problem.topology.fixed_node_indices;

    // 1. nf_fixed dense buffer
    for (i, &node) in fixed.iter().enumerate() {
        for d in 0..3 {
            cache.nf_fixed[[i, d]] = cache.nf[[node, d]];
        }
    }

    // 2. cf_nf = Cf * nf_fixed   (sparse × dense, column by column)
    spmm_into(&cache.cf, &cache.nf_fixed, &mut cache.cf_nf);

    // 3. q_cf_nf = diag(q) * cf_nf
    let ne = cache.q.len();
    for i in 0..ne {
        let qi = cache.q[i];
        for d in 0..3 {
            cache.q_cf_nf[[i, d]] = qi * cache.cf_nf[[i, d]];
        }
    }

    // 4. rhs = Pn − Cn^T * q_cf_nf
    cache.rhs.assign(&cache.pn);
    let cn_t = cache.cn.transpose();
    spmm_sub_into(&cn_t, &cache.q_cf_nf, &mut cache.rhs);
}

// ─────────────────────────────────────────────────────────────
//  Linear solve  (dense fallback — will be replaced with LDL)
// ─────────────────────────────────────────────────────────────

/// Factor A via sparse Cholesky or LDL^T and solve A x = rhs for all 3 columns.
///
/// On first call, performs a fresh factorization (symbolic + numeric).
/// On subsequent calls, reuses the symbolic structure via `Factorization::update()`
/// — only numeric values change.
///
/// If the preferred strategy (Cholesky) fails because the matrix is no longer
/// SPD (e.g. q values drifted negative during optimisation), we automatically
/// fall back to LDL and rebuild the factorization from scratch.
pub fn factor_and_solve(cache: &mut FdmCache, perturbation: f64) -> Result<(), TheseusError> {
    // Add diagonal perturbation if requested
    if perturbation > 0.0 {
        cache.a_matrix.add_diagonal(perturbation);
    }

    // Try the current factorization, falling back to LDL on Cholesky failure.
    let mut need_ldl_fallback = false;

    match &mut cache.factorization {
        Some(fac) => {
            if let Err(_e) = fac.update(&cache.a_matrix) {
                if fac.strategy() == FactorizationStrategy::Cholesky {
                    need_ldl_fallback = true;
                } else {
                    return Err(_e.into());
                }
            }
        }
        None => {
            match Factorization::new(&cache.a_matrix, cache.strategy) {
                Ok(fac) => {
                    cache.factorization = Some(fac);
                }
                Err(_e) if cache.strategy == FactorizationStrategy::Cholesky => {
                    need_ldl_fallback = true;
                }
                Err(e) => return Err(e.into()),
            }
        }
    }

    if need_ldl_fallback {
        cache.strategy = FactorizationStrategy::LDL;
        cache.factorization = None;
        cache.factorization = Some(Factorization::new(&cache.a_matrix, FactorizationStrategy::LDL)?);
    }

    // Batch solve for all 3 coordinate columns
    let n = cache.a_matrix.nrows;
    let mut rhs_flat = Vec::with_capacity(n * 3);
    for d in 0..3 {
        for i in 0..n {
            rhs_flat.push(cache.rhs[[i, d]]);
        }
    }
    let fac = cache.factorization.as_ref()
        .ok_or(TheseusError::MissingFactorization)?;
    let solutions = fac.solve_batch(&rhs_flat, 3);
    for (d, x) in solutions.into_iter().enumerate() {
        for i in 0..n {
            cache.x[[i, d]] = x[i];
        }
        if x.iter().any(|v| !v.is_finite()) {
            return Err(TheseusError::Solver(
                "FDM linear solve produced non-finite solution (singular or ill-conditioned equilibrium matrix). \
                 Check network connectivity, supports, and initial force densities.".into(),
            ));
        }
    }

    Ok(())
}

// ─────────────────────────────────────────────────────────────
//  Top-level forward solve
// ─────────────────────────────────────────────────────────────

/// Full forward FDM solve.  Updates `cache.x`, `cache.nf`,
/// `cache.member_lengths`, `cache.member_forces`, `cache.reactions`.
pub fn solve_fdm(
    cache: &mut FdmCache,
    q: &[f64],
    problem: &Problem,
    anchor_positions: &Array2<f64>,
    perturbation: f64,
) -> Result<(), TheseusError> {
    // 0. Sync q
    cache.q.copy_from_slice(q);

    // 1. Assemble A
    assemble_a(cache);

    // 2. Update fixed positions in Nf
    update_fixed_positions(cache, problem, anchor_positions);

    // 3. Assemble RHS
    assemble_rhs(cache, problem);

    // 4. Factor A and solve A x = rhs
    factor_and_solve(cache, perturbation)?;

    // 5. Write free-node positions back to Nf
    for (i, &node) in problem.topology.free_node_indices.iter().enumerate() {
        for d in 0..3 {
            cache.nf[[node, d]] = cache.x[[i, d]];
        }
    }

    // 6. Compute derived geometry
    compute_geometry(cache, problem);

    // 7. Ensure geometry is finite (catches NaN/Inf from positions or downstream)
    for (i, &len) in cache.member_lengths.iter().enumerate() {
        if !len.is_finite() {
            return Err(TheseusError::Solver(format!(
                "FDM geometry produced non-finite member length at edge {i} (node positions may be NaN/Inf; \
                 check network topology and initial q for singularity or ill-conditioning).",
            )));
        }
    }
    for (i, &f) in cache.member_forces.iter().enumerate() {
        if !f.is_finite() {
            return Err(TheseusError::Solver(format!(
                "FDM geometry produced non-finite member force at edge {i} (check network and initial q).",
            )));
        }
    }

    Ok(())
}

/// Compute member lengths, forces, and reactions from current positions and q.
/// Uses max(0, …) before sqrt to avoid NaN from floating-point negative squared length.
pub fn compute_geometry(cache: &mut FdmCache, problem: &Problem) {
    let ne = problem.topology.num_edges;

    // Per-edge member length and force (embarrassingly parallel — no write conflicts)
    let nf = &cache.nf;
    let edge_starts = &cache.edge_starts;
    let edge_ends = &cache.edge_ends;
    let q = &cache.q;

    cache.member_lengths.par_iter_mut().zip(cache.member_forces.par_iter_mut())
        .enumerate()
        .for_each(|(i, (len_out, force_out))| {
            let s = edge_starts[i];
            let e = edge_ends[i];

            let dx = nf[[e, 0]] - nf[[s, 0]];
            let dy = nf[[e, 1]] - nf[[s, 1]];
            let dz = nf[[e, 2]] - nf[[s, 2]];

            let len_sq = dx * dx + dy * dy + dz * dz;
            let len = len_sq.max(0.0).sqrt();
            *len_out = len;
            *force_out = q[i] * len;
        });

    // Reactions: fold/reduce per-thread buffers to avoid write conflicts
    let nn = cache.reactions.nrows();
    let reaction_sum = (0..ne).into_par_iter()
        .fold(
            || vec![0.0f64; nn * 3],
            |mut buf, i| {
                let s = edge_starts[i];
                let e = edge_ends[i];
                let qi = q[i];

                let rx = (nf[[e, 0]] - nf[[s, 0]]) * qi;
                let ry = (nf[[e, 1]] - nf[[s, 1]]) * qi;
                let rz = (nf[[e, 2]] - nf[[s, 2]]) * qi;

                buf[s * 3]     += rx;
                buf[s * 3 + 1] += ry;
                buf[s * 3 + 2] += rz;

                buf[e * 3]     -= rx;
                buf[e * 3 + 1] -= ry;
                buf[e * 3 + 2] -= rz;

                buf
            },
        )
        .reduce(
            || vec![0.0f64; nn * 3],
            |mut a, b| {
                for (ai, bi) in a.iter_mut().zip(b.iter()) {
                    *ai += *bi;
                }
                a
            },
        );

    for node in 0..nn {
        cache.reactions[[node, 0]] = reaction_sum[node * 3];
        cache.reactions[[node, 1]] = reaction_sum[node * 3 + 1];
        cache.reactions[[node, 2]] = reaction_sum[node * 3 + 2];
    }
}

// ─────────────────────────────────────────────────────────────
//  Sparse × dense helpers
// ─────────────────────────────────────────────────────────────

use crate::sparse::SparseColMatOwned;

/// out = A * B   where A is CSC (m × k), B is dense (k × 3), out is dense (m × 3).
fn spmm_into(a: &SparseColMatOwned, b: &Array2<f64>, out: &mut Array2<f64>) {
    out.fill(0.0);
    let ncols_a = a.ncols;
    for col in 0..ncols_a {
        let start = a.col_ptrs[col] as usize;
        let end_ = a.col_ptrs[col + 1] as usize;
        for nz in start..end_ {
            let row = a.row_indices[nz] as usize;
            let val = a.values[nz];
            for d in 0..3 {
                out[[row, d]] += val * b[[col, d]];
            }
        }
    }
}

/// out -= A * B   (subtract sparse-dense product from existing out).
fn spmm_sub_into(a: &SparseColMatOwned, b: &Array2<f64>, out: &mut Array2<f64>) {
    let ncols_a = a.ncols;
    for col in 0..ncols_a {
        let start = a.col_ptrs[col] as usize;
        let end_ = a.col_ptrs[col + 1] as usize;
        for nz in start..end_ {
            let row = a.row_indices[nz] as usize;
            let val = a.values[nz];
            for d in 0..3 {
                out[[row, d]] -= val * b[[col, d]];
            }
        }
    }
}
