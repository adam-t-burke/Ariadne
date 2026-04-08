//! Forward FDM solver: assemble A(q), build RHS, factorise, triangular solve.
//!
//! Mirrors `src/FDM.jl` from the Julia code.

use crate::types::{FdmCache, Factorization, FactorizationStrategy, Problem, SelfWeightParams, PressureParams, TheseusError};
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
//  Self-weight load computation
// ─────────────────────────────────────────────────────────────

/// Update `cache.sw_mu` from current forces (sizing mode only).
fn update_sizing_mu(cache: &mut FdmCache, rho: f64, sigma: f64) {
    let ne = cache.member_lengths.len();
    for k in 0..ne {
        let a_k = cache.member_forces[k].abs() / sigma;
        cache.sw_mu[k] = rho * a_k;
        cache.cross_section_areas[k] = a_k;
    }
}

/// Add lumped self-weight loads into `cache.pn` (which should already
/// contain base loads).  Does not reset `pn` -- caller handles that.
fn accumulate_self_weight_loads(
    cache: &mut FdmCache,
    gravity: &[f64; 3],
) {
    let ne = cache.member_lengths.len();
    for k in 0..ne {
        let mu_k = cache.sw_mu[k];
        let l_k = cache.member_lengths[k];
        let w_half = 0.5 * mu_k * l_k;

        let s = cache.edge_starts[k];
        let e = cache.edge_ends[k];

        if let Some(sf) = cache.node_to_free_idx[s] {
            for d in 0..3 {
                cache.pn[[sf, d]] += w_half * gravity[d];
            }
        }
        if let Some(ef) = cache.node_to_free_idx[e] {
            for d in 0..3 {
                cache.pn[[ef, d]] += w_half * gravity[d];
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Pressure load computation
// ─────────────────────────────────────────────────────────────

/// Compute Newell area-weighted normal for a face (un-normalized).
/// Returns `(nx, ny, nz)` where `||(nx,ny,nz)|| = face area`.
pub fn newell_normal(face: &[usize], nf: &Array2<f64>) -> [f64; 3] {
    let nv = face.len();
    let mut n = [0.0f64; 3];
    for i in 0..nv {
        let j = (i + 1) % nv;
        let vi = face[i];
        let vj = face[j];
        // Cross product contribution: (v_i × v_j)
        n[0] += (nf[[vi, 1]] - nf[[vj, 1]]) * (nf[[vi, 2]] + nf[[vj, 2]]);
        n[1] += (nf[[vi, 2]] - nf[[vj, 2]]) * (nf[[vi, 0]] + nf[[vj, 0]]);
        n[2] += (nf[[vi, 0]] - nf[[vj, 0]]) * (nf[[vi, 1]] + nf[[vj, 1]]);
    }
    // Newell gives 2× the area-weighted normal
    n[0] *= 0.5;
    n[1] *= 0.5;
    n[2] *= 0.5;
    n
}

/// Add pressure loads into `cache.pn` (which should already contain
/// base loads and any self-weight).  Does not reset `pn`.
/// Dispatches on the pressure mode (Normal / Hydrostatic / Directional).
fn accumulate_pressure_loads(
    cache: &mut FdmCache,
    pressure: &PressureParams,
) {
    let faces = &pressure.face_topology().faces;
    match pressure {
        PressureParams::Normal { pressures, .. } => {
            for (f_idx, face) in faces.iter().enumerate() {
                let p_f = pressures[f_idx];
                let n = newell_normal(face, &cache.nf);
                let nv = face.len() as f64;
                for &vi in face {
                    if let Some(fi) = cache.node_to_free_idx[vi] {
                        for d in 0..3 {
                            cache.pn[[fi, d]] += p_f * n[d] / nv;
                        }
                    }
                }
            }
        }
        PressureParams::Hydrostatic {
            rho_fluid, g_magnitude, z_datum, up_direction, ..
        } => {
            for face in faces {
                let nv = face.len() as f64;
                if nv < 3.0 { continue; }

                // Face centroid projected onto the "up" direction
                let mut centroid_up = 0.0;
                for &vi in face {
                    for d in 0..3 {
                        centroid_up += cache.nf[[vi, d]] * up_direction[d];
                    }
                }
                centroid_up /= nv;

                let depth = z_datum - centroid_up;
                if depth <= 0.0 { continue; }

                let p_f = rho_fluid * g_magnitude * depth;
                let n = newell_normal(face, &cache.nf);
                for &vi in face {
                    if let Some(fi) = cache.node_to_free_idx[vi] {
                        for d in 0..3 {
                            cache.pn[[fi, d]] += p_f * n[d] / nv;
                        }
                    }
                }
            }
        }
        PressureParams::Directional {
            pressures, direction, ..
        } => {
            for (f_idx, face) in faces.iter().enumerate() {
                let p_f = pressures[f_idx];
                let n = newell_normal(face, &cache.nf);
                let nv = face.len() as f64;

                // Projected area = n_f · d_hat
                let a_proj: f64 = (0..3).map(|d| n[d] * direction[d]).sum();
                if a_proj <= 0.0 { continue; }

                let load_mag = p_f * a_proj / nv;
                for &vi in face {
                    if let Some(fi) = cache.node_to_free_idx[vi] {
                        for d in 0..3 {
                            cache.pn[[fi, d]] += load_mag * direction[d];
                        }
                    }
                }
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Universal forward solve with geometry-dependent loads
// ─────────────────────────────────────────────────────────────

/// Forward solve that handles self-weight and/or pressure loads via iteration.
///
/// When neither `problem.self_weight` nor `problem.pressure` is set, this
/// delegates directly to [`solve_fdm`] with zero overhead.
pub fn solve_fdm_with_loads(
    cache: &mut FdmCache,
    q: &[f64],
    problem: &Problem,
    anchor_positions: &Array2<f64>,
    perturbation: f64,
) -> Result<(), TheseusError> {
    let has_sw = problem.self_weight.is_some();
    let has_pressure = problem.pressure.is_some();

    if !has_sw && !has_pressure {
        return solve_fdm(cache, q, problem, anchor_positions, perturbation);
    }

    // Determine iteration parameters from whichever load type is active
    let max_iters = match (&problem.self_weight, &problem.pressure) {
        (Some(sw), Some(pr)) => sw.max_iters().max(pr.max_iters()),
        (Some(sw), None) => sw.max_iters(),
        (None, Some(pr)) => pr.max_iters(),
        (None, None) => unreachable!(),
    };
    let tolerance = match (&problem.self_weight, &problem.pressure) {
        (Some(sw), Some(pr)) => sw.tolerance().min(pr.tolerance()),
        (Some(sw), None) => sw.tolerance(),
        (None, Some(pr)) => pr.tolerance(),
        (None, None) => unreachable!(),
    };
    let relaxation = match &problem.self_weight {
        Some(sw) => sw.relaxation(),
        None => 1.0,
    };

    // Reset pn to base loads
    cache.pn.assign(&cache.pn_base);

    let nn_free = problem.topology.free_node_indices.len();

    for _iter in 0..max_iters {
        // Inner linear FDM solve with current loads
        solve_fdm(cache, q, problem, anchor_positions, perturbation)?;

        // Update sizing mu from current forces if applicable
        if let Some(SelfWeightParams::Sizing { rho, sigma, .. }) = &problem.self_weight {
            update_sizing_mu(cache, *rho, *sigma);
        }

        // Save current pn, then rebuild from base + geometry-dependent loads
        let pn_old = cache.pn.clone();
        cache.pn.assign(&cache.pn_base);

        if let Some(sw) = &problem.self_weight {
            accumulate_self_weight_loads(cache, sw.gravity());
        }
        if let Some(pr) = &problem.pressure {
            accumulate_pressure_loads(cache, pr);
        }

        // Convergence: relative change in load vector
        let mut delta_sq = 0.0;
        let mut pn_norm_sq = 0.0;
        for i in 0..nn_free {
            for d in 0..3 {
                let diff = cache.pn[[i, d]] - pn_old[[i, d]];
                delta_sq += diff * diff;
                pn_norm_sq += cache.pn[[i, d]] * cache.pn[[i, d]];
            }
        }
        if pn_norm_sq > 0.0 && delta_sq / pn_norm_sq < tolerance * tolerance {
            break;
        }

        // Relaxation: blend new loads with old
        if relaxation < 1.0 {
            for i in 0..nn_free {
                for d in 0..3 {
                    cache.pn[[i, d]] = relaxation * cache.pn[[i, d]]
                        + (1.0 - relaxation) * pn_old[[i, d]];
                }
            }
        }
    }

    // Final solve with converged loads
    solve_fdm(cache, q, problem, anchor_positions, perturbation)
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
