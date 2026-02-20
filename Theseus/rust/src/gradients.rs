//! Hand-coded gradients for all 13 objectives + adjoint solve.
//!
//! The adjoint method computes dJ/dq via:
//!   1. Accumulate explicit dJ/dx̂ from each objective
//!   2. Solve adjoint system A^T λ = dJ/dx̂  (reuse A Factorization)
//!   3. Implicit gradient:  dJ/dq_k = −Δλ_k · ΔN_k
//!   4. Add explicit dJ/dq_k terms (SumForceLength, barrier, etc.)
//!
//! All gradients derived analytically — no AD framework needed.

use crate::objectives::{softplus_grad, bounds_penalty_grad};
use crate::types::{FdmCache, GeometrySnapshot, Problem, TheseusError};
use ndarray::Array2;

// ─────────────────────────────────────────────────────────────
//  Adjoint solve  (reuses LDL Factorization from forward solve)
// ─────────────────────────────────────────────────────────────

/// Solve A λ = dJ/dx̂ for each coordinate column.
///
/// Since A is symmetric (A = Aᵀ), we reuse the **same** factorization
/// (Cholesky or LDL) from the forward solve — no refactoring needed.
pub fn solve_adjoint(cache: &mut FdmCache) -> Result<(), TheseusError> {
    let n = cache.a_matrix.cols();

    for d in 0..3 {
        let rhs: Vec<f64> = (0..n).map(|i| cache.grad_x[[i, d]]).collect();
        let x = cache.factorization.as_ref()
            .ok_or(TheseusError::MissingFactorization)?
            .solve(&rhs);
        for i in 0..n {
            cache.lambda[[i, d]] = x[i];
        }
    }

    Ok(())
}

// ─────────────────────────────────────────────────────────────
//  Implicit gradient:  dJ/dq_k (from adjoint)
// ─────────────────────────────────────────────────────────────

/// Accumulate implicit gradient dJ/dq_k = −Δλ_k · ΔN_k
/// and dJ/dNf contributions for variable anchors.
pub fn accumulate_implicit_gradients(cache: &mut FdmCache, problem: &Problem) {
    let ne = problem.topology.num_edges;

    for k in 0..ne {
        let u = cache.edge_starts[k];
        let v = cache.edge_ends[k];

        let u_free = cache.node_to_free_idx[u];
        let v_free = cache.node_to_free_idx[v];

        for d in 0..3 {
            let lam_u = if let Some(uf) = u_free { cache.lambda[[uf, d]] } else { 0.0 };
            let lam_v = if let Some(vf) = v_free { cache.lambda[[vf, d]] } else { 0.0 };
            let d_lam = lam_v - lam_u;

            let d_n = cache.nf[[v, d]] - cache.nf[[u, d]];

            // Implicit ∂J/∂q_k
            cache.grad_q[k] -= d_lam * d_n;

            // ∂J/∂Nf  (fixed-node contributions)
            let term = -cache.q[k] * d_lam;
            if v_free.is_none() {
                cache.grad_nf[[v, d]] += term;
            }
            if u_free.is_none() {
                cache.grad_nf[[u, d]] -= term;
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Explicit  dJ/dx̂  for each objective type
// ─────────────────────────────────────────────────────────────

/// Zero and accumulate dJ/dx̂ (grad_x) and explicit dJ/dq from all objectives.
///
/// After this, `cache.grad_x` is ready for the adjoint solve.
pub fn accumulate_explicit_gradients(
    cache: &mut FdmCache,
    problem: &Problem,
) {
    cache.grad_x.fill(0.0);

    // We also need to zero grad_q here because the adjoint adds to it later.
    // Explicit dJ/dq terms are accumulated inline.
    // (grad_q will be zeroed in value_and_gradient before this is called)

    for obj in &problem.objectives {
        obj.accumulate_gradient(cache, problem);
    }
}

// ─────────────────────────────────────────────────────────────
//  Per-objective gradient implementations
// ─────────────────────────────────────────────────────────────

/// TargetXYZ:  L = w Σ_i ‖xyz[idx] − t_i‖²
/// dL/dx̂[j,d] = 2w (xyz[idx,d] − t[i,d])  if idx is a free node at position j
pub(crate) fn grad_target_xyz(
    cache: &mut FdmCache,
    weight: f64,
    node_indices: &[usize],
    target: &Array2<f64>,
    _free_node_indices: &[usize],
) {
    for (i, &idx) in node_indices.iter().enumerate() {
        if let Some(j) = cache.node_to_free_idx[idx] {
            for d in 0..3 {
                cache.grad_x[[j, d]] += 2.0 * weight * (cache.nf[[idx, d]] - target[[i, d]]);
            }
        }
        // Fixed node contributions go to grad_nf (handled by adjoint accumulation)
    }
}

/// TargetXY:  only x,y dimensions
pub(crate) fn grad_target_xy(
    cache: &mut FdmCache,
    weight: f64,
    node_indices: &[usize],
    target: &Array2<f64>,
    _free_node_indices: &[usize],
) {
    for (i, &idx) in node_indices.iter().enumerate() {
        if let Some(j) = cache.node_to_free_idx[idx] {
            for d in 0..2 {
                cache.grad_x[[j, d]] += 2.0 * weight * (cache.nf[[idx, d]] - target[[i, d]]);
            }
        }
    }
}

/// TargetLength:  L = w Σ (ℓ_k − t_k)²
/// dL/dx̂ via chain rule through ℓ_k = ‖ΔN_k‖
///   dℓ/dx̂[j,d] = ΔN_k[d] / ℓ_k  ×  (±1 depending on edge orientation)
///   dL/dx̂ += 2w (ℓ_k − t_k) · dℓ/dx̂
pub(crate) fn grad_target_length(
    cache: &mut FdmCache,
    weight: f64,
    edge_indices: &[usize],
    target: &[f64],
) {
    for (i, &k) in edge_indices.iter().enumerate() {
        let s = cache.edge_starts[k];
        let e = cache.edge_ends[k];
        let len = cache.member_lengths[k];
        if len < f64::EPSILON { continue; }

        let scale = 2.0 * weight * (len - target[i]) / len;

        let s_free = cache.node_to_free_idx[s];
        let e_free = cache.node_to_free_idx[e];

        for d in 0..3 {
            let delta = cache.nf[[e, d]] - cache.nf[[s, d]];
            if let Some(ef) = e_free {
                cache.grad_x[[ef, d]] += scale * delta;
            }
            if let Some(sf) = s_free {
                cache.grad_x[[sf, d]] -= scale * delta;
            }
        }
    }
}

/// Softmax weights: w_i = exp(β(v_i − m)) / Σ exp(β(v_j − m))
/// Returns one weight per edge in `edge_indices`.
fn softmax_weights(values: &[f64], edge_indices: &[usize], beta: f64) -> Vec<f64> {
    let m = edge_indices.iter().map(|&i| values[i]).fold(f64::NEG_INFINITY, f64::max);
    let exps: Vec<f64> = edge_indices.iter().map(|&i| ((values[i] - m) * beta).exp()).collect();
    let sum: f64 = exps.iter().sum();
    exps.into_iter().map(|e| e / sum).collect()
}

/// LengthVariation:  L = w (smooth_max(ℓ) − smooth_min(ℓ))
///   dL/dℓ_i = w (softmax_i(β ℓ) − softmax_i(−β ℓ))
/// Then chain through ℓ → x̂.
pub(crate) fn grad_length_variation(
    cache: &mut FdmCache,
    weight: f64,
    edge_indices: &[usize],
    beta: f64,
) {
    if edge_indices.is_empty() { return; }

    // softmax for smooth_max  (positive β)
    let w_max = softmax_weights(&cache.member_lengths, edge_indices, beta);
    // softmax for smooth_min = −smooth_max(−v):  d(smooth_min)/dv_i = softmax_i(−β v)
    let w_min = softmax_weights(&cache.member_lengths, edge_indices, -beta);

    for (j, &k) in edge_indices.iter().enumerate() {
        let dl_dl = weight * (w_max[j] - w_min[j]);
        add_length_grad_to_x(cache, k, dl_dl);
    }
}

/// ForceVariation:  L = w (smooth_max(f) − smooth_min(f))
///   dL/df_i = w (softmax_i(β f) − softmax_i(−β f))
/// Then chain through f → (x̂, q).
pub(crate) fn grad_force_variation(
    cache: &mut FdmCache,
    weight: f64,
    edge_indices: &[usize],
    beta: f64,
) {
    if edge_indices.is_empty() { return; }

    let w_max = softmax_weights(&cache.member_forces, edge_indices, beta);
    let w_min = softmax_weights(&cache.member_forces, edge_indices, -beta);

    for (j, &k) in edge_indices.iter().enumerate() {
        let dl_df = weight * (w_max[j] - w_min[j]);
        add_force_grad(cache, k, dl_df);
    }
}

/// SumForceLength:  L = w Σ ℓ_k · f_k = w Σ q_k ℓ_k²
/// dL/dx̂ via ℓ_k:  2w q_k ℓ_k · dℓ_k/dx̂
/// dL/dq_k = w ℓ_k²  (explicit)
pub(crate) fn grad_sum_force_length(
    cache: &mut FdmCache,
    weight: f64,
    edge_indices: &[usize],
) {
    for &k in edge_indices {
        let len = cache.member_lengths[k];
        let qk = cache.q[k];

        // Explicit dJ/dq_k
        cache.grad_q[k] += weight * len * len;

        // dJ/dx̂ through ℓ_k
        let scale = 2.0 * weight * qk;  // factor out ℓ_k/ℓ_k cancel with direction
        add_length_grad_to_x(cache, k, scale * len);
    }
}

/// MinLength barrier:  L = w Σ softplus(ℓ_k, threshold_k, −k_sharp)
pub(crate) fn grad_min_length(
    cache: &mut FdmCache,
    weight: f64,
    edge_indices: &[usize],
    threshold: &[f64],
    sharpness: f64,
) {
    for (i, &k) in edge_indices.iter().enumerate() {
        if !threshold[i].is_finite() { continue; }
        let dsp = softplus_grad(cache.member_lengths[k], threshold[i], -sharpness);
        add_length_grad_to_x(cache, k, weight * dsp);
    }
}

/// MaxLength barrier
pub(crate) fn grad_max_length(
    cache: &mut FdmCache,
    weight: f64,
    edge_indices: &[usize],
    threshold: &[f64],
    sharpness: f64,
) {
    for (i, &k) in edge_indices.iter().enumerate() {
        if !threshold[i].is_finite() { continue; }
        let dsp = softplus_grad(cache.member_lengths[k], threshold[i], sharpness);
        add_length_grad_to_x(cache, k, weight * dsp);
    }
}

/// MinForce barrier:  chain through f_k = q_k ℓ_k
pub(crate) fn grad_min_force(
    cache: &mut FdmCache,
    weight: f64,
    edge_indices: &[usize],
    threshold: &[f64],
    sharpness: f64,
) {
    for (i, &k) in edge_indices.iter().enumerate() {
        if !threshold[i].is_finite() { continue; }
        let dsp = softplus_grad(cache.member_forces[k], threshold[i], -sharpness);
        add_force_grad(cache, k, weight * dsp);
    }
}

/// MaxForce barrier
pub(crate) fn grad_max_force(
    cache: &mut FdmCache,
    weight: f64,
    edge_indices: &[usize],
    threshold: &[f64],
    sharpness: f64,
) {
    for (i, &k) in edge_indices.iter().enumerate() {
        if !threshold[i].is_finite() { continue; }
        let dsp = softplus_grad(cache.member_forces[k], threshold[i], sharpness);
        add_force_grad(cache, k, weight * dsp);
    }
}

/// RigidSetCompare:  L = w Σ_{i<j} (d_tgt − d_net)²
pub(crate) fn grad_rigid_set_compare(
    cache: &mut FdmCache,
    weight: f64,
    node_indices: &[usize],
    target: &Array2<f64>,
    _free_node_indices: &[usize],
) {
    let n = node_indices.len();
    for i in 0..n {
        let idx_i = node_indices[i];
        for j in 0..i {
            let idx_j = node_indices[j];

            let mut delta = [0.0f64; 3];
            let mut t_delta = [0.0f64; 3];
            for d in 0..3 {
                delta[d] = cache.nf[[idx_i, d]] - cache.nf[[idx_j, d]];
                t_delta[d] = target[[i, d]] - target[[j, d]];
            }
            let d_net = (delta[0].powi(2) + delta[1].powi(2) + delta[2].powi(2)).sqrt();
            let d_tgt = (t_delta[0].powi(2) + t_delta[1].powi(2) + t_delta[2].powi(2)).sqrt();

            if d_net < f64::EPSILON { continue; }

            // dL/d(d_net) = −2w (d_tgt − d_net)
            // d(d_net)/dx̂[v, d] = delta[d] / d_net  (for node idx_i), negative for idx_j
            let scale = -2.0 * weight * (d_tgt - d_net) / d_net;

            let fi = cache.node_to_free_idx[idx_i];
            let fj = cache.node_to_free_idx[idx_j];

            for d in 0..3 {
                if let Some(fi_val) = fi {
                    cache.grad_x[[fi_val, d]] += scale * delta[d];
                }
                if let Some(fj_val) = fj {
                    cache.grad_x[[fj_val, d]] -= scale * delta[d];
                }
            }
        }
    }
}

/// ReactionDirection:  L = w Σ (1 − r̂·d̂)
/// 
/// Gradient through reactions → q and x̂ is complex.
/// Reactions depend on q and geometry: R_v = Σ_{k∈star(v)} q_k (x_e − x_s)
/// This contribution flows through the adjoint automatically since reactions
/// are computed from q and x̂.  We push dL/dReactions into grad_x via the
/// chain rule.
pub(crate) fn grad_reaction_direction(
    cache: &mut FdmCache,
    problem: &Problem,
    weight: f64,
    anchor_indices: &[usize],
    target_directions: &Array2<f64>,
) {
    // Compute dL/dR[node, d] for each anchor, then chain to dL/dx̂ and dL/dq
    for (row, &node) in anchor_indices.iter().enumerate() {
        let r = [cache.reactions[[node, 0]], cache.reactions[[node, 1]], cache.reactions[[node, 2]]];
        let d_hat = [target_directions[[row, 0]], target_directions[[row, 1]], target_directions[[row, 2]]];

        let r_norm = (r[0].powi(2) + r[1].powi(2) + r[2].powi(2)).sqrt();
        if r_norm < f64::EPSILON { continue; }

        let dot: f64 = (0..3).map(|d| r[d] * d_hat[d]).sum();
        let cos = dot / r_norm;

        // d(1 − cos)/dR[d] = −(d̂[d]/‖R‖ − cos·R[d]/‖R‖²)
        //                   = −(d̂[d] − cos·R̂[d]) / ‖R‖
        let mut dl_dr = [0.0f64; 3];
        for d in 0..3 {
            dl_dr[d] = -weight * (d_hat[d] - cos * r[d] / r_norm) / r_norm;
        }

        // Chain dL/dR → dL/dx̂ and dL/dq through reaction formula:
        //   R[node, d] = Σ_k  sign_k · q_k · (x_e[d] − x_s[d])
        accumulate_reaction_grad(cache, problem, node, &dl_dr);
    }
}

/// ReactionDirectionMagnitude:  adds magnitude penalty
pub(crate) fn grad_reaction_direction_magnitude(
    cache: &mut FdmCache,
    problem: &Problem,
    weight: f64,
    anchor_indices: &[usize],
    target_directions: &Array2<f64>,
    target_magnitudes: &[f64],
) {
    for (row, &node) in anchor_indices.iter().enumerate() {
        let r = [cache.reactions[[node, 0]], cache.reactions[[node, 1]], cache.reactions[[node, 2]]];
        let d_hat = [target_directions[[row, 0]], target_directions[[row, 1]], target_directions[[row, 2]]];

        let r_norm = (r[0].powi(2) + r[1].powi(2) + r[2].powi(2)).sqrt();
        if r_norm < f64::EPSILON { continue; }

        let dot: f64 = (0..3).map(|d| r[d] * d_hat[d]).sum();
        let cos = dot / r_norm;

        // Direction gradient
        let mut dl_dr = [0.0f64; 3];
        for d in 0..3 {
            dl_dr[d] = -weight * (d_hat[d] - cos * r[d] / r_norm) / r_norm;
        }

        // Magnitude gradient:  d max(‖R‖ − m, 0) / dR  =  R̂  if ‖R‖ > m
        if r_norm > target_magnitudes[row] {
            for d in 0..3 {
                dl_dr[d] += weight * r[d] / r_norm;
            }
        }

        accumulate_reaction_grad(cache, problem, node, &dl_dr);
    }
}

// ─────────────────────────────────────────────────────────────
//  Helpers
// ─────────────────────────────────────────────────────────────

/// Add dL/dℓ_k · (dℓ_k/dx̂) into cache.grad_x.
///   dℓ_k/dx̂[e_free, d] = ΔN_k[d] / ℓ_k
///   dℓ_k/dx̂[s_free, d] = −ΔN_k[d] / ℓ_k
fn add_length_grad_to_x(cache: &mut FdmCache, k: usize, dl_dl: f64) {
    let s = cache.edge_starts[k];
    let e = cache.edge_ends[k];
    let len = cache.member_lengths[k];
    if len < f64::EPSILON { return; }

    let s_free = cache.node_to_free_idx[s];
    let e_free = cache.node_to_free_idx[e];

    for d in 0..3 {
        let delta = cache.nf[[e, d]] - cache.nf[[s, d]];
        let g = dl_dl * delta / len;
        if let Some(ef) = e_free {
            cache.grad_x[[ef, d]] += g;
        }
        if let Some(sf) = s_free {
            cache.grad_x[[sf, d]] -= g;
        }
    }
}

/// Add dL/df_k · (df_k/dx̂, df_k/dq_k) into grad_x and grad_q.
/// f_k = q_k ℓ_k  →  df/dx̂ = q_k dℓ/dx̂,  df/dq_k = ℓ_k
fn add_force_grad(cache: &mut FdmCache, k: usize, dl_df: f64) {
    // Explicit dJ/dq_k
    cache.grad_q[k] += dl_df * cache.member_lengths[k];
    // dJ/dx̂ through ℓ
    add_length_grad_to_x(cache, k, dl_df * cache.q[k]);
}

/// Chain dL/dR[node] through the reaction formula to dL/dx̂ and dL/dq.
///
/// Reaction at `node` from edge k:
///   if node == edge_start:  R += q_k * (x_end − x_start)  →  sign = +1
///   if node == edge_end:    R -= q_k * (x_end − x_start)  →  sign = −1
fn accumulate_reaction_grad(
    cache: &mut FdmCache,
    problem: &Problem,
    node: usize,
    dl_dr: &[f64; 3],
) {
    let ne = problem.topology.num_edges;
    for k in 0..ne {
        let s = cache.edge_starts[k];
        let e = cache.edge_ends[k];
        let qi = cache.q[k];

        let sign = if s == node {
            1.0
        } else if e == node {
            -1.0
        } else {
            continue;
        };

        // dR/dq_k  = sign * (x_e − x_s)
        let mut dq_contrib = 0.0;
        for d in 0..3 {
            let delta = cache.nf[[e, d]] - cache.nf[[s, d]];
            dq_contrib += dl_dr[d] * sign * delta;
        }
        cache.grad_q[k] += dq_contrib;

        // dR/dx̂  = sign * q_k * I  (for end node),  −sign*q_k*I (for start)
        let e_free = cache.node_to_free_idx[e];
        let s_free = cache.node_to_free_idx[s];
        for d in 0..3 {
            let g = dl_dr[d] * sign * qi;
            if let Some(ef) = e_free {
                cache.grad_x[[ef, d]] += g;
            }
            if let Some(sf) = s_free {
                cache.grad_x[[sf, d]] -= g;
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Full value-and-gradient  (the L-BFGS entry point)
// ─────────────────────────────────────────────────────────────

/// Compute both J(θ) and ∇J(θ) in one pass.
///
/// θ = [q₁..qₙ, anchor_x₁, anchor_y₁, anchor_z₁, …]
///
/// Steps:
///   1. Unpack θ into q and anchor positions
///   2. Forward solve  → geometry snapshot
///   3. Evaluate total loss J
///   4. Accumulate explicit dJ/dx̂ from objectives
///   5. Adjoint solve  A λ = dJ/dx̂
///   6. Implicit dJ/dq  += −Δλ · ΔN
///   7. Barrier gradient on θ
///   8. Pack grad_q + grad_anchors → grad vector
pub fn value_and_gradient(
    cache: &mut FdmCache,
    problem: &Problem,
    theta: &[f64],
    grad: &mut [f64],
    lb: &[f64],
    ub: &[f64],
    lb_idx: &[usize],
    ub_idx: &[usize],
) -> Result<f64, TheseusError> {
    let ne = problem.topology.num_edges;
    let nvar = problem.anchors.variable_indices.len();

    // 1. Unpack
    let q = &theta[..ne];
    let anchor_data = &theta[ne..];
    let anchor_positions = if nvar > 0 {
        let mut a = Array2::<f64>::zeros((nvar, 3));
        for i in 0..nvar {
            a[[i, 0]] = anchor_data[i * 3];
            a[[i, 1]] = anchor_data[i * 3 + 1];
            a[[i, 2]] = anchor_data[i * 3 + 2];
        }
        a
    } else {
        Array2::<f64>::zeros((0, 3))
    };

    // 2. Forward solve
    crate::fdm::solve_fdm(cache, q, problem, &anchor_positions, 1e-12)?;

    // 3. Build snapshot and evaluate loss
    let snap = GeometrySnapshot {
        xyz_full: &cache.nf,
        member_lengths: &cache.member_lengths,
        member_forces: &cache.member_forces,
        reactions: &cache.reactions,
    };
    let geometric_loss = crate::objectives::total_loss(&problem.objectives, &snap);
    let barrier_loss = crate::objectives::bounds_penalty(
        theta, lb, ub, lb_idx, ub_idx, problem.solver.barrier_sharpness,
    );
    let total = geometric_loss + barrier_loss * problem.solver.barrier_weight;

    // 4. Explicit gradients (fills grad_x, partial grad_q)
    cache.grad_q.fill(0.0);
    cache.grad_nf.fill(0.0);
    accumulate_explicit_gradients(cache, problem);

    // 5. Adjoint solve
    solve_adjoint(cache)?;

    // 6. Implicit gradients
    accumulate_implicit_gradients(cache, problem);

    // 7. Pack into output gradient
    grad.fill(0.0);
    grad[..ne].copy_from_slice(&cache.grad_q);

    // Anchor gradients
    if nvar > 0 {
        for i in 0..nvar {
            let node = problem.anchors.variable_indices[i];
            grad[ne + i * 3] += cache.grad_nf[[node, 0]];
            grad[ne + i * 3 + 1] += cache.grad_nf[[node, 1]];
            grad[ne + i * 3 + 2] += cache.grad_nf[[node, 2]];
        }
    }

    // 8. Barrier gradient
    bounds_penalty_grad(
        grad, theta, lb, ub, lb_idx, ub_idx,
        problem.solver.barrier_sharpness,
        problem.solver.barrier_weight,
    );

    Ok(total)
}
