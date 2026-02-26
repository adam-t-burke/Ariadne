//! Objective loss functions — pure ℝ → ℝ math, no AD.
//!
//! Mirrors `src/objectives.jl` from the Julia code.  Each function computes
//! a scalar loss from the current geometry snapshot.  The corresponding
//! hand-coded gradients live in `gradients.rs`.

use crate::types::{
    GeometrySnapshot, ObjectiveTrait, FdmCache, Problem,
    TargetXYZ, TargetXY, TargetPlane, TargetLength, LengthVariation, ForceVariation,
    SumForceLength, MinLength, MaxLength, MinForce, MaxForce,
    RigidSetCompare, ReactionDirection, ReactionDirectionMagnitude,
};
use crate::gradients;
use ndarray::Array2;

// ─────────────────────────────────────────────────────────────
//  Softplus barrier
// ─────────────────────────────────────────────────────────────

/// Numerically stable log(1 + exp(z)).
#[inline]
fn log1pexp(z: f64) -> f64 {
    if z > 0.0 {
        z + (-z).exp().ln_1p()
    } else {
        z.exp().ln_1p()
    }
}

/// Smooth one-sided barrier.
/// `k < 0` ⟹  penalise x < b  (min barrier).
/// `k > 0` ⟹  penalise x > b  (max barrier).
#[inline]
pub fn softplus(x: f64, b: f64, k: f64) -> f64 {
    let z = -k * (b - x) - 1.0;
    log1pexp(z)
}

/// Derivative of `softplus` w.r.t. `x`.
/// d/dx softplus = k · σ(z)  where z = −k(b−x)−1 and σ is the logistic fn.
#[inline]
pub fn softplus_grad(x: f64, b: f64, k: f64) -> f64 {
    let z = -k * (b - x) - 1.0;
    let sigma = 1.0 / (1.0 + (-z).exp());
    k * sigma
}

// ─────────────────────────────────────────────────────────────
//  Bound penalties  (barrier on the packed θ vector)
// ─────────────────────────────────────────────────────────────

/// Smooth barrier loss for lower and upper bounds on θ.
pub fn bounds_penalty(theta: &[f64], lb: &[f64], ub: &[f64], lb_idx: &[usize], ub_idx: &[usize], sharpness: f64) -> f64 {
    let mut loss = 0.0;
    for &i in lb_idx {
        if lb[i].is_finite() {
            loss += softplus(theta[i], lb[i], -sharpness);
        }
    }
    for &i in ub_idx {
        if ub[i].is_finite() {
            loss += softplus(theta[i], ub[i], sharpness);
        }
    }
    loss
}

/// Gradient of `bounds_penalty` w.r.t. θ.  Accumulates into `grad`.
pub fn bounds_penalty_grad(
    grad: &mut [f64],
    theta: &[f64],
    lb: &[f64],
    ub: &[f64],
    lb_idx: &[usize],
    ub_idx: &[usize],
    sharpness: f64,
    barrier_weight: f64,
) {
    for &i in lb_idx {
        if lb[i].is_finite() {
            grad[i] += barrier_weight * softplus_grad(theta[i], lb[i], -sharpness);
        }
    }
    for &i in ub_idx {
        if ub[i].is_finite() {
            grad[i] += barrier_weight * softplus_grad(theta[i], ub[i], sharpness);
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Individual objective losses
// ─────────────────────────────────────────────────────────────

/// TargetXYZ:  Σ_i ‖xyz[idx_i] − target_i‖²
fn target_xyz_loss(xyz: &Array2<f64>, node_indices: &[usize], target: &Array2<f64>) -> f64 {
    let mut loss = 0.0;
    for (i, &idx) in node_indices.iter().enumerate() {
        for d in 0..3 {
            let diff = xyz[[idx, d]] - target[[i, d]];
            loss += diff * diff;
        }
    }
    loss
}

/// TargetXY:  Σ_i (Δx² + Δy²)
fn target_xy_loss(xyz: &Array2<f64>, node_indices: &[usize], target: &Array2<f64>) -> f64 {
    let mut loss = 0.0;
    for (i, &idx) in node_indices.iter().enumerate() {
        for d in 0..2 {
            let diff = xyz[[idx, d]] - target[[i, d]];
            loss += diff * diff;
        }
    }
    loss
}

/// TargetPlane:  Σ_i ((u_p − u_t)² + (v_p − v_t)²) in plane coords.
/// u = (p − origin)·x_axis, v = (p − origin)·y_axis.
fn target_plane_loss(
    xyz: &Array2<f64>,
    node_indices: &[usize],
    target: &Array2<f64>,
    origin: &[f64; 3],
    x_axis: &[f64; 3],
    y_axis: &[f64; 3],
) -> f64 {
    let mut loss = 0.0;
    for (i, &idx) in node_indices.iter().enumerate() {
        let u_p = (xyz[[idx, 0]] - origin[0]) * x_axis[0]
            + (xyz[[idx, 1]] - origin[1]) * x_axis[1]
            + (xyz[[idx, 2]] - origin[2]) * x_axis[2];
        let v_p = (xyz[[idx, 0]] - origin[0]) * y_axis[0]
            + (xyz[[idx, 1]] - origin[1]) * y_axis[1]
            + (xyz[[idx, 2]] - origin[2]) * y_axis[2];
        let u_t = (target[[i, 0]] - origin[0]) * x_axis[0]
            + (target[[i, 1]] - origin[1]) * x_axis[1]
            + (target[[i, 2]] - origin[2]) * x_axis[2];
        let v_t = (target[[i, 0]] - origin[0]) * y_axis[0]
            + (target[[i, 1]] - origin[1]) * y_axis[1]
            + (target[[i, 2]] - origin[2]) * y_axis[2];
        let du = u_p - u_t;
        let dv = v_p - v_t;
        loss += du * du + dv * dv;
    }
    loss
}

/// TargetLength:  Σ_i (ℓ[idx_i] − target_i)²
fn target_length_loss(lengths: &[f64], edge_indices: &[usize], target: &[f64]) -> f64 {
    let mut loss = 0.0;
    for (i, &idx) in edge_indices.iter().enumerate() {
        let diff = lengths[idx] - target[i];
        loss += diff * diff;
    }
    loss
}

/// Numerically stable log-sum-exp smooth maximum:
///   smooth_max_β(v) = m + (1/β) ln Σ exp(β(v_i − m)),  m = max(v)
/// As β → ∞ this converges to the true max.
fn smooth_max(values: &[f64], indices: &[usize], beta: f64) -> f64 {
    let m = indices.iter().map(|&i| values[i]).fold(f64::NEG_INFINITY, f64::max);
    let sum: f64 = indices.iter().map(|&i| ((values[i] - m) * beta).exp()).sum();
    m + sum.ln() / beta
}

/// Numerically stable log-sum-exp smooth minimum:
///   smooth_min_β(v) = −smooth_max_β(−v)
fn smooth_min(values: &[f64], indices: &[usize], beta: f64) -> f64 {
    let m = indices.iter().map(|&i| values[i]).fold(f64::INFINITY, f64::min);
    let sum: f64 = indices.iter().map(|&i| ((m - values[i]) * beta).exp()).sum();
    m - sum.ln() / beta
}

/// Minimum sharpness to avoid division-by-zero in smooth_max/smooth_min (1/β term).
const MIN_VARIATION_SHARPNESS: f64 = 1e-10;

/// LengthVariation:  smooth_max(ℓ) − smooth_min(ℓ)  over selected edges.
fn length_variation_loss(lengths: &[f64], edge_indices: &[usize], beta: f64) -> f64 {
    if edge_indices.is_empty() {
        return 0.0;
    }
    let beta = beta.abs().max(MIN_VARIATION_SHARPNESS);
    smooth_max(lengths, edge_indices, beta) - smooth_min(lengths, edge_indices, beta)
}

/// ForceVariation:  smooth_max(f) − smooth_min(f)  over selected edges.
fn force_variation_loss(forces: &[f64], edge_indices: &[usize], beta: f64) -> f64 {
    if edge_indices.is_empty() {
        return 0.0;
    }
    let beta = beta.abs().max(MIN_VARIATION_SHARPNESS);
    smooth_max(forces, edge_indices, beta) - smooth_min(forces, edge_indices, beta)
}

/// SumForceLength:  Σ_i ℓ_i · f_i  =  Σ_i q_i · ℓ_i²
fn sum_force_length_loss(lengths: &[f64], forces: &[f64], edge_indices: &[usize]) -> f64 {
    let mut loss = 0.0;
    for &idx in edge_indices {
        loss += lengths[idx] * forces[idx];
    }
    loss
}

/// MinLength / MinForce barrier:  Σ softplus(x_i, threshold_i, −k)
fn min_penalty(values: &[f64], edge_indices: &[usize], threshold: &[f64], k: f64) -> f64 {
    let mut loss = 0.0;
    for (i, &idx) in edge_indices.iter().enumerate() {
        if threshold[i].is_finite() {
            loss += softplus(values[idx], threshold[i], -k);
        }
    }
    loss
}

/// MaxLength / MaxForce barrier:  Σ softplus(x_i, threshold_i, +k)
fn max_penalty(values: &[f64], edge_indices: &[usize], threshold: &[f64], k: f64) -> f64 {
    let mut loss = 0.0;
    for (i, &idx) in edge_indices.iter().enumerate() {
        if threshold[i].is_finite() {
            loss += softplus(values[idx], threshold[i], k);
        }
    }
    loss
}

/// RigidSetCompare: Σ_{i<j} (d_target − d_network)²
/// Uses max(0, …) before sqrt to avoid NaN from floating-point negative squared distance.
fn rigid_set_compare_loss(xyz: &Array2<f64>, node_indices: &[usize], target: &Array2<f64>) -> f64 {
    let n = node_indices.len();
    let mut loss = 0.0;
    for i in 0..n {
        let idx_i = node_indices[i];
        for j in 0..i {
            let idx_j = node_indices[j];

            let dx = xyz[[idx_i, 0]] - xyz[[idx_j, 0]];
            let dy = xyz[[idx_i, 1]] - xyz[[idx_j, 1]];
            let dz = xyz[[idx_i, 2]] - xyz[[idx_j, 2]];
            let d_net = (dx * dx + dy * dy + dz * dz).max(0.0).sqrt();

            let tx = target[[i, 0]] - target[[j, 0]];
            let ty = target[[i, 1]] - target[[j, 1]];
            let tz = target[[i, 2]] - target[[j, 2]];
            let d_tgt = (tx * tx + ty * ty + tz * tz).max(0.0).sqrt();

            let diff = d_tgt - d_net;
            loss += diff * diff;
        }
    }
    loss
}

/// Cosine-based misalignment:  1 − (r̂ · d̂).
/// Returns 1.0 for zero reaction (full penalty).
#[inline]
fn direction_misalignment(reaction: [f64; 3], target_dir: [f64; 3]) -> f64 {
    let r_norm = (reaction[0].powi(2) + reaction[1].powi(2) + reaction[2].powi(2)).sqrt();
    if r_norm < f64::EPSILON {
        return 1.0;
    }
    let dot: f64 = (0..3).map(|d| reaction[d] * target_dir[d]).sum();
    let cos = (dot / r_norm).clamp(-1.0, 1.0);
    1.0 - cos
}

/// ReactionDirection:  Σ_i  (1 − r̂_i · d̂_i)
fn reaction_direction_loss(
    reactions: &Array2<f64>,
    anchor_indices: &[usize],
    target_directions: &Array2<f64>,
) -> f64 {
    let mut total = 0.0;
    for (row, &node) in anchor_indices.iter().enumerate() {
        let r = [reactions[[node, 0]], reactions[[node, 1]], reactions[[node, 2]]];
        let d = [target_directions[[row, 0]], target_directions[[row, 1]], target_directions[[row, 2]]];
        total += direction_misalignment(r, d);
    }
    total
}

/// ReactionDirectionMagnitude:  Σ_i [ (1 − r̂·d̂) + max(‖r‖ − m_i, 0) ]
fn reaction_direction_magnitude_loss(
    reactions: &Array2<f64>,
    anchor_indices: &[usize],
    target_directions: &Array2<f64>,
    target_magnitudes: &[f64],
) -> f64 {
    let mut total = 0.0;
    for (row, &node) in anchor_indices.iter().enumerate() {
        let r = [reactions[[node, 0]], reactions[[node, 1]], reactions[[node, 2]]];
        let d = [target_directions[[row, 0]], target_directions[[row, 1]], target_directions[[row, 2]]];
        let dir_loss = direction_misalignment(r, d);
        let r_norm = (r[0].powi(2) + r[1].powi(2) + r[2].powi(2)).sqrt();
        let mag_loss = (r_norm - target_magnitudes[row]).max(0.0);
        total += dir_loss + mag_loss;
    }
    total
}

// ─────────────────────────────────────────────────────────────
//  ObjectiveTrait implementations for all 13 built-in types
// ─────────────────────────────────────────────────────────────

impl ObjectiveTrait for TargetXYZ {
    fn loss(&self, snap: &GeometrySnapshot) -> f64 {
        self.weight * target_xyz_loss(snap.xyz_full, &self.node_indices, &self.target)
    }
    fn accumulate_gradient(&self, cache: &mut FdmCache, problem: &Problem) {
        gradients::grad_target_xyz(cache, self.weight, &self.node_indices, &self.target, &problem.topology.free_node_indices);
    }
    fn weight(&self) -> f64 { self.weight }
}

impl ObjectiveTrait for TargetXY {
    fn loss(&self, snap: &GeometrySnapshot) -> f64 {
        self.weight * target_xy_loss(snap.xyz_full, &self.node_indices, &self.target)
    }
    fn accumulate_gradient(&self, cache: &mut FdmCache, problem: &Problem) {
        gradients::grad_target_xy(cache, self.weight, &self.node_indices, &self.target, &problem.topology.free_node_indices);
    }
    fn weight(&self) -> f64 { self.weight }
}

impl ObjectiveTrait for TargetPlane {
    fn loss(&self, snap: &GeometrySnapshot) -> f64 {
        self.weight * target_plane_loss(
            snap.xyz_full,
            &self.node_indices,
            &self.target,
            &self.origin,
            &self.x_axis,
            &self.y_axis,
        )
    }
    fn accumulate_gradient(&self, cache: &mut FdmCache, problem: &Problem) {
        gradients::grad_target_plane(
            cache,
            self.weight,
            &self.node_indices,
            &self.target,
            &self.origin,
            &self.x_axis,
            &self.y_axis,
            &problem.topology.free_node_indices,
        );
    }
    fn weight(&self) -> f64 {
        self.weight
    }
}

impl ObjectiveTrait for TargetLength {
    fn loss(&self, snap: &GeometrySnapshot) -> f64 {
        self.weight * target_length_loss(snap.member_lengths, &self.edge_indices, &self.target)
    }
    fn accumulate_gradient(&self, cache: &mut FdmCache, _problem: &Problem) {
        gradients::grad_target_length(cache, self.weight, &self.edge_indices, &self.target);
    }
    fn weight(&self) -> f64 { self.weight }
}

impl ObjectiveTrait for LengthVariation {
    fn loss(&self, snap: &GeometrySnapshot) -> f64 {
        self.weight * length_variation_loss(snap.member_lengths, &self.edge_indices, self.sharpness)
    }
    fn accumulate_gradient(&self, cache: &mut FdmCache, _problem: &Problem) {
        gradients::grad_length_variation(cache, self.weight, &self.edge_indices, self.sharpness);
    }
    fn weight(&self) -> f64 { self.weight }
}

impl ObjectiveTrait for ForceVariation {
    fn loss(&self, snap: &GeometrySnapshot) -> f64 {
        self.weight * force_variation_loss(snap.member_forces, &self.edge_indices, self.sharpness)
    }
    fn accumulate_gradient(&self, cache: &mut FdmCache, _problem: &Problem) {
        gradients::grad_force_variation(cache, self.weight, &self.edge_indices, self.sharpness);
    }
    fn weight(&self) -> f64 { self.weight }
}

impl ObjectiveTrait for SumForceLength {
    fn loss(&self, snap: &GeometrySnapshot) -> f64 {
        self.weight * sum_force_length_loss(snap.member_lengths, snap.member_forces, &self.edge_indices)
    }
    fn accumulate_gradient(&self, cache: &mut FdmCache, _problem: &Problem) {
        gradients::grad_sum_force_length(cache, self.weight, &self.edge_indices);
    }
    fn weight(&self) -> f64 { self.weight }
}

impl ObjectiveTrait for MinLength {
    fn loss(&self, snap: &GeometrySnapshot) -> f64 {
        self.weight * min_penalty(snap.member_lengths, &self.edge_indices, &self.threshold, self.sharpness)
    }
    fn accumulate_gradient(&self, cache: &mut FdmCache, _problem: &Problem) {
        gradients::grad_min_length(cache, self.weight, &self.edge_indices, &self.threshold, self.sharpness);
    }
    fn weight(&self) -> f64 { self.weight }
}

impl ObjectiveTrait for MaxLength {
    fn loss(&self, snap: &GeometrySnapshot) -> f64 {
        self.weight * max_penalty(snap.member_lengths, &self.edge_indices, &self.threshold, self.sharpness)
    }
    fn accumulate_gradient(&self, cache: &mut FdmCache, _problem: &Problem) {
        gradients::grad_max_length(cache, self.weight, &self.edge_indices, &self.threshold, self.sharpness);
    }
    fn weight(&self) -> f64 { self.weight }
}

impl ObjectiveTrait for MinForce {
    fn loss(&self, snap: &GeometrySnapshot) -> f64 {
        self.weight * min_penalty(snap.member_forces, &self.edge_indices, &self.threshold, self.sharpness)
    }
    fn accumulate_gradient(&self, cache: &mut FdmCache, _problem: &Problem) {
        gradients::grad_min_force(cache, self.weight, &self.edge_indices, &self.threshold, self.sharpness);
    }
    fn weight(&self) -> f64 { self.weight }
}

impl ObjectiveTrait for MaxForce {
    fn loss(&self, snap: &GeometrySnapshot) -> f64 {
        self.weight * max_penalty(snap.member_forces, &self.edge_indices, &self.threshold, self.sharpness)
    }
    fn accumulate_gradient(&self, cache: &mut FdmCache, _problem: &Problem) {
        gradients::grad_max_force(cache, self.weight, &self.edge_indices, &self.threshold, self.sharpness);
    }
    fn weight(&self) -> f64 { self.weight }
}

impl ObjectiveTrait for RigidSetCompare {
    fn loss(&self, snap: &GeometrySnapshot) -> f64 {
        self.weight * rigid_set_compare_loss(snap.xyz_full, &self.node_indices, &self.target)
    }
    fn accumulate_gradient(&self, cache: &mut FdmCache, problem: &Problem) {
        gradients::grad_rigid_set_compare(cache, self.weight, &self.node_indices, &self.target, &problem.topology.free_node_indices);
    }
    fn weight(&self) -> f64 { self.weight }
}

impl ObjectiveTrait for ReactionDirection {
    fn loss(&self, snap: &GeometrySnapshot) -> f64 {
        self.weight * reaction_direction_loss(snap.reactions, &self.anchor_indices, &self.target_directions)
    }
    fn accumulate_gradient(&self, cache: &mut FdmCache, problem: &Problem) {
        gradients::grad_reaction_direction(cache, problem, self.weight, &self.anchor_indices, &self.target_directions);
    }
    fn weight(&self) -> f64 { self.weight }
}

impl ObjectiveTrait for ReactionDirectionMagnitude {
    fn loss(&self, snap: &GeometrySnapshot) -> f64 {
        self.weight * reaction_direction_magnitude_loss(snap.reactions, &self.anchor_indices, &self.target_directions, &self.target_magnitudes)
    }
    fn accumulate_gradient(&self, cache: &mut FdmCache, problem: &Problem) {
        gradients::grad_reaction_direction_magnitude(cache, problem, self.weight, &self.anchor_indices, &self.target_directions, &self.target_magnitudes);
    }
    fn weight(&self) -> f64 { self.weight }
}

// ─────────────────────────────────────────────────────────────
//  Dispatch:  trait-based total loss
// ─────────────────────────────────────────────────────────────

/// Evaluate total geometric loss (sum of all objectives).
pub fn total_loss(objectives: &[Box<dyn ObjectiveTrait>], snap: &GeometrySnapshot) -> f64 {
    objectives.iter().map(|obj| obj.loss(snap)).sum()
}

