//! FEA objective functions for optimization.
//!
//! Each objective implements FeaObjectiveTrait with loss and gradient methods.
//! Gradients are accumulated into grad_u (dJ/du for the adjoint) and optionally
//! into grad_area or grad_pos for explicit design-variable contributions.

use crate::fea_types::{FeaObjectiveTrait, FeaSnapshot};
use crate::objectives::softplus;

// ─────────────────────────────────────────────────────────────
//  Compliance: minimize strain energy u^T K u = u^T f
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FeaCompliance {
    pub weight: f64,
}

impl FeaObjectiveTrait for FeaCompliance {
    fn loss(&self, snap: &FeaSnapshot) -> f64 {
        // Compliance = u^T f, but we have u and can compute u^T K u
        // For linear elastic: u^T K u = u^T f = 2 * strain energy
        // We use the displacement-based form: sum over elements of F_e * delta_e
        let mut compliance = 0.0;
        for (e, &force) in snap.axial_forces.iter().enumerate() {
            let l = snap.elem_lengths[e];
            if l < f64::EPSILON { continue; }
            compliance += force * snap.strains[e] * l;
        }
        self.weight * compliance
    }

    fn accumulate_grad_u(&self, grad_u: &mut [f64], snap: &FeaSnapshot) {
        // For compliance J = u^T K u = u^T f, dJ/du = 2f = 2Ku.
        // We accumulate 2 * u into grad_u. The adjoint solve K*lambda = 2*u
        // gives lambda = 2*u (since K*u = f), and the implicit gradient
        // captures the full sensitivity.
        // Here we set dJ/du = 2*u (in full-DOF space), which after the adjoint
        // solve K*lambda = 2*u gives lambda = 2*u, and the implicit gradient
        // -lambda^T dK/dtheta u gives the correct total derivative.
        // Actually: dJ/du = 2*K*u = 2*f. We don't have K here, but we can
        // note that for the adjoint, K*lambda = dJ/du, so lambda = K^{-1} dJ/du.
        // If dJ/du = 2*f = 2*K*u, then lambda = 2*u.
        // We set grad_u[i] = 2 * displacements[i] * (something to make K*lambda = 2*K*u).
        // Simplest: just set grad_u = 2*u and note that K*(2u) = 2*f = dJ/du. Wait no.
        // K * lambda = grad_u. If grad_u = 2*K*u, then lambda = 2*u. But we can't
        // compute K*u here. Instead, we set grad_u to a proxy that the adjoint will
        // handle correctly. The cleanest: set grad_u[i] = 2 * snap.displacements[i].
        // Then K*lambda = 2*u => lambda = 2*K^{-1}*u, which is NOT 2*u.
        // The correct approach: we need grad_u = 2*f in free-DOF space.
        // But f is the RHS which we don't have in the snapshot.
        // The gradient module will handle this by using cache.rhs directly.
        // So this is intentionally a no-op for compliance.
        let _ = (grad_u, snap);
    }

    fn weight(&self) -> f64 { self.weight }

    fn has_area_gradient(&self) -> bool { true }

    fn accumulate_grad_area(&self, _grad_a: &mut [f64], _snap: &FeaSnapshot) {
        // For compliance, the area gradient is fully handled by the adjoint:
        // dJ/dA_e = -lambda^T (dke/dA_e) u, computed in fea_gradients.rs.
        // No explicit contribution needed here.
    }
}

// ─────────────────────────────────────────────────────────────
//  MaxDisplacement: smooth-max of displacement magnitudes
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FeaMaxDisplacement {
    pub weight: f64,
    pub node_indices: Vec<usize>,
}

impl FeaObjectiveTrait for FeaMaxDisplacement {
    fn loss(&self, snap: &FeaSnapshot) -> f64 {
        let beta = 20.0;
        let mut disp_mags: Vec<f64> = Vec::with_capacity(self.node_indices.len());
        for &ni in &self.node_indices {
            let ux = snap.displacements[ni * 3];
            let uy = snap.displacements[ni * 3 + 1];
            let uz = snap.displacements[ni * 3 + 2];
            disp_mags.push((ux * ux + uy * uy + uz * uz).sqrt());
        }
        if disp_mags.is_empty() { return 0.0; }
        let m = disp_mags.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let sum: f64 = disp_mags.iter().map(|&d| ((d - m) * beta).exp()).sum();
        self.weight * (m + sum.ln() / beta)
    }

    fn accumulate_grad_u(&self, grad_u: &mut [f64], snap: &FeaSnapshot) {
        let beta = 20.0;
        let mut disp_mags: Vec<f64> = Vec::with_capacity(self.node_indices.len());
        for &ni in &self.node_indices {
            let ux = snap.displacements[ni * 3];
            let uy = snap.displacements[ni * 3 + 1];
            let uz = snap.displacements[ni * 3 + 2];
            disp_mags.push((ux * ux + uy * uy + uz * uz).sqrt());
        }
        if disp_mags.is_empty() { return; }
        let m = disp_mags.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let exp_vals: Vec<f64> = disp_mags.iter().map(|&d| ((d - m) * beta).exp()).collect();
        let sum_exp: f64 = exp_vals.iter().sum();

        for (k, &ni) in self.node_indices.iter().enumerate() {
            let mag = disp_mags[k];
            if mag < f64::EPSILON { continue; }
            let softmax_weight = exp_vals[k] / sum_exp;
            for d in 0..3 {
                let u_d = snap.displacements[ni * 3 + d];
                // d(smooth_max)/du_d = softmax_weight * d(mag)/du_d
                // d(mag)/du_d = u_d / mag
                let g = self.weight * softmax_weight * u_d / mag;
                // Map to free DOF index
                grad_u[ni * 3 + d] += g;
            }
        }
    }

    fn weight(&self) -> f64 { self.weight }
}

// ─────────────────────────────────────────────────────────────
//  TargetDisplacement: L2 distance to target displacements
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FeaTargetDisplacement {
    pub weight: f64,
    pub node_indices: Vec<usize>,
    pub targets: Vec<[f64; 3]>,
}

impl FeaObjectiveTrait for FeaTargetDisplacement {
    fn loss(&self, snap: &FeaSnapshot) -> f64 {
        let mut loss = 0.0;
        for (k, &ni) in self.node_indices.iter().enumerate() {
            for d in 0..3 {
                let diff = snap.displacements[ni * 3 + d] - self.targets[k][d];
                loss += diff * diff;
            }
        }
        self.weight * loss
    }

    fn accumulate_grad_u(&self, grad_u: &mut [f64], snap: &FeaSnapshot) {
        for (k, &ni) in self.node_indices.iter().enumerate() {
            for d in 0..3 {
                let diff = snap.displacements[ni * 3 + d] - self.targets[k][d];
                grad_u[ni * 3 + d] += self.weight * 2.0 * diff;
            }
        }
    }

    fn weight(&self) -> f64 { self.weight }
}

// ─────────────────────────────────────────────────────────────
//  MinWeight: minimize total structural weight
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FeaMinWeight {
    pub weight: f64,
}

impl FeaObjectiveTrait for FeaMinWeight {
    fn loss(&self, snap: &FeaSnapshot) -> f64 {
        let mut total = 0.0;
        for e in 0..snap.elem_lengths.len() {
            total += snap.densities[e] * snap.areas[e] * snap.elem_lengths[e];
        }
        self.weight * total
    }

    fn accumulate_grad_u(&self, _grad_u: &mut [f64], _snap: &FeaSnapshot) {
        // Weight does not depend on displacements directly
    }

    fn weight(&self) -> f64 { self.weight }

    fn has_area_gradient(&self) -> bool { true }

    fn accumulate_grad_area(&self, grad_a: &mut [f64], snap: &FeaSnapshot) {
        for e in 0..snap.elem_lengths.len() {
            grad_a[e] += self.weight * snap.densities[e] * snap.elem_lengths[e];
        }
    }

    fn has_position_gradient(&self) -> bool { true }

    fn accumulate_grad_pos(&self, _grad_pos: &mut [f64], _snap: &FeaSnapshot) {
        // d(rho * A * L)/d(node_pos) requires edge connectivity which is
        // not available in the snapshot. This gradient is handled in
        // fea_gradients.rs where we have access to the problem definition.
    }
}

// ─────────────────────────────────────────────────────────────
//  MaxStress: softplus penalty on |sigma| - threshold
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FeaMaxStress {
    pub weight: f64,
    pub edge_indices: Vec<usize>,
    pub threshold: Vec<f64>,
    pub sharpness: f64,
}

impl FeaObjectiveTrait for FeaMaxStress {
    fn loss(&self, snap: &FeaSnapshot) -> f64 {
        let mut loss = 0.0;
        for (k, &e) in self.edge_indices.iter().enumerate() {
            let abs_stress = snap.stresses[e].abs();
            loss += softplus(abs_stress, self.threshold[k], self.sharpness);
        }
        self.weight * loss
    }

    fn accumulate_grad_u(&self, grad_u: &mut [f64], snap: &FeaSnapshot) {
        // d softplus(|sigma|, threshold, k) / du
        // = softplus_grad(|sigma|, threshold, k) * d|sigma|/du
        // d|sigma|/du = sign(sigma) * dsigma/du
        // dsigma/du = E * d(strain)/du = E * d(c^T(u_j-u_i)/L)/du
        // = E/L * c[d] * (±1)
        // This needs element connectivity which we don't have in the snapshot.
        // The gradient module handles this via the adjoint.
        // For the adjoint: dJ/du is accumulated here as:
        // dJ/du_i[d] = weight * softplus_grad * sign(sigma) * E * c[d] * (-1/L)
        // dJ/du_j[d] = weight * softplus_grad * sign(sigma) * E * c[d] * (+1/L)
        // But we need edge_nodes and direction cosines from the problem.
        // These are handled in fea_gradients.rs.
        let _ = (grad_u, snap);
    }

    fn weight(&self) -> f64 { self.weight }
}

// ─────────────────────────────────────────────────────────────
//  TargetGeometry: L2 distance of deformed positions to targets
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FeaTargetGeometry {
    pub weight: f64,
    pub node_indices: Vec<usize>,
    pub targets: Vec<[f64; 3]>,
}

impl FeaObjectiveTrait for FeaTargetGeometry {
    fn loss(&self, snap: &FeaSnapshot) -> f64 {
        let mut loss = 0.0;
        for (k, &ni) in self.node_indices.iter().enumerate() {
            for d in 0..3 {
                let deformed = snap.deformed_xyz[[ni, d]];
                let diff = deformed - self.targets[k][d];
                loss += diff * diff;
            }
        }
        self.weight * loss
    }

    fn accumulate_grad_u(&self, grad_u: &mut [f64], snap: &FeaSnapshot) {
        // deformed = original + u, so d(deformed)/du = I
        // dJ/du_i[d] = 2 * weight * (deformed_i[d] - target_i[d])
        for (k, &ni) in self.node_indices.iter().enumerate() {
            for d in 0..3 {
                let deformed = snap.deformed_xyz[[ni, d]];
                let diff = deformed - self.targets[k][d];
                grad_u[ni * 3 + d] += self.weight * 2.0 * diff;
            }
        }
    }

    fn weight(&self) -> f64 { self.weight }

    fn has_position_gradient(&self) -> bool { true }

    fn accumulate_grad_pos(&self, grad_pos: &mut [f64], snap: &FeaSnapshot) {
        // deformed = original + u, so d(deformed)/d(original) = I
        // dJ/d(pos_i[d]) = 2 * weight * (deformed_i[d] - target_i[d])
        for (k, &ni) in self.node_indices.iter().enumerate() {
            for d in 0..3 {
                let deformed = snap.deformed_xyz[[ni, d]];
                let diff = deformed - self.targets[k][d];
                grad_pos[ni * 3 + d] += self.weight * 2.0 * diff;
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Bound penalties (reused from FDM objectives)
// ─────────────────────────────────────────────────────────────

pub fn fea_bounds_penalty(
    theta: &[f64], lb: &[f64], ub: &[f64],
    lb_idx: &[usize], ub_idx: &[usize], sharpness: f64,
) -> f64 {
    crate::objectives::bounds_penalty(theta, lb, ub, lb_idx, ub_idx, sharpness)
}

pub fn fea_bounds_penalty_grad(
    grad: &mut [f64], theta: &[f64], lb: &[f64], ub: &[f64],
    lb_idx: &[usize], ub_idx: &[usize], sharpness: f64, barrier_weight: f64,
) {
    crate::objectives::bounds_penalty_grad(grad, theta, lb, ub, lb_idx, ub_idx, sharpness, barrier_weight);
}

// ─────────────────────────────────────────────────────────────
//  Total loss
// ─────────────────────────────────────────────────────────────

pub fn fea_total_loss(objectives: &[Box<dyn FeaObjectiveTrait>], snap: &FeaSnapshot) -> f64 {
    objectives.iter().map(|obj| obj.loss(snap)).sum()
}
