//! FEA objective functions for optimization.
//!
//! Each objective implements FeaObjectiveTrait with loss and gradient methods.
//! Gradients are accumulated into grad_u (dJ/du for the adjoint) and optionally
//! into grad_area or grad_pos for explicit design-variable contributions.

use crate::fea_types::{FeaObjectiveTrait, FeaSnapshot};
use crate::objectives::{softplus, softplus_grad};

// ─────────────────────────────────────────────────────────────
//  Compliance: minimize strain energy u^T K u = u^T f
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FeaCompliance {
    pub weight: f64,
}

impl FeaObjectiveTrait for FeaCompliance {
    fn loss(&self, snap: &FeaSnapshot) -> f64 {
        // J = u^T K u = Σ_e (EA/L) * δ_e^2  where δ_e = c^T(u_j - u_i)
        // Equivalently: Σ_e force_e * strain_e * L_e = Σ_e (EA*strain^2*L)
        let mut compliance = 0.0;
        for e in 0..snap.elem_lengths.len() {
            let l = snap.elem_lengths[e];
            if l < f64::EPSILON { continue; }
            compliance += snap.axial_forces[e] * snap.strains[e] * l;
        }
        self.weight * compliance
    }

    fn accumulate_grad_u(&self, _grad_u: &mut [f64], _snap: &FeaSnapshot) {
        // For compliance J = u^T K u, dJ/du = 2Ku = 2f.
        // This requires the RHS vector f which is in cache.rhs (free-DOF space).
        // The gradient module handles this via the is_compliance() flag,
        // injecting 2*cache.rhs into grad_u_full directly.
    }

    fn weight(&self) -> f64 { self.weight }
    fn is_compliance(&self) -> bool { true }

    fn has_area_gradient(&self) -> bool { true }

    fn accumulate_grad_area(&self, _grad_a: &mut [f64], _snap: &FeaSnapshot) {
        // For compliance, the area gradient is fully handled by the adjoint:
        // dJ/dA_s = -λ^T (dK/dA_s) u, computed in fea_gradients.rs.
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
            let w_k = exp_vals[k] / sum_exp;
            for d in 0..3 {
                let u_d = snap.displacements[ni * 3 + d];
                // d(smooth_max)/du_d = w_k * u_d / mag
                grad_u[ni * 3 + d] += self.weight * w_k * u_d / mag;
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
            // areas is per-section; look up via section_indices
            let a = snap.areas[snap.section_indices[e]];
            total += snap.densities[e] * a * snap.elem_lengths[e];
        }
        self.weight * total
    }

    fn accumulate_grad_u(&self, _grad_u: &mut [f64], _snap: &FeaSnapshot) {
        // Weight does not depend on displacements directly.
    }

    fn weight(&self) -> f64 { self.weight }

    fn has_area_gradient(&self) -> bool { true }

    fn accumulate_grad_area(&self, grad_a: &mut [f64], snap: &FeaSnapshot) {
        // dJ/dA_s = Σ_{e : section_indices[e]==s} weight * ρ_e * L_e
        for e in 0..snap.elem_lengths.len() {
            let s = snap.section_indices[e];
            grad_a[s] += self.weight * snap.densities[e] * snap.elem_lengths[e];
        }
    }

    fn has_position_gradient(&self) -> bool { true }

    fn accumulate_grad_pos(&self, grad_pos: &mut [f64], snap: &FeaSnapshot) {
        // dJ/d(pos) = Σ_e weight * ρ_e * A_e * dL_e/d(pos)
        // dL/d(xi_d) = -c[d],  dL/d(xj_d) = +c[d]
        for e in 0..snap.elem_lengths.len() {
            let l = snap.elem_lengths[e];
            if l < f64::EPSILON { continue; }
            let (ni, nj) = snap.edge_nodes[e];
            let c = snap.elem_cos[e];
            let a = snap.areas[snap.section_indices[e]];
            let rho_a_w = self.weight * snap.densities[e] * a;
            for d in 0..3 {
                grad_pos[ni * 3 + d] += rho_a_w * (-c[d]);
                grad_pos[nj * 3 + d] += rho_a_w * c[d];
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  MaxStress: softplus penalty on |σ| - threshold
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
        // σ_e = E_e * c_e^T (u_j - u_i) / L_e
        // d|σ|/du = sign(σ) * dσ/du
        // dσ/du_i[d] = -E/L * c[d],  dσ/du_j[d] = +E/L * c[d]
        // dJ/du = weight * softplus_grad(|σ|, threshold, k) * sign(σ) * dσ/du
        for (k, &e) in self.edge_indices.iter().enumerate() {
            let l = snap.elem_lengths[e];
            if l < f64::EPSILON { continue; }
            let sigma = snap.stresses[e];
            let abs_sigma = sigma.abs();
            let sign_sigma = if sigma >= 0.0 { 1.0 } else { -1.0 };
            let sp_grad = softplus_grad(abs_sigma, self.threshold[k], self.sharpness);

            let (ni, nj) = snap.edge_nodes[e];
            let c = snap.elem_cos[e];
            let e_mod = snap.materials_e[e];
            let e_over_l = e_mod / l;

            let factor = self.weight * sp_grad * sign_sigma * e_over_l;
            for d in 0..3 {
                grad_u[ni * 3 + d] += factor * (-c[d]);
                grad_u[nj * 3 + d] += factor * c[d];
            }
        }
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
                let diff = snap.deformed_xyz[[ni, d]] - self.targets[k][d];
                loss += diff * diff;
            }
        }
        self.weight * loss
    }

    fn accumulate_grad_u(&self, grad_u: &mut [f64], snap: &FeaSnapshot) {
        // deformed = original + u  =>  d(deformed)/du = I
        for (k, &ni) in self.node_indices.iter().enumerate() {
            for d in 0..3 {
                let diff = snap.deformed_xyz[[ni, d]] - self.targets[k][d];
                grad_u[ni * 3 + d] += self.weight * 2.0 * diff;
            }
        }
    }

    fn weight(&self) -> f64 { self.weight }

    fn has_position_gradient(&self) -> bool { true }

    fn accumulate_grad_pos(&self, grad_pos: &mut [f64], snap: &FeaSnapshot) {
        // deformed = original + u  =>  d(deformed)/d(original) = I
        for (k, &ni) in self.node_indices.iter().enumerate() {
            for d in 0..3 {
                let diff = snap.deformed_xyz[[ni, d]] - self.targets[k][d];
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
