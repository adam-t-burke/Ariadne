//! Stiffness matrix assembly for 3D bar/truss elements.
//!
//! Each bar element has 2 nodes × 3 DOF = 6 DOF total.
//! Element stiffness in global coordinates:
//!
//!   ke = (EA/L) * [  T  -T ]
//!                  [ -T   T ]
//!
//! where T = c * c^T, c = [cx, cy, cz] (direction cosines).

use crate::fea_types::{FeaCache, FeaProblem};
use crate::types::TheseusError;

// ─────────────────────────────────────────────────────────────
//  Element stiffness
// ─────────────────────────────────────────────────────────────

/// Compute the 6×6 element stiffness matrix for a 3D bar in global coordinates.
/// Returns row-major [6][6].
#[inline]
pub fn bar_element_stiffness(ea_over_l: f64, c: &[f64; 3]) -> [[f64; 6]; 6] {
    let mut ke = [[0.0f64; 6]; 6];
    for i in 0..3 {
        for j in 0..3 {
            let t = ea_over_l * c[i] * c[j];
            ke[i][j] = t;
            ke[i][j + 3] = -t;
            ke[i + 3][j] = -t;
            ke[i + 3][j + 3] = t;
        }
    }
    ke
}

// ─────────────────────────────────────────────────────────────
//  Update element geometry from current node positions
// ─────────────────────────────────────────────────────────────

pub fn update_element_geometry(cache: &mut FeaCache, problem: &FeaProblem, positions: &ndarray::Array2<f64>) {
    let ne = problem.num_elements;
    for e in 0..ne {
        let (ni, nj) = problem.edge_nodes[e];
        let dx = positions[[nj, 0]] - positions[[ni, 0]];
        let dy = positions[[nj, 1]] - positions[[ni, 1]];
        let dz = positions[[nj, 2]] - positions[[ni, 2]];
        let l = (dx * dx + dy * dy + dz * dz).max(0.0).sqrt();
        cache.elem_lengths[e] = l;
        if l > f64::EPSILON {
            cache.elem_cos[e] = [dx / l, dy / l, dz / l];
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Global stiffness assembly (in-place via element_to_nz)
// ─────────────────────────────────────────────────────────────

/// Zero-allocation in-place update of K's values from current element properties.
pub fn assemble_global_k(
    cache: &mut FeaCache,
    problem: &FeaProblem,
    areas: &[f64],
) -> Result<(), TheseusError> {
    for v in cache.k_matrix.values.iter_mut() {
        *v = 0.0;
    }

    let ne = problem.num_elements;
    for e in 0..ne {
        let props = &problem.element_props[e];
        let mat = &problem.materials[props.material_idx];
        let a = areas[props.section_idx];
        let l = cache.elem_lengths[e];
        if l < f64::EPSILON {
            continue;
        }
        let ea_over_l = mat.e * a / l;
        let c = &cache.elem_cos[e];
        let ke = bar_element_stiffness(ea_over_l, c);

        for &(nz_idx, li, lj) in &cache.element_to_nz.entries[e] {
            cache.k_matrix.values[nz_idx] += ke[li][lj];
        }
    }

    Ok(())
}

// ─────────────────────────────────────────────────────────────
//  Load assembly
// ─────────────────────────────────────────────────────────────

/// Assemble the RHS force vector (free DOFs only).
pub fn assemble_loads(
    cache: &mut FeaCache,
    problem: &FeaProblem,
    areas: &[f64],
    _positions: &ndarray::Array2<f64>,
) {
    let n_free = cache.dof_map.num_free_dofs;
    for i in 0..n_free {
        cache.rhs[i] = 0.0;
    }

    // Point loads
    for load in &problem.loads {
        for d in 0..3 {
            let global_dof = load.node_idx * 3 + d;
            if let Some(free_idx) = cache.dof_map.global_to_free[global_dof] {
                cache.rhs[free_idx] += load.force[d];
            }
        }
    }

    // Self-weight (gravity loads distributed to nodes)
    if problem.include_self_weight {
        let ne = problem.num_elements;
        for e in 0..ne {
            let props = &problem.element_props[e];
            let mat = &problem.materials[props.material_idx];
            let a = areas[props.section_idx];
            let l = cache.elem_lengths[e];
            let weight = mat.density * a * l;
            let half_weight = weight * 0.5;

            let (ni, nj) = problem.edge_nodes[e];
            for d in 0..3 {
                let f = half_weight * problem.gravity[d];
                let gi = ni * 3 + d;
                let gj = nj * 3 + d;
                if let Some(fi) = cache.dof_map.global_to_free[gi] {
                    cache.rhs[fi] += f;
                }
                if let Some(fj) = cache.dof_map.global_to_free[gj] {
                    cache.rhs[fj] += f;
                }
            }
        }
    }
}
