//! Stiffness matrix assembly for 3D truss and beam elements.
//!
//! Truss: 2 nodes × 3 DOF = 6 DOF per element (axial only).
//! Beam:  2 nodes × 6 DOF = 12 DOF per element (axial + shear + bending + torsion).

use crate::fea_types::{BeamFormulation, FeaCache, FeaProblem};
use crate::types::TheseusError;

// ─────────────────────────────────────────────────────────────
//  Truss element stiffness (6×6)
// ─────────────────────────────────────────────────────────────

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
//  Beam local frame
// ─────────────────────────────────────────────────────────────

/// Build a 3×3 rotation matrix R whose rows are (e_x, e_y, e_z) in global coords.
/// `e_x` = beam axis (from node i to node j), `e_y`/`e_z` in the cross-section plane.
pub fn beam_local_frame(c: &[f64; 3]) -> [[f64; 3]; 3] {
    let ex = *c;

    // Choose reference vector not parallel to beam axis
    let v_ref = if ex[0].abs() < 0.9 && ex[1].abs() < 0.9 {
        [0.0, 0.0, 1.0]
    } else {
        [1.0, 0.0, 0.0]
    };

    // e_y = normalize(v_ref - (v_ref . e_x) * e_x)
    let dot = v_ref[0] * ex[0] + v_ref[1] * ex[1] + v_ref[2] * ex[2];
    let mut ey = [v_ref[0] - dot * ex[0], v_ref[1] - dot * ex[1], v_ref[2] - dot * ex[2]];
    let ey_len = (ey[0] * ey[0] + ey[1] * ey[1] + ey[2] * ey[2]).sqrt();
    ey[0] /= ey_len;
    ey[1] /= ey_len;
    ey[2] /= ey_len;

    // e_z = e_x × e_y
    let ez = [
        ex[1] * ey[2] - ex[2] * ey[1],
        ex[2] * ey[0] - ex[0] * ey[2],
        ex[0] * ey[1] - ex[1] * ey[0],
    ];

    [ex, ey, ez]
}

// ─────────────────────────────────────────────────────────────
//  Beam element stiffness (12×12)
// ─────────────────────────────────────────────────────────────

/// Compute the 12×12 element stiffness for an Euler-Bernoulli or Timoshenko beam
/// in local coordinates, then transform to global.
///
/// DOF order per node: [u, v, w, θx, θy, θz].
/// Euler-Bernoulli: phi_y = phi_z = 0.
/// Timoshenko: phi_y, phi_z from shear deformation.
pub fn beam_element_stiffness(
    e_mod: f64, nu: f64, sec: &crate::fea_types::FeaSection,
    length: f64, c: &[f64; 3], formulation: BeamFormulation,
) -> Vec<f64> {
    let a = sec.area;
    let iy = sec.iy;
    let iz = sec.iz;
    let j = sec.j;
    let l = length;
    let l2 = l * l;
    let l3 = l2 * l;
    let g = e_mod / (2.0 * (1.0 + nu));

    // Shear deformation parameters (0 for Euler-Bernoulli)
    let (phi_y, phi_z) = if formulation == BeamFormulation::Timoshenko {
        let asy = if sec.asy > 0.0 { sec.asy } else { (5.0 / 6.0) * a };
        let asz = if sec.asz > 0.0 { sec.asz } else { (5.0 / 6.0) * a };
        let py = if g * asy * l2 > 1e-30 { 12.0 * e_mod * iz / (g * asy * l2) } else { 0.0 };
        let pz = if g * asz * l2 > 1e-30 { 12.0 * e_mod * iy / (g * asz * l2) } else { 0.0 };
        (py, pz)
    } else {
        (0.0, 0.0)
    };

    // Build local 12×12 stiffness (symmetric, upper triangle then mirror)
    let mut kl = vec![0.0f64; 144];
    let set = |kl: &mut Vec<f64>, r: usize, c: usize, v: f64| {
        kl[r * 12 + c] = v;
        kl[c * 12 + r] = v;
    };

    let ea_l = e_mod * a / l;
    let gj_l = g * j / l;

    // Axial
    set(&mut kl, 0, 0, ea_l);
    set(&mut kl, 0, 6, -ea_l);
    set(&mut kl, 6, 6, ea_l);

    // Torsion
    set(&mut kl, 3, 3, gj_l);
    set(&mut kl, 3, 9, -gj_l);
    set(&mut kl, 9, 9, gj_l);

    // Bending in x-y plane (uses Iz, phi_y)
    {
        let d = 1.0 + phi_y;
        let a12 = 12.0 * e_mod * iz / (d * l3);
        let a6  = 6.0 * e_mod * iz / (d * l2);
        let a4  = (4.0 + phi_y) * e_mod * iz / (d * l);
        let a2  = (2.0 - phi_y) * e_mod * iz / (d * l);

        set(&mut kl, 1, 1, a12);
        set(&mut kl, 1, 5, a6);
        set(&mut kl, 1, 7, -a12);
        set(&mut kl, 1, 11, a6);
        set(&mut kl, 5, 5, a4);
        set(&mut kl, 5, 7, -a6);
        set(&mut kl, 5, 11, a2);
        set(&mut kl, 7, 7, a12);
        set(&mut kl, 7, 11, -a6);
        set(&mut kl, 11, 11, a4);
    }

    // Bending in x-z plane (uses Iy, phi_z)
    {
        let d = 1.0 + phi_z;
        let a12 = 12.0 * e_mod * iy / (d * l3);
        let a6  = 6.0 * e_mod * iy / (d * l2);
        let a4  = (4.0 + phi_z) * e_mod * iy / (d * l);
        let a2  = (2.0 - phi_z) * e_mod * iy / (d * l);

        set(&mut kl, 2, 2, a12);
        set(&mut kl, 2, 4, -a6);
        set(&mut kl, 2, 8, -a12);
        set(&mut kl, 2, 10, -a6);
        set(&mut kl, 4, 4, a4);
        set(&mut kl, 4, 8, a6);
        set(&mut kl, 4, 10, a2);
        set(&mut kl, 8, 8, a12);
        set(&mut kl, 8, 10, a6);
        set(&mut kl, 10, 10, a4);
    }

    // Transform to global: K_global = T^T * K_local * T
    // T = block_diag(R, R, R, R) where R is the 3×3 rotation matrix
    let rot = beam_local_frame(c);

    let mut kg = vec![0.0f64; 144];
    // temp = K_local * T
    let mut temp = vec![0.0f64; 144];
    for i in 0..12 {
        for j in 0..12 {
            let block_j = j / 3;
            let local_j = j % 3;
            let mut s = 0.0;
            for k in 0..12 {
                let block_k = k / 3;
                let local_k = k % 3;
                if block_k == block_j {
                    s += kl[i * 12 + k] * rot[local_j][local_k];
                }
            }
            temp[i * 12 + j] = s;
        }
    }

    // K_global = T^T * temp
    for i in 0..12 {
        for j in 0..12 {
            let block_i = i / 3;
            let local_i = i % 3;
            let mut s = 0.0;
            for k in 0..12 {
                let block_k = k / 3;
                let local_k = k % 3;
                if block_k == block_i {
                    s += rot[local_i][local_k] * temp[k * 12 + j];
                }
            }
            kg[i * 12 + j] = s;
        }
    }

    kg
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
//  Global stiffness assembly
// ─────────────────────────────────────────────────────────────

pub fn assemble_global_k(
    cache: &mut FeaCache,
    problem: &FeaProblem,
    areas: &[f64],
) -> Result<(), TheseusError> {
    for v in cache.k_matrix.values.iter_mut() {
        *v = 0.0;
    }

    let ne = problem.num_elements;
    let is_beam = problem.beam_formulation != BeamFormulation::Truss;

    for e in 0..ne {
        let props = &problem.element_props[e];
        let mat = &problem.materials[props.material_idx];
        let sec = &problem.sections[props.section_idx];
        let l = cache.elem_lengths[e];
        if l < f64::EPSILON {
            continue;
        }
        let c = &cache.elem_cos[e];

        if is_beam {
            let kg = beam_element_stiffness(
                mat.e, 0.3, sec, l, c, problem.beam_formulation,
            );
            for &(nz_idx, li, lj) in &cache.element_to_nz.entries[e] {
                cache.k_matrix.values[nz_idx] += kg[li * 12 + lj];
            }
        } else {
            let a = areas[props.section_idx];
            let ea_over_l = mat.e * a / l;
            let ke = bar_element_stiffness(ea_over_l, c);
            for &(nz_idx, li, lj) in &cache.element_to_nz.entries[e] {
                cache.k_matrix.values[nz_idx] += ke[li][lj];
            }
        }
    }

    Ok(())
}

// ─────────────────────────────────────────────────────────────
//  Load assembly
// ─────────────────────────────────────────────────────────────

pub fn assemble_loads(
    cache: &mut FeaCache,
    problem: &FeaProblem,
    areas: &[f64],
    _positions: &ndarray::Array2<f64>,
) {
    let n_free = cache.dof_map.num_free_dofs;
    let dpn = problem.dofs_per_node;

    for i in 0..n_free {
        cache.rhs[i] = 0.0;
    }

    for load in &problem.loads {
        for d in 0..3 {
            let global_dof = load.node_idx * dpn + d;
            if let Some(free_idx) = cache.dof_map.global_to_free[global_dof] {
                cache.rhs[free_idx] += load.force[d];
            }
        }
        if dpn == 6 {
            for d in 0..3 {
                let global_dof = load.node_idx * dpn + 3 + d;
                if let Some(free_idx) = cache.dof_map.global_to_free[global_dof] {
                    cache.rhs[free_idx] += load.moment[d];
                }
            }
        }
    }

    if problem.include_self_weight {
        let ne = problem.num_elements;
        for e in 0..ne {
            let props = &problem.element_props[e];
            let mat = &problem.materials[props.material_idx];
            let a = areas[props.section_idx];
            let l = cache.elem_lengths[e];
            let half_weight = mat.density * a * l * 0.5;

            let (ni, nj) = problem.edge_nodes[e];
            for d in 0..3 {
                let f = half_weight * problem.gravity[d];
                let gi = ni * dpn + d;
                let gj = nj * dpn + d;
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
