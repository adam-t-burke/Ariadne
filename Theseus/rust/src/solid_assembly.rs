//! Stiffness matrix assembly for 3D tet4 solid elements.
//!
//! Each tet4 element has 4 nodes × 3 DOF = 12 DOF total.
//! Linear shape functions give constant strain/stress per element.
//! Element stiffness: K_e = V * B^T * D * B  (12×12).

use crate::solid_types::{SolidCache, SolidProblem};
use crate::types::TheseusError;
use rayon::prelude::*;

// ─────────────────────────────────────────────────────────────
//  3×3 matrix helpers (stack-allocated)
// ─────────────────────────────────────────────────────────────

#[inline]
fn det3(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
  - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
  + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

#[inline]
fn inv3(m: &[[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let d = det3(m);
    if d.abs() < 1e-30 {
        return None;
    }
    let inv_d = 1.0 / d;
    Some([
        [
            (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv_d,
            (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv_d,
            (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv_d,
        ],
        [
            (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv_d,
            (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv_d,
            (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv_d,
        ],
        [
            (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv_d,
            (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv_d,
            (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv_d,
        ],
    ])
}

#[inline]
fn transpose3(m: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}

// ─────────────────────────────────────────────────────────────
//  Constitutive matrix D (6×6, isotropic linear elasticity)
// ─────────────────────────────────────────────────────────────

#[inline]
pub fn constitutive_matrix(e: f64, nu: f64) -> [[f64; 6]; 6] {
    let c = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
    let g = (1.0 - 2.0 * nu) / 2.0;
    [
        [c * (1.0 - nu), c * nu,          c * nu,          0.0,   0.0,   0.0  ],
        [c * nu,          c * (1.0 - nu), c * nu,          0.0,   0.0,   0.0  ],
        [c * nu,          c * nu,          c * (1.0 - nu), 0.0,   0.0,   0.0  ],
        [0.0,             0.0,             0.0,             c * g, 0.0,   0.0  ],
        [0.0,             0.0,             0.0,             0.0,   c * g, 0.0  ],
        [0.0,             0.0,             0.0,             0.0,   0.0,   c * g],
    ]
}

// ─────────────────────────────────────────────────────────────
//  Tet4 element: B matrix, volume, stiffness
// ─────────────────────────────────────────────────────────────

/// Compute the 6×12 strain-displacement matrix B and element volume for a tet4.
/// `coords` is [[x0,y0,z0], [x1,y1,z1], [x2,y2,z2], [x3,y3,z3]].
/// Returns (B, volume) or None if the element is degenerate.
#[inline]
pub fn tet4_b_matrix(coords: &[[f64; 3]; 4]) -> Option<([[f64; 12]; 6], f64)> {
    let j = [
        [coords[1][0] - coords[0][0], coords[2][0] - coords[0][0], coords[3][0] - coords[0][0]],
        [coords[1][1] - coords[0][1], coords[2][1] - coords[0][1], coords[3][1] - coords[0][1]],
        [coords[1][2] - coords[0][2], coords[2][2] - coords[0][2], coords[3][2] - coords[0][2]],
    ];

    let det_j = det3(&j);
    let volume = det_j.abs() / 6.0;
    if volume < 1e-30 {
        return None;
    }

    let j_inv_t = transpose3(&inv3(&j)?);

    // Shape function derivatives in parametric space:
    //   dN0/d(xi) = [-1, -1, -1]
    //   dN1/d(xi) = [ 1,  0,  0]
    //   dN2/d(xi) = [ 0,  1,  0]
    //   dN3/d(xi) = [ 0,  0,  1]
    //
    // Physical derivatives: dN/dx = J^{-T} * dN/d(xi)
    // J^{-T} is 3×3, dN/d(xi) is 3×1 per node → dN/dx is 3×1 per node
    let dn_dxi: [[f64; 3]; 4] = [
        [-1.0, -1.0, -1.0],
        [ 1.0,  0.0,  0.0],
        [ 0.0,  1.0,  0.0],
        [ 0.0,  0.0,  1.0],
    ];

    let mut dn_dx = [[0.0f64; 3]; 4]; // [node][x,y,z]
    for n in 0..4 {
        for i in 0..3 {
            dn_dx[n][i] = j_inv_t[i][0] * dn_dxi[n][0]
                        + j_inv_t[i][1] * dn_dxi[n][1]
                        + j_inv_t[i][2] * dn_dxi[n][2];
        }
    }

    // Build B (6×12)
    let mut b = [[0.0f64; 12]; 6];
    for n in 0..4 {
        let col = n * 3;
        let (dx, dy, dz) = (dn_dx[n][0], dn_dx[n][1], dn_dx[n][2]);

        b[0][col]     = dx;                       // eps_xx
        b[1][col + 1] = dy;                       // eps_yy
        b[2][col + 2] = dz;                       // eps_zz
        b[3][col]     = dy;  b[3][col + 1] = dx;  // gamma_xy
        b[4][col + 1] = dz;  b[4][col + 2] = dy;  // gamma_yz
        b[5][col]     = dz;  b[5][col + 2] = dx;  // gamma_xz
    }

    Some((b, volume))
}

/// Compute the 12×12 element stiffness matrix: K_e = V * B^T * D * B.
#[inline]
pub fn tet4_element_stiffness(
    b: &[[f64; 12]; 6],
    d: &[[f64; 6]; 6],
    volume: f64,
) -> [[f64; 12]; 12] {
    // DB = D * B  (6×12)
    let mut db = [[0.0f64; 12]; 6];
    for i in 0..6 {
        for j in 0..12 {
            let mut s = 0.0;
            for k in 0..6 {
                s += d[i][k] * b[k][j];
            }
            db[i][j] = s;
        }
    }

    // K_e = V * B^T * DB  (12×12)
    let mut ke = [[0.0f64; 12]; 12];
    for i in 0..12 {
        for j in 0..12 {
            let mut s = 0.0;
            for k in 0..6 {
                s += b[k][i] * db[k][j];
            }
            ke[i][j] = volume * s;
        }
    }
    ke
}

// ─────────────────────────────────────────────────────────────
//  Global stiffness assembly (parallel via rayon fold/reduce)
// ─────────────────────────────────────────────────────────────

/// Zero-allocation in-place update of K's values from tet4 element properties.
/// Uses rayon parallel fold/reduce with per-thread value buffers.
pub fn assemble_global_k(
    cache: &mut SolidCache,
    problem: &SolidProblem,
) -> Result<(), TheseusError> {
    let nnz = cache.k_matrix.values.len();
    let ne = problem.num_elements;

    let positions = &problem.node_positions;
    let materials = &problem.materials;
    let element_props = &problem.element_props;
    let elements = &problem.elements;
    let elem_to_nz = &cache.element_to_nz;

    let merged = (0..ne)
        .into_par_iter()
        .fold(
            || vec![0.0f64; nnz],
            |mut buf, e| {
                let nodes = &elements[e];
                let mat = &materials[element_props[e].material_idx];

                let coords = [
                    [positions[[nodes[0], 0]], positions[[nodes[0], 1]], positions[[nodes[0], 2]]],
                    [positions[[nodes[1], 0]], positions[[nodes[1], 1]], positions[[nodes[1], 2]]],
                    [positions[[nodes[2], 0]], positions[[nodes[2], 1]], positions[[nodes[2], 2]]],
                    [positions[[nodes[3], 0]], positions[[nodes[3], 1]], positions[[nodes[3], 2]]],
                ];

                if let Some((b, vol)) = tet4_b_matrix(&coords) {
                    let d = constitutive_matrix(mat.e, mat.nu);
                    let ke = tet4_element_stiffness(&b, &d, vol);

                    for &(nz_idx, li, lj) in &elem_to_nz.entries[e] {
                        buf[nz_idx] += ke[li][lj];
                    }
                }
                buf
            },
        )
        .reduce(
            || vec![0.0f64; nnz],
            |mut a, b| {
                for (ai, bi) in a.iter_mut().zip(b.iter()) {
                    *ai += *bi;
                }
                a
            },
        );

    cache.k_matrix.values.copy_from_slice(&merged);

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::tet4_b_matrix;

    #[test]
    fn skewed_tet_gradients_match_closed_form() {
        let coords = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.0, 3.0, 0.0],
            [0.0, 0.0, 4.0],
        ];

        let (b, volume) = tet4_b_matrix(&coords).expect("valid tet");
        assert!((volume - 4.0).abs() < 1e-12);

        let expected = [
            [-0.5, -1.0 / 6.0, -0.25],
            [ 0.5, -1.0 / 6.0,  0.0],
            [ 0.0,  1.0 / 3.0,  0.0],
            [ 0.0,  0.0,        0.25],
        ];

        for n in 0..4 {
            let col = n * 3;
            assert!((b[0][col] - expected[n][0]).abs() < 1e-12);
            assert!((b[1][col + 1] - expected[n][1]).abs() < 1e-12);
            assert!((b[2][col + 2] - expected[n][2]).abs() < 1e-12);
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Load assembly
// ─────────────────────────────────────────────────────────────

/// Assemble the RHS force vector (free DOFs only).
pub fn assemble_loads(
    cache: &mut SolidCache,
    problem: &SolidProblem,
) {
    let n_free = cache.dof_map.num_free_dofs;
    for i in 0..n_free {
        cache.rhs[i] = 0.0;
    }

    for load in &problem.loads {
        for d in 0..3 {
            let global_dof = load.node_idx * 3 + d;
            if let Some(free_idx) = cache.dof_map.global_to_free[global_dof] {
                cache.rhs[free_idx] += load.force[d];
            }
        }
    }

    if problem.include_self_weight {
        let ne = problem.num_elements;
        let positions = &problem.node_positions;
        for e in 0..ne {
            let nodes = &problem.elements[e];
            let mat = &problem.materials[problem.element_props[e].material_idx];

            let coords = [
                [positions[[nodes[0], 0]], positions[[nodes[0], 1]], positions[[nodes[0], 2]]],
                [positions[[nodes[1], 0]], positions[[nodes[1], 1]], positions[[nodes[1], 2]]],
                [positions[[nodes[2], 0]], positions[[nodes[2], 1]], positions[[nodes[2], 2]]],
                [positions[[nodes[3], 0]], positions[[nodes[3], 1]], positions[[nodes[3], 2]]],
            ];

            if let Some((_, vol)) = tet4_b_matrix(&coords) {
                let weight = mat.density * vol;
                let quarter_weight = weight * 0.25;

                for &ni in nodes {
                    for d in 0..3 {
                        let f = quarter_weight * problem.gravity[d];
                        let gi = ni * 3 + d;
                        if let Some(fi) = cache.dof_map.global_to_free[gi] {
                            cache.rhs[fi] += f;
                        }
                    }
                }
            }
        }
    }
}
