//! Stiffness matrix assembly for 3D tet4 solid elements with 4-point
//! Gauss quadrature.
//!
//! Each tet4 element has 4 nodes × 3 DOF = 12 DOF total.
//! Element stiffness: K_e = Σ w_i |J(ξ_i)| B(ξ_i)^T D B(ξ_i)  (12×12).

use crate::solid_types::{SolidCache, SolidProblem};
use crate::types::TheseusError;
use ndarray::Array2;
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
//  Coordinate extraction helper
// ─────────────────────────────────────────────────────────────

#[inline]
pub fn extract_tet4_coords(positions: &Array2<f64>, nodes: &[usize]) -> [[f64; 3]; 4] {
    [
        [positions[[nodes[0], 0]], positions[[nodes[0], 1]], positions[[nodes[0], 2]]],
        [positions[[nodes[1], 0]], positions[[nodes[1], 1]], positions[[nodes[1], 2]]],
        [positions[[nodes[2], 0]], positions[[nodes[2], 1]], positions[[nodes[2], 2]]],
        [positions[[nodes[3], 0]], positions[[nodes[3], 1]], positions[[nodes[3], 2]]],
    ]
}

#[inline]
pub fn extract_tet10_coords(positions: &Array2<f64>, nodes: &[usize]) -> [[f64; 3]; 10] {
    let mut coords = [[0.0; 3]; 10];
    for i in 0..10 {
        coords[i] = [positions[[nodes[i], 0]], positions[[nodes[i], 1]], positions[[nodes[i], 2]]];
    }
    coords
}

// ─────────────────────────────────────────────────────────────
//  4-point Gauss quadrature on reference tetrahedron
// ─────────────────────────────────────────────────────────────

const GP_A: f64 = 0.138_196_601_125_010_5;  // (5 - sqrt(5)) / 20
const GP_B: f64 = 0.585_410_196_624_968_5;  // (5 + 3*sqrt(5)) / 20
const GP_W: f64 = 1.0 / 24.0;              // each weight (sums to 1/6 = ref tet volume)

pub const NUM_GP: usize = 4;

pub const TET4_GAUSS_POINTS: [([f64; 3], f64); NUM_GP] = [
    ([GP_A, GP_A, GP_A], GP_W),
    ([GP_B, GP_A, GP_A], GP_W),
    ([GP_A, GP_B, GP_A], GP_W),
    ([GP_A, GP_A, GP_B], GP_W),
];

// ─────────────────────────────────────────────────────────────
//  Shape functions and derivatives
// ─────────────────────────────────────────────────────────────

#[inline]
pub fn tet4_shape_functions(xi: f64, eta: f64, zeta: f64) -> [f64; 4] {
    [1.0 - xi - eta - zeta, xi, eta, zeta]
}

/// Shape function derivatives in parametric space (constant for linear tet4).
#[inline]
pub fn tet4_shape_derivs(_xi: f64, _eta: f64, _zeta: f64) -> [[f64; 3]; 4] {
    [
        [-1.0, -1.0, -1.0],
        [ 1.0,  0.0,  0.0],
        [ 0.0,  1.0,  0.0],
        [ 0.0,  0.0,  1.0],
    ]
}

#[inline]
pub fn tet10_shape_functions(xi: f64, eta: f64, zeta: f64) -> [f64; 10] {
    let l1 = 1.0 - xi - eta - zeta;
    let l2 = xi;
    let l3 = eta;
    let l4 = zeta;
    [
        l1 * (2.0 * l1 - 1.0),
        l2 * (2.0 * l2 - 1.0),
        l3 * (2.0 * l3 - 1.0),
        l4 * (2.0 * l4 - 1.0),
        4.0 * l1 * l2,
        4.0 * l2 * l3,
        4.0 * l3 * l1,
        4.0 * l1 * l4,
        4.0 * l2 * l4,
        4.0 * l3 * l4,
    ]
}

#[inline]
pub fn tet10_shape_derivs(xi: f64, eta: f64, zeta: f64) -> [[f64; 3]; 10] {
    let l1 = 1.0 - xi - eta - zeta;
    let l2 = xi;
    let l3 = eta;
    let l4 = zeta;
    
    let mut dn = [[0.0; 3]; 10];
    
    let dn0_dl1 = 4.0 * l1 - 1.0;
    dn[0][0] = -dn0_dl1; dn[0][1] = -dn0_dl1; dn[0][2] = -dn0_dl1;
    
    let dn1_dl2 = 4.0 * l2 - 1.0;
    dn[1][0] = dn1_dl2;
    
    let dn2_dl3 = 4.0 * l3 - 1.0;
    dn[2][1] = dn2_dl3;
    
    let dn3_dl4 = 4.0 * l4 - 1.0;
    dn[3][2] = dn3_dl4;
    
    dn[4][0] = 4.0 * (l1 - l2); dn[4][1] = -4.0 * l2; dn[4][2] = -4.0 * l2;
    dn[5][0] = 4.0 * l3;        dn[5][1] = 4.0 * l2;
    dn[6][0] = -4.0 * l3;       dn[6][1] = 4.0 * (l1 - l3); dn[6][2] = -4.0 * l3;
    dn[7][0] = -4.0 * l4;       dn[7][1] = -4.0 * l4;       dn[7][2] = 4.0 * (l1 - l4);
    dn[8][0] = 4.0 * l4;        dn[8][2] = 4.0 * l2;
    dn[9][1] = 4.0 * l4;        dn[9][2] = 4.0 * l3;
    
    dn
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

/// Precompute D for each material so we don't rebuild it per-element.
pub fn precompute_d_matrices(problem: &SolidProblem) -> Vec<[[f64; 6]; 6]> {
    problem.materials.iter()
        .map(|m| constitutive_matrix(m.e, m.nu))
        .collect()
}

// ─────────────────────────────────────────────────────────────
//  D-sparsity-aware D*B multiply
// ─────────────────────────────────────────────────────────────

/// Compute DB = D * B (6×12) exploiting the block structure of the
/// isotropic elasticity D: dense 3×3 upper-left, diagonal 3×3 lower-right.
#[inline]
pub fn d_times_b(d: &[[f64; 6]; 6], b: &[[f64; 12]; 6]) -> [[f64; 12]; 6] {
    let mut db = [[0.0f64; 12]; 6];
    for j in 0..12 {
        for i in 0..3 {
            db[i][j] = d[i][0] * b[0][j] + d[i][1] * b[1][j] + d[i][2] * b[2][j];
        }
        db[3][j] = d[3][3] * b[3][j];
        db[4][j] = d[4][4] * b[4][j];
        db[5][j] = d[5][5] * b[5][j];
    }
    db
}

#[inline]
pub fn d_times_b_30(d: &[[f64; 6]; 6], b: &[[f64; 30]; 6]) -> [[f64; 30]; 6] {
    let mut db = [[0.0f64; 30]; 6];
    for j in 0..30 {
        for i in 0..3 {
            db[i][j] = d[i][0] * b[0][j] + d[i][1] * b[1][j] + d[i][2] * b[2][j];
        }
        db[3][j] = d[3][3] * b[3][j];
        db[4][j] = d[4][4] * b[4][j];
        db[5][j] = d[5][5] * b[5][j];
    }
    db
}

// ─────────────────────────────────────────────────────────────
//  Tet4 element: B matrix at a parametric point
// ─────────────────────────────────────────────────────────────

/// Compute the 6×12 strain-displacement matrix B and |det(J)| at parametric
/// point (xi, eta, zeta). Returns None if the element is degenerate.
#[inline]
pub fn tet4_b_matrix_at(
    coords: &[[f64; 3]; 4],
    xi: f64, eta: f64, zeta: f64,
) -> Option<([[f64; 12]; 6], f64)> {
    let j = [
        [coords[1][0] - coords[0][0], coords[2][0] - coords[0][0], coords[3][0] - coords[0][0]],
        [coords[1][1] - coords[0][1], coords[2][1] - coords[0][1], coords[3][1] - coords[0][1]],
        [coords[1][2] - coords[0][2], coords[2][2] - coords[0][2], coords[3][2] - coords[0][2]],
    ];

    let det_j = det3(&j);
    if det_j.abs() < 1e-30 {
        return None;
    }

    let j_inv_t = transpose3(&inv3(&j)?);
    let dn_dxi = tet4_shape_derivs(xi, eta, zeta);

    let mut dn_dx = [[0.0f64; 3]; 4];
    for n in 0..4 {
        for i in 0..3 {
            dn_dx[n][i] = j_inv_t[i][0] * dn_dxi[n][0]
                        + j_inv_t[i][1] * dn_dxi[n][1]
                        + j_inv_t[i][2] * dn_dxi[n][2];
        }
    }

    let mut b = [[0.0f64; 12]; 6];
    for n in 0..4 {
        let col = n * 3;
        let (dx, dy, dz) = (dn_dx[n][0], dn_dx[n][1], dn_dx[n][2]);

        b[0][col]     = dx;
        b[1][col + 1] = dy;
        b[2][col + 2] = dz;
        b[3][col]     = dy;  b[3][col + 1] = dx;
        b[4][col + 1] = dz;  b[4][col + 2] = dy;
        b[5][col]     = dz;  b[5][col + 2] = dx;
    }

    Some((b, det_j.abs()))
}

#[inline]
pub fn tet10_b_matrix_at(
    coords: &[[f64; 3]; 10],
    xi: f64, eta: f64, zeta: f64,
) -> Option<([[f64; 30]; 6], f64)> {
    let dn_dxi = tet10_shape_derivs(xi, eta, zeta);
    
    let mut j = [[0.0f64; 3]; 3];
    for n in 0..10 {
        for r in 0..3 {
            for c in 0..3 {
                j[r][c] += coords[n][r] * dn_dxi[n][c];
            }
        }
    }

    let det_j = det3(&j);
    if det_j.abs() < 1e-30 {
        return None;
    }

    let j_inv_t = transpose3(&inv3(&j)?);

    let mut dn_dx = [[0.0f64; 3]; 10];
    for n in 0..10 {
        for i in 0..3 {
            dn_dx[n][i] = j_inv_t[i][0] * dn_dxi[n][0]
                        + j_inv_t[i][1] * dn_dxi[n][1]
                        + j_inv_t[i][2] * dn_dxi[n][2];
        }
    }

    let mut b = [[0.0f64; 30]; 6];
    for n in 0..10 {
        let col = n * 3;
        let (dx, dy, dz) = (dn_dx[n][0], dn_dx[n][1], dn_dx[n][2]);

        b[0][col]     = dx;
        b[1][col + 1] = dy;
        b[2][col + 2] = dz;
        b[3][col]     = dy;  b[3][col + 1] = dx;
        b[4][col + 1] = dz;  b[4][col + 2] = dy;
        b[5][col]     = dz;  b[5][col + 2] = dx;
    }

    Some((b, det_j.abs()))
}

/// Backward-compatible wrapper: B matrix and volume (single evaluation).
#[inline]
pub fn tet4_b_matrix(coords: &[[f64; 3]; 4]) -> Option<([[f64; 12]; 6], f64)> {
    let (b, abs_det_j) = tet4_b_matrix_at(coords, 0.25, 0.25, 0.25)?;
    Some((b, abs_det_j / 6.0))
}

// ─────────────────────────────────────────────────────────────
//  Tet4 element stiffness via 4-point quadrature
// ─────────────────────────────────────────────────────────────

/// Compute the 12×12 element stiffness matrix using 4-point Gauss quadrature.
/// Exploits K_e symmetry (upper triangle + mirror) and D block sparsity.
#[inline]
pub fn tet4_element_stiffness_quad(
    coords: &[[f64; 3]; 4],
    d: &[[f64; 6]; 6],
) -> Option<[[f64; 12]; 12]> {
    let mut ke = [[0.0f64; 12]; 12];

    for &(pt, w) in &TET4_GAUSS_POINTS {
        let (b, abs_det_j) = tet4_b_matrix_at(coords, pt[0], pt[1], pt[2])?;
        let db = d_times_b(d, &b);
        let factor = w * abs_det_j;

        for i in 0..12 {
            for j in i..12 {
                let mut s = 0.0;
                for k in 0..6 {
                    s += b[k][i] * db[k][j];
                }
                ke[i][j] += factor * s;
            }
        }
    }

    for i in 0..12 {
        for j in 0..i {
            ke[i][j] = ke[j][i];
        }
    }

    Some(ke)
}

#[inline]
pub fn tet10_element_stiffness_quad(
    coords: &[[f64; 3]; 10],
    d: &[[f64; 6]; 6],
) -> Option<[[f64; 30]; 30]> {
    let mut ke = [[0.0f64; 30]; 30];

    for &(pt, w) in &TET4_GAUSS_POINTS {
        let (b, abs_det_j) = tet10_b_matrix_at(coords, pt[0], pt[1], pt[2])?;
        let db = d_times_b_30(d, &b);
        let factor = w * abs_det_j;

        for i in 0..30 {
            for j in i..30 {
                let mut s = 0.0;
                for k in 0..6 {
                    s += b[k][i] * db[k][j];
                }
                ke[i][j] += factor * s;
            }
        }
    }

    for i in 0..30 {
        for j in 0..i {
            ke[i][j] = ke[j][i];
        }
    }

    Some(ke)
}

/// Legacy single-point stiffness (kept for tests).
#[inline]
pub fn tet4_element_stiffness(
    b: &[[f64; 12]; 6],
    d: &[[f64; 6]; 6],
    volume: f64,
) -> [[f64; 12]; 12] {
    let db = d_times_b(d, b);
    let mut ke = [[0.0f64; 12]; 12];
    for i in 0..12 {
        for j in i..12 {
            let mut s = 0.0;
            for k in 0..6 {
                s += b[k][i] * db[k][j];
            }
            ke[i][j] = volume * s;
        }
    }
    for i in 0..12 {
        for j in 0..i {
            ke[i][j] = ke[j][i];
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
    let element_props = &problem.element_props;
    let elements = &problem.elements;
    let elem_to_nz = &cache.element_to_nz;

    let d_matrices = precompute_d_matrices(problem);

    let merged = (0..ne)
        .into_par_iter()
        .fold(
            || vec![0.0f64; nnz],
            |mut buf, e| {
                let d = &d_matrices[element_props[e].material_idx];
                let nodes = &elements[e];

                if nodes.len() == 10 {
                    let coords = extract_tet10_coords(positions, nodes);
                    if let Some(ke) = tet10_element_stiffness_quad(&coords, d) {
                        for &(nz_idx, li, lj) in &elem_to_nz.entries[e] {
                            buf[nz_idx] += ke[li][lj];
                        }
                    }
                } else {
                    let coords = extract_tet4_coords(positions, nodes);
                    if let Some(ke) = tet4_element_stiffness_quad(&coords, d) {
                        for &(nz_idx, li, lj) in &elem_to_nz.entries[e] {
                            buf[nz_idx] += ke[li][lj];
                        }
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

            if nodes.len() == 10 {
                let coords = extract_tet10_coords(positions, nodes);
                for &(pt, w) in &TET4_GAUSS_POINTS {
                    let n_vals = tet10_shape_functions(pt[0], pt[1], pt[2]);
                    if let Some((_, abs_det_j)) = tet10_b_matrix_at(&coords, pt[0], pt[1], pt[2]) {
                        let factor = w * abs_det_j * mat.density;
                        for node_local in 0..10 {
                            let ni = nodes[node_local];
                            for d in 0..3 {
                                let gi = ni * 3 + d;
                                if let Some(fi) = cache.dof_map.global_to_free[gi] {
                                    cache.rhs[fi] += factor * n_vals[node_local] * problem.gravity[d];
                                }
                            }
                        }
                    }
                }
            } else {
                let coords = extract_tet4_coords(positions, nodes);
                for &(pt, w) in &TET4_GAUSS_POINTS {
                    let n_vals = tet4_shape_functions(pt[0], pt[1], pt[2]);
                    if let Some((_, abs_det_j)) = tet4_b_matrix_at(&coords, pt[0], pt[1], pt[2]) {
                        let factor = w * abs_det_j * mat.density;
                        for node_local in 0..4 {
                            let ni = nodes[node_local];
                            for d in 0..3 {
                                let gi = ni * 3 + d;
                                if let Some(fi) = cache.dof_map.global_to_free[gi] {
                                    cache.rhs[fi] += factor * n_vals[node_local] * problem.gravity[d];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Tests
// ─────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

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

    #[test]
    fn quad_stiffness_matches_analytical_for_tet4() {
        let coords = [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [1.0, 3.0, 0.0],
            [0.0, 0.0, 4.0],
        ];

        let d = constitutive_matrix(200e9, 0.3);

        let (b, vol) = tet4_b_matrix(&coords).expect("valid tet");
        let ke_analytical = tet4_element_stiffness(&b, &d, vol);
        let ke_quad = tet4_element_stiffness_quad(&coords, &d).expect("valid tet");

        for i in 0..12 {
            for j in 0..12 {
                assert!(
                    (ke_analytical[i][j] - ke_quad[i][j]).abs() < 1e-6 * ke_analytical[i][j].abs().max(1.0),
                    "mismatch at [{i}][{j}]: analytical={}, quad={}",
                    ke_analytical[i][j], ke_quad[i][j]
                );
            }
        }
    }

    #[test]
    fn quad_stiffness_is_symmetric() {
        let coords = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ];
        let d = constitutive_matrix(200e9, 0.3);
        let ke = tet4_element_stiffness_quad(&coords, &d).expect("valid tet");

        for i in 0..12 {
            for j in 0..12 {
                assert!(
                    (ke[i][j] - ke[j][i]).abs() < 1e-10,
                    "not symmetric at [{i}][{j}]: {} vs {}",
                    ke[i][j], ke[j][i]
                );
            }
        }
    }

    #[test]
    fn shape_functions_partition_of_unity() {
        for &(pt, _) in &TET4_GAUSS_POINTS {
            let n = tet4_shape_functions(pt[0], pt[1], pt[2]);
            let sum: f64 = n.iter().sum();
            assert!((sum - 1.0).abs() < 1e-15, "N sum = {sum} at {pt:?}");
        }
    }
}
