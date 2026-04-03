//! Stiffness matrix assembly for flat shell triangle elements with 6 DOFs per node.
//!
//! Element formulation: CST membrane + DKT bending + Hughes-Brezzi drilling.
//! Reference: Batoz, Bathe, Ho (1980) for DKT; Hughes & Brezzi (1989) for drilling.

use crate::shell_types::{ShellCache, ShellProblem};
use crate::types::TheseusError;
use rayon::prelude::*;

/// 3-point mid-side Hammer quadrature on the reference triangle.
/// Points are at edge midpoints; weights are 1/3 each.
/// Exact for quadratic polynomials (sufficient for DKT since B is linear in xi,eta).
const TRI_GAUSS_PTS: [([f64; 2], f64); 3] = [
    ([0.5, 0.0], 1.0 / 3.0),
    ([0.5, 0.5], 1.0 / 3.0),
    ([0.0, 0.5], 1.0 / 3.0),
];

// ─────────────────────────────────────────────────────────────
//  Local coordinate system for a triangle in 3D
// ─────────────────────────────────────────────────────────────

/// Returns (ex, ey, ez, local_coords, area) for a triangle defined by 3 points in 3D.
/// ex is along edge 1→2, ez is the unit normal, ey = ez × ex.
/// local_coords are (x1,y1), (x2,y2), (x3,y3) in the local frame (x1=y1=y2=0).
pub fn tri_local_frame(
    p1: &[f64], p2: &[f64], p3: &[f64],
) -> ([f64; 3], [f64; 3], [f64; 3], [f64; 6], f64) {
    let v12 = [p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]];
    let v13 = [p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]];

    let l12 = (v12[0] * v12[0] + v12[1] * v12[1] + v12[2] * v12[2]).sqrt();
    let ex = [v12[0] / l12, v12[1] / l12, v12[2] / l12];

    let mut ez = [
        ex[1] * v13[2] - ex[2] * v13[1],
        ex[2] * v13[0] - ex[0] * v13[2],
        ex[0] * v13[1] - ex[1] * v13[0],
    ];
    let lz = (ez[0] * ez[0] + ez[1] * ez[1] + ez[2] * ez[2]).sqrt();
    ez[0] /= lz;
    ez[1] /= lz;
    ez[2] /= lz;

    let ey = [
        ez[1] * ex[2] - ez[2] * ex[1],
        ez[2] * ex[0] - ez[0] * ex[2],
        ez[0] * ex[1] - ez[1] * ex[0],
    ];

    let x1 = 0.0;
    let y1 = 0.0;
    let x2 = l12;
    let y2 = 0.0;
    let x3 = v13[0] * ex[0] + v13[1] * ex[1] + v13[2] * ex[2];
    let y3 = v13[0] * ey[0] + v13[1] * ey[1] + v13[2] * ey[2];

    let area = 0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)).abs();

    (ex, ey, ez, [x1, y1, x2, y2, x3, y3], area)
}

// ─────────────────────────────────────────────────────────────
//  CST Membrane stiffness (6×6)
// ─────────────────────────────────────────────────────────────

/// Constant Strain Triangle membrane stiffness in local coordinates.
/// DOFs: [u1, v1, u2, v2, u3, v3].
fn compute_cst_membrane_ke(
    lc: &[f64; 6], area: f64, e_mod: f64, nu: f64, t: f64,
) -> [[f64; 6]; 6] {
    let (x1, y1, x2, y2, x3, y3) = (lc[0], lc[1], lc[2], lc[3], lc[4], lc[5]);

    let d_mem = [
        [e_mod / (1.0 - nu * nu), nu * e_mod / (1.0 - nu * nu), 0.0],
        [nu * e_mod / (1.0 - nu * nu), e_mod / (1.0 - nu * nu), 0.0],
        [0.0, 0.0, 0.5 * e_mod / (1.0 + nu)],
    ];

    let b_mem = [
        [y2 - y3, 0.0, y3 - y1, 0.0, y1 - y2, 0.0],
        [0.0, x3 - x2, 0.0, x1 - x3, 0.0, x2 - x1],
        [x3 - x2, y2 - y3, x1 - x3, y3 - y1, x2 - x1, y1 - y2],
    ];
    let inv_2a = 1.0 / (2.0 * area);

    // K_mem = t * A * (1/2A)^2 * B^T D B = t / (4A) * B^T D B
    let factor = t * inv_2a * inv_2a * area;
    let mut ke = [[0.0f64; 6]; 6];
    for i in 0..6 {
        for j in 0..6 {
            let mut s = 0.0;
            for k in 0..3 {
                for l in 0..3 {
                    s += b_mem[k][i] * d_mem[k][l] * b_mem[l][j];
                }
            }
            ke[i][j] = s * factor;
        }
    }
    ke
}

// ─────────────────────────────────────────────────────────────
//  DKT Bending B-matrix and stiffness (9×9)
// ─────────────────────────────────────────────────────────────

/// Evaluate the DKT curvature-displacement B-matrix (3×9) at parametric point (xi, eta).
///
/// Returns `bb[row][col]` where rows are `[kappa_x, kappa_y, kappa_xy]` and
/// columns correspond to DOFs `[w1, thx1, thy1, w2, thx2, thy2, w3, thx3, thy3]`.
///
/// `lc` = local coordinates `[x1,y1, x2,y2, x3,y3]` from `tri_local_frame`.
/// `two_a` = signed twice-area = `(x3-x1)*(y1-y2) - (x1-x2)*(y3-y1)`.
///
/// Reference: Batoz, Bathe, Ho (1980).
pub fn dkt_bending_b_matrix(
    lc: &[f64; 6], two_a: f64, xi: f64, eta: f64,
) -> [[f64; 9]; 3] {
    let (x1, y1, x2, y2, x3, y3) = (lc[0], lc[1], lc[2], lc[3], lc[4], lc[5]);

    let x23 = x2 - x3;
    let x31 = x3 - x1;
    let x12 = x1 - x2;
    let y23 = y2 - y3;
    let y31 = y3 - y1;
    let y12 = y1 - y2;

    let l23_sq = x23 * x23 + y23 * y23;
    let l31_sq = x31 * x31 + y31 * y31;
    let l12_sq = x12 * x12 + y12 * y12;

    let p4 = -6.0 * x23 / l23_sq;
    let p5 = -6.0 * x31 / l31_sq;
    let p6 = -6.0 * x12 / l12_sq;
    let q4 = 3.0 * x23 * y23 / l23_sq;
    let q5 = 3.0 * x31 * y31 / l31_sq;
    let q6 = 3.0 * x12 * y12 / l12_sq;
    let r4 = 3.0 * y23 * y23 / l23_sq;
    let r5 = 3.0 * y31 * y31 / l31_sq;
    let r6 = 3.0 * y12 * y12 / l12_sq;
    let t4 = -6.0 * y23 / l23_sq;
    let t5 = -6.0 * y31 / l31_sq;
    let t6 = -6.0 * y12 / l12_sq;

    let dhx_dxi = [
        p6 * (1.0 - 2.0 * xi) + (p5 - p6) * eta,
        q6 * (1.0 - 2.0 * xi) - (q5 + q6) * eta,
        -4.0 + 6.0 * (xi + eta) + r6 * (1.0 - 2.0 * xi) - eta * (r5 + r6),
        -p6 * (1.0 - 2.0 * xi) + eta * (p4 + p6),
        q6 * (1.0 - 2.0 * xi) - eta * (q6 - q4),
        -2.0 + 6.0 * xi + r6 * (1.0 - 2.0 * xi) + eta * (r4 - r6),
        -eta * (p5 + p4),
        eta * (q4 - q5),
        -eta * (r5 - r4),
    ];

    let dhy_dxi = [
        t6 * (1.0 - 2.0 * xi) + eta * (t5 - t6),
        1.0 + r6 * (1.0 - 2.0 * xi) - eta * (r5 + r6),
        -q6 * (1.0 - 2.0 * xi) + eta * (q5 + q6),
        -t6 * (1.0 - 2.0 * xi) + eta * (t4 + t6),
        -1.0 + r6 * (1.0 - 2.0 * xi) + eta * (r4 - r6),
        -q6 * (1.0 - 2.0 * xi) - eta * (q4 - q6),
        -eta * (t4 + t5),
        eta * (r4 - r5),
        -eta * (q4 - q5),
    ];

    let dhx_deta = [
        -p5 * (1.0 - 2.0 * eta) - xi * (p6 - p5),
        q5 * (1.0 - 2.0 * eta) - xi * (q5 + q6),
        -4.0 + 6.0 * (xi + eta) + r5 * (1.0 - 2.0 * eta) - xi * (r5 + r6),
        xi * (p4 + p6),
        xi * (q4 - q6),
        -xi * (r6 - r4),
        p5 * (1.0 - 2.0 * eta) - xi * (p4 + p5),
        q5 * (1.0 - 2.0 * eta) + xi * (q4 - q5),
        -2.0 + 6.0 * eta + r5 * (1.0 - 2.0 * eta) + xi * (r4 - r5),
    ];

    let dhy_deta = [
        -t5 * (1.0 - 2.0 * eta) - xi * (t6 - t5),
        1.0 + r5 * (1.0 - 2.0 * eta) - xi * (r5 + r6),
        -q5 * (1.0 - 2.0 * eta) + xi * (q5 + q6),
        xi * (t4 + t6),
        xi * (r4 - r6),
        -xi * (q4 - q6),
        t5 * (1.0 - 2.0 * eta) - xi * (t4 + t5),
        -1.0 + r5 * (1.0 - 2.0 * eta) + xi * (r4 - r5),
        -q5 * (1.0 - 2.0 * eta) - xi * (q4 - q5),
    ];

    let inv_2a = 1.0 / two_a;
    let mut bb = [[0.0f64; 9]; 3];
    for j in 0..9 {
        bb[0][j] = (y31 * dhx_dxi[j] + y12 * dhx_deta[j]) * inv_2a;
        bb[1][j] = (-x31 * dhy_dxi[j] - x12 * dhy_deta[j]) * inv_2a;
        bb[2][j] = (-x31 * dhx_dxi[j] - x12 * dhx_deta[j]
            + y31 * dhy_dxi[j] + y12 * dhy_deta[j])
            * inv_2a;
    }
    bb
}

/// DKT (Discrete Kirchhoff Triangle) bending stiffness in local coordinates.
/// DOFs: [w1, thx1, thy1, w2, thx2, thy2, w3, thx3, thy3].
/// Reference: Batoz, Bathe, Ho (1980).
fn compute_dkt_bending_ke(
    lc: &[f64; 6], _area: f64, e_mod: f64, nu: f64, t: f64,
) -> [[f64; 9]; 9] {
    let (x1, y1, x2, y2, x3, y3) = (lc[0], lc[1], lc[2], lc[3], lc[4], lc[5]);
    let x31 = x3 - x1;
    let x12 = x1 - x2;
    let y31 = y3 - y1;
    let y12 = y1 - y2;
    let two_a = x31 * y12 - x12 * y31;

    let d_coeff = e_mod * t * t * t / (12.0 * (1.0 - nu * nu));
    let db = [
        [d_coeff, nu * d_coeff, 0.0],
        [nu * d_coeff, d_coeff, 0.0],
        [0.0, 0.0, d_coeff * (1.0 - nu) / 2.0],
    ];

    let mut ke = [[0.0f64; 9]; 9];

    for &([xi, eta], wt) in &TRI_GAUSS_PTS {
        let bb = dkt_bending_b_matrix(lc, two_a, xi, eta);

        let factor = wt * two_a.abs();
        for i in 0..9 {
            for j in i..9 {
                let mut s = 0.0;
                for k in 0..3 {
                    for l in 0..3 {
                        s += bb[k][i] * db[k][l] * bb[l][j];
                    }
                }
                let val = s * factor;
                ke[i][j] += val;
                if i != j {
                    ke[j][i] += val;
                }
            }
        }
    }

    ke
}

// ─────────────────────────────────────────────────────────────
//  Hughes-Brezzi drilling stiffness (9×9)
// ─────────────────────────────────────────────────────────────

/// Hughes-Brezzi drilling DOF stiffness in local coordinates.
/// DOFs: [u1, v1, thz1, u2, v2, thz2, u3, v3, thz3].
/// gamma = G = E / (2(1+nu)) (shear modulus).
/// Reference: Hughes & Brezzi (1989).
fn compute_drilling_ke(
    lc: &[f64; 6], area: f64, e_mod: f64, nu: f64, t: f64,
) -> [[f64; 9]; 9] {
    let (_x1, _y1, x2, y2, x3, y3) = (lc[0], lc[1], lc[2], lc[3], lc[4], lc[5]);

    // CST shape function gradients (bi = yj - yk, ci = xk - xj, cyclic 1→2→3)
    let b1 = y2 - y3;
    let b2 = y3; // y3 - y1, but y1 = 0
    let b3 = -y2; // y1 - y2, but y1 = 0
    let c1 = x3 - x2;
    let c2 = -x3; // x1 - x3, but x1 = 0
    let c3 = x2; // x2 - x1, but x1 = 0

    let gamma = e_mod / (2.0 * (1.0 + nu));
    let inv_4a = 1.0 / (4.0 * area);

    // B3 is a 1×9 row vector (function of position via shape functions Ni).
    // B3 = [c1/4A, -b1/4A, N1, c2/4A, -b2/4A, N2, c3/4A, -b3/4A, N3]
    // The constant part (from displacement gradients) is position-independent for CST.
    // The Ni part varies linearly. We need integral(B3^T * B3 * dA).
    //
    // Split B3 = B_const + B_shape where:
    //   B_const[3*i]   = ci / (4A)      (constant)
    //   B_const[3*i+1] = -bi / (4A)     (constant)
    //   B_const[3*i+2] = 0              (shape function part)
    //   B_shape[3*i+2] = Ni             (linear)
    //
    // integral(B3^T B3 dA) has three types of terms:
    //   const*const: integral = A * val
    //   const*Ni:    integral = A/3 * val  (since integral(Ni dA) = A/3)
    //   Ni*Nj:       integral = A/12*(1+delta_ij)

    let bc = [
        c1 * inv_4a, -b1 * inv_4a,
        c2 * inv_4a, -b2 * inv_4a,
        c3 * inv_4a, -b3 * inv_4a,
    ];

    let mut ke = [[0.0f64; 9]; 9];
    let scale = gamma * t;

    for i in 0..3 {
        for j in 0..3 {
            let ci_val = bc[2 * i];
            let bi_val = bc[2 * i + 1];
            let cj_val = bc[2 * j];
            let bj_val = bc[2 * j + 1];

            // (u_i, u_j): const*const → A * ci*cj
            ke[3 * i][3 * j] += scale * area * ci_val * cj_val;
            // (u_i, v_j): const*const → A * ci*bj
            ke[3 * i][3 * j + 1] += scale * area * ci_val * bj_val;
            // (v_i, u_j): const*const → A * bi*cj
            ke[3 * i + 1][3 * j] += scale * area * bi_val * cj_val;
            // (v_i, v_j): const*const → A * bi*bj
            ke[3 * i + 1][3 * j + 1] += scale * area * bi_val * bj_val;

            // (u_i, thz_j): const*Nj → A/3 * ci
            ke[3 * i][3 * j + 2] += scale * (area / 3.0) * ci_val;
            // (thz_j, u_i): symmetric
            ke[3 * j + 2][3 * i] += scale * (area / 3.0) * ci_val;

            // (v_i, thz_j): const*Nj → A/3 * bi
            ke[3 * i + 1][3 * j + 2] += scale * (area / 3.0) * bi_val;
            // (thz_j, v_i): symmetric
            ke[3 * j + 2][3 * i + 1] += scale * (area / 3.0) * bi_val;

            // (thz_i, thz_j): Ni*Nj → A/12*(1+delta_ij)
            let mass_ij = if i == j {
                area / 6.0
            } else {
                area / 12.0
            };
            ke[3 * i + 2][3 * j + 2] += scale * mass_ij;
        }
    }

    ke
}

// ─────────────────────────────────────────────────────────────
//  Full element stiffness (18×18) and global assembly
// ─────────────────────────────────────────────────────────────

/// Compute the full 18×18 global element stiffness for a flat shell triangle.
/// Combines CST membrane + DKT bending + Hughes-Brezzi drilling,
/// then transforms from local to global coordinates.
pub fn compute_tri_shell_ke_global(e: usize, problem: &ShellProblem, ke: &mut [f64]) {
    let nodes = &problem.elements[e];
    let p1 = problem.node_positions.row(nodes[0]);
    let p2 = problem.node_positions.row(nodes[1]);
    let p3 = problem.node_positions.row(nodes[2]);

    let props = &problem.element_props[e];
    let mat = &problem.materials[props.material_idx];

    let t = (problem.node_thicknesses[nodes[0]]
        + problem.node_thicknesses[nodes[1]]
        + problem.node_thicknesses[nodes[2]])
        / 3.0;

    let (ex, ey, ez, lc, area) =
        tri_local_frame(p1.as_slice().unwrap(), p2.as_slice().unwrap(), p3.as_slice().unwrap());

    let e_mod = mat.e;
    let nu = mat.nu;

    // Sub-element stiffness matrices in local coordinates
    let ke_mem = compute_cst_membrane_ke(&lc, area, e_mod, nu, t);
    let ke_bend = compute_dkt_bending_ke(&lc, area, e_mod, nu, t);
    let ke_drill = compute_drilling_ke(&lc, area, e_mod, nu, t);

    // Assemble into 18×18 local Ke
    // Local DOF order per node: [u, v, w, thx, thy, thz]
    let mut ke_local = [[0.0f64; 18]; 18];

    // Membrane: DOFs [u1,v1, u2,v2, u3,v3] → indices [0,1, 6,7, 12,13]
    let m_map: [usize; 6] = [0, 1, 6, 7, 12, 13];
    for i in 0..6 {
        for j in 0..6 {
            ke_local[m_map[i]][m_map[j]] += ke_mem[i][j];
        }
    }

    // Bending: DOFs [w1,thx1,thy1, w2,thx2,thy2, w3,thx3,thy3] → indices [2,3,4, 8,9,10, 14,15,16]
    let b_map: [usize; 9] = [2, 3, 4, 8, 9, 10, 14, 15, 16];
    for i in 0..9 {
        for j in 0..9 {
            ke_local[b_map[i]][b_map[j]] += ke_bend[i][j];
        }
    }

    // Drilling: DOFs [u1,v1,thz1, u2,v2,thz2, u3,v3,thz3] → indices [0,1,5, 6,7,11, 12,13,17]
    let d_map: [usize; 9] = [0, 1, 5, 6, 7, 11, 12, 13, 17];
    for i in 0..9 {
        for j in 0..9 {
            ke_local[d_map[i]][d_map[j]] += ke_drill[i][j];
        }
    }

    // Build 18×18 rotation matrix T (block-diagonal: 6 blocks of 3×3 rotation)
    let rot = [
        [ex[0], ex[1], ex[2]],
        [ey[0], ey[1], ey[2]],
        [ez[0], ez[1], ez[2]],
    ];

    let mut t_mat = [[0.0f64; 18]; 18];
    for n in 0..3 {
        for i in 0..3 {
            for j in 0..3 {
                t_mat[n * 6 + i][n * 6 + j] = rot[i][j];
                t_mat[n * 6 + 3 + i][n * 6 + 3 + j] = rot[i][j];
            }
        }
    }

    // Global Ke = T^T * Ke_local * T
    // Two-step: temp = Ke_local * T, then Ke_global = T^T * temp
    let mut temp = [[0.0f64; 18]; 18];
    for i in 0..18 {
        for j in 0..18 {
            let mut s = 0.0;
            for k in 0..18 {
                s += ke_local[i][k] * t_mat[k][j];
            }
            temp[i][j] = s;
        }
    }

    for i in 0..18 {
        for j in 0..18 {
            let mut s = 0.0;
            for k in 0..18 {
                s += t_mat[k][i] * temp[k][j];
            }
            ke[i * 18 + j] = s;
        }
    }
}

/// Assemble the global stiffness matrix K for shell elements.
pub fn assemble_shell_k(
    cache: &mut ShellCache,
    problem: &ShellProblem,
) -> Result<(), TheseusError> {
    let nnz = cache.k_matrix.values.len();
    let ne = problem.num_elements;
    let elem_to_nz = &cache.element_to_nz;

    let merged = (0..ne)
        .into_par_iter()
        .fold(
            || vec![0.0f64; nnz],
            |mut buf, e| {
                let nodes = &problem.elements[e];
                let num_nodes = nodes.len();
                let num_dofs = num_nodes * 6;

                let mut ke_flat = vec![0.0f64; num_dofs * num_dofs];

                if num_nodes == 3 {
                    compute_tri_shell_ke_global(e, problem, &mut ke_flat);
                }

                for &(nz_idx, li, lj) in &elem_to_nz.entries[e] {
                    buf[nz_idx] += ke_flat[li * num_dofs + lj];
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
//  Self-weight: consistent nodal forces via Gauss quadrature
// ─────────────────────────────────────────────────────────────

/// Consistent self-weight nodal forces for a triangle shell element.
/// Returns `[f_node0, f_node1, f_node2]` where each is a 3D force vector.
///
/// Uses 3-point Hammer quadrature with linear shape functions and linearly
/// interpolated thickness, giving **exact** integration for the quadratic
/// integrand `N_i(xi,eta) * t(xi,eta)`.
pub fn shell_element_selfweight(
    problem: &ShellProblem, e: usize,
) -> [[f64; 3]; 3] {
    let nodes = &problem.elements[e];
    let mat = &problem.materials[problem.element_props[e].material_idx];
    let rho = mat.density;

    let p1 = problem.node_positions.row(nodes[0]);
    let p2 = problem.node_positions.row(nodes[1]);
    let p3 = problem.node_positions.row(nodes[2]);
    let (_, _, _, _, area) = tri_local_frame(
        p1.as_slice().unwrap(), p2.as_slice().unwrap(), p3.as_slice().unwrap());

    let two_area = 2.0 * area;
    let t_nodes = [
        problem.node_thicknesses[nodes[0]],
        problem.node_thicknesses[nodes[1]],
        problem.node_thicknesses[nodes[2]],
    ];

    let mut forces = [[0.0f64; 3]; 3];

    for &([xi, eta], wt) in &TRI_GAUSS_PTS {
        let n_vals = [1.0 - xi - eta, xi, eta];
        let t_gp = n_vals[0] * t_nodes[0] + n_vals[1] * t_nodes[1] + n_vals[2] * t_nodes[2];
        let factor = wt * two_area * rho * t_gp;

        for node_local in 0..3 {
            let f = factor * n_vals[node_local];
            for d in 0..3 {
                forces[node_local][d] += f * problem.gravity[d];
            }
        }
    }
    forces
}

// ─────────────────────────────────────────────────────────────
//  Load assembly
// ─────────────────────────────────────────────────────────────

/// Assemble the RHS force vector (free DOFs only) including forces, moments, and self-weight.
pub fn assemble_shell_loads(cache: &mut ShellCache, problem: &ShellProblem) {
    let n_free = cache.dof_map.num_free_dofs;
    for i in 0..n_free {
        cache.rhs[i] = 0.0;
    }

    for load in &problem.loads {
        for d in 0..3 {
            let global_dof = load.node_idx * 6 + d;
            if let Some(free_idx) = cache.dof_map.global_to_free[global_dof] {
                cache.rhs[free_idx] += load.force[d];
            }
        }
        for d in 0..3 {
            let global_dof = load.node_idx * 6 + 3 + d;
            if let Some(free_idx) = cache.dof_map.global_to_free[global_dof] {
                cache.rhs[free_idx] += load.moment[d];
            }
        }
    }

    if problem.include_self_weight {
        for e in 0..problem.num_elements {
            if problem.elements[e].len() != 3 { continue; }
            let forces = shell_element_selfweight(problem, e);
            let nodes = &problem.elements[e];
            for node_local in 0..3 {
                let ni = nodes[node_local];
                for d in 0..3 {
                    let gi = ni * 6 + d;
                    if let Some(fi) = cache.dof_map.global_to_free[gi] {
                        cache.rhs[fi] += forces[node_local][d];
                    }
                }
            }
        }
    }
}
