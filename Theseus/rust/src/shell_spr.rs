//! Zienkiewicz-Zhu Superconvergent Patch Recovery (SPR) for shell elements
//! using a 2D tangent-plane polynomial basis.
//!
//! Shell element centroids are (nearly) coplanar, making the solid SPR's 3D
//! linear basis ill-conditioned. This module uses a 2D basis [1, s, t] in
//! each node's tangent plane, requiring only 4 adjacent elements for a
//! well-conditioned fit (vs 6 for the solid 3D basis).

use crate::shell_assembly::tri_local_frame;
use crate::solid_spr::{principal_stresses, von_mises_from_principals};
use ndarray::Array2;
use rayon::prelude::*;

// ─────────────────────────────────────────────────────────────
//  2D tangent-plane SPR recovery
// ─────────────────────────────────────────────────────────────

/// Recover smooth nodal stresses from per-element shell stresses using
/// Zienkiewicz-Zhu SPR with a 2D tangent-plane polynomial basis.
///
/// For each node, collects adjacent element stresses and centroids, builds a
/// tangent plane from area-weighted element normals, projects centroids to 2D,
/// and fits a centered linear polynomial sigma*(s,t) = a0 + a1*s + a2*t.
///
/// Falls back to area-weighted averaging when the patch has fewer than 4
/// elements, the 3x3 normal matrix is ill-conditioned, or the fit produces
/// outlier values.
///
/// `element_stresses` must be 6-component Voigt tensors in the global frame.
///
/// Returns `(nodal_stresses, element_errors)`.
pub fn spr_shell_recover_nodal_stresses(
    node_positions: &Array2<f64>,
    elements: &[Vec<usize>],
    element_stresses: &[[f64; 6]],
) -> (Vec<[f64; 6]>, Vec<f64>) {
    let nn = node_positions.nrows();
    let ne = elements.len();

    // 1. Build node-to-element adjacency
    let mut node_to_elems: Vec<Vec<usize>> = vec![Vec::new(); nn];
    for (e, elem) in elements.iter().enumerate() {
        for &ni in elem {
            node_to_elems[ni].push(e);
        }
    }

    // 2. Compute element centroids and areas
    let mut centroids: Vec<[f64; 3]> = Vec::with_capacity(ne);
    let mut areas: Vec<f64> = Vec::with_capacity(ne);
    let mut normals: Vec<[f64; 3]> = Vec::with_capacity(ne);

    for elem in elements {
        if elem.len() == 3 {
            let p1 = [
                node_positions[[elem[0], 0]],
                node_positions[[elem[0], 1]],
                node_positions[[elem[0], 2]],
            ];
            let p2 = [
                node_positions[[elem[1], 0]],
                node_positions[[elem[1], 1]],
                node_positions[[elem[1], 2]],
            ];
            let p3 = [
                node_positions[[elem[2], 0]],
                node_positions[[elem[2], 1]],
                node_positions[[elem[2], 2]],
            ];
            let (_ex, _ey, ez, _lc, area) = tri_local_frame(&p1, &p2, &p3);
            centroids.push([
                (p1[0] + p2[0] + p3[0]) / 3.0,
                (p1[1] + p2[1] + p3[1]) / 3.0,
                (p1[2] + p2[2] + p3[2]) / 3.0,
            ]);
            areas.push(area);
            normals.push(ez);
        } else {
            centroids.push([0.0; 3]);
            areas.push(0.0);
            normals.push([0.0, 0.0, 1.0]);
        }
    }

    // 3. Recover nodal stresses in parallel
    let nodal_stresses: Vec<[f64; 6]> = (0..nn)
        .into_par_iter()
        .map(|ni| {
            let adj = &node_to_elems[ni];
            let n_adj = adj.len();

            if n_adj == 0 {
                return [0.0; 6];
            }

            let node_pos = [
                node_positions[[ni, 0]],
                node_positions[[ni, 1]],
                node_positions[[ni, 2]],
            ];

            if n_adj < 4 {
                return area_weighted_average(&areas, adj, element_stresses);
            }

            // Build tangent plane from area-weighted average normal
            let (t1, t2) = build_tangent_basis(&normals, &areas, adj);

            // Project centroids to tangent-plane coordinates (centered at node)
            let mut ptp = [[0.0f64; 3]; 3];
            let mut pts = [[0.0f64; 3]; 6];
            let mut max_stress_sq = 0.0f64;

            for &e in adj {
                let ds = [
                    centroids[e][0] - node_pos[0],
                    centroids[e][1] - node_pos[1],
                    centroids[e][2] - node_pos[2],
                ];
                let s = ds[0] * t1[0] + ds[1] * t1[1] + ds[2] * t1[2];
                let t = ds[0] * t2[0] + ds[1] * t2[1] + ds[2] * t2[2];

                let p = [1.0, s, t];
                for i in 0..3 {
                    for j in 0..3 {
                        ptp[i][j] += p[i] * p[j];
                    }
                    for c in 0..6 {
                        pts[c][i] += p[i] * element_stresses[e][c];
                    }
                }
                let s_sq: f64 = element_stresses[e].iter().map(|v| v * v).sum();
                if s_sq > max_stress_sq {
                    max_stress_sq = s_sq;
                }
            }

            let max_stress = max_stress_sq.sqrt();

            let coeffs = match solve_3x3_multi(&ptp, &pts) {
                Some(c) => c,
                None => return area_weighted_average(&areas, adj, element_stresses),
            };

            // With centered basis, P_node = [1, 0, 0] → result is a0
            let mut result = [0.0; 6];
            for c in 0..6 {
                result[c] = coeffs[c][0];
            }

            let threshold = 3.0 * max_stress;
            if threshold > 1e-30 {
                for c in 0..6 {
                    if result[c].abs() > threshold {
                        return area_weighted_average(&areas, adj, element_stresses);
                    }
                }
            }

            result
        })
        .collect();

    // 4. Compute SPR error estimator per element
    let element_errors: Vec<f64> = (0..ne)
        .into_par_iter()
        .map(|e| {
            let elem = &elements[e];
            let num_nodes = elem.len();
            let mut spr_avg = [0.0; 6];
            for &ni in elem {
                for c in 0..6 {
                    spr_avg[c] += nodal_stresses[ni][c];
                }
            }
            let inv_n = 1.0 / (num_nodes as f64);
            for c in 0..6 {
                spr_avg[c] *= inv_n;
            }

            let mut diff_norm_sq = 0.0;
            let mut spr_norm_sq = 0.0;
            for c in 0..6 {
                let diff = spr_avg[c] - element_stresses[e][c];
                diff_norm_sq += diff * diff;
                spr_norm_sq += spr_avg[c] * spr_avg[c];
            }

            let spr_norm = spr_norm_sq.sqrt();
            if spr_norm < 1e-30 {
                0.0
            } else {
                diff_norm_sq.sqrt() / spr_norm
            }
        })
        .collect();

    (nodal_stresses, element_errors)
}

// ─────────────────────────────────────────────────────────────
//  Tangent-plane construction
// ─────────────────────────────────────────────────────────────

/// Build an orthonormal tangent basis (t1, t2) at a node from area-weighted
/// average of adjacent element normals.
fn build_tangent_basis(
    normals: &[[f64; 3]],
    areas: &[f64],
    adj: &[usize],
) -> ([f64; 3], [f64; 3]) {
    let mut avg_n = [0.0f64; 3];
    for &e in adj {
        let w = areas[e];
        for d in 0..3 {
            avg_n[d] += w * normals[e][d];
        }
    }
    let len = (avg_n[0] * avg_n[0] + avg_n[1] * avg_n[1] + avg_n[2] * avg_n[2]).sqrt();
    if len < 1e-30 {
        return ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    }
    let n = [avg_n[0] / len, avg_n[1] / len, avg_n[2] / len];

    // Choose a seed vector not parallel to n
    let seed = if n[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };

    // t1 = normalize(seed - (seed . n) * n)
    let dot = seed[0] * n[0] + seed[1] * n[1] + seed[2] * n[2];
    let mut t1 = [seed[0] - dot * n[0], seed[1] - dot * n[1], seed[2] - dot * n[2]];
    let t1_len = (t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]).sqrt();
    if t1_len < 1e-30 {
        return ([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    }
    t1[0] /= t1_len;
    t1[1] /= t1_len;
    t1[2] /= t1_len;

    // t2 = n × t1
    let t2 = [
        n[1] * t1[2] - n[2] * t1[1],
        n[2] * t1[0] - n[0] * t1[2],
        n[0] * t1[1] - n[1] * t1[0],
    ];

    (t1, t2)
}

// ─────────────────────────────────────────────────────────────
//  Area-weighted average fallback
// ─────────────────────────────────────────────────────────────

fn area_weighted_average(
    areas: &[f64],
    adj: &[usize],
    element_stresses: &[[f64; 6]],
) -> [f64; 6] {
    let mut result = [0.0; 6];
    let mut total_a = 0.0;
    for &e in adj {
        let a = areas[e].max(1e-30);
        total_a += a;
        for c in 0..6 {
            result[c] += a * element_stresses[e][c];
        }
    }
    let inv_a = 1.0 / total_a;
    for c in 0..6 {
        result[c] *= inv_a;
    }
    result
}

// ─────────────────────────────────────────────────────────────
//  3×3 linear system solver for 6 RHS (Gaussian elimination)
// ─────────────────────────────────────────────────────────────

/// Solve 3×3 linear system A x = b for 6 right-hand sides simultaneously.
/// Returns `None` if the system is near-singular.
fn solve_3x3_multi(a: &[[f64; 3]; 3], rhs: &[[f64; 3]; 6]) -> Option<[[f64; 3]; 6]> {
    let mut aug_a = *a;
    let mut aug_rhs = *rhs;
    let mut max_pivot = 0.0f64;

    for col in 0..3 {
        let mut max_val = aug_a[col][col].abs();
        let mut max_row = col;
        for row in (col + 1)..3 {
            let v = aug_a[row][col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }

        if max_val > max_pivot {
            max_pivot = max_val;
        }

        let tol = if max_pivot > 1e-30 { 1e-10 * max_pivot } else { 1e-30 };
        if max_val < tol {
            return None;
        }

        if max_row != col {
            aug_a.swap(col, max_row);
            for c in 0..6 {
                let tmp = aug_rhs[c][col];
                aug_rhs[c][col] = aug_rhs[c][max_row];
                aug_rhs[c][max_row] = tmp;
            }
        }

        let pivot = aug_a[col][col];
        for row in (col + 1)..3 {
            let factor = aug_a[row][col] / pivot;
            for j in col..3 {
                aug_a[row][j] -= factor * aug_a[col][j];
            }
            for c in 0..6 {
                aug_rhs[c][row] -= factor * aug_rhs[c][col];
            }
        }
    }

    let mut result = [[0.0f64; 3]; 6];
    for c in 0..6 {
        for i in (0..3).rev() {
            let mut s = aug_rhs[c][i];
            for j in (i + 1)..3 {
                s -= aug_a[i][j] * result[c][j];
            }
            result[c][i] = s / aug_a[i][i];
        }
    }
    Some(result)
}

// ─────────────────────────────────────────────────────────────
//  Full shell SPR orchestration: recover + principals + von Mises
// ─────────────────────────────────────────────────────────────

/// Result of the full shell SPR pipeline: smoothed nodal stresses, principal
/// decomposition, and von Mises for top and bottom fibers.
pub struct ShellSprResult {
    pub nodal_membrane: Vec<[f64; 6]>,
    pub nodal_top: Vec<[f64; 6]>,
    pub nodal_bottom: Vec<[f64; 6]>,
    pub element_errors: Vec<f64>,
    pub principal_values_top: Vec<[f64; 3]>,
    pub principal_values_bot: Vec<[f64; 3]>,
    pub principal_vectors_top: Vec<[[f64; 3]; 3]>,
    pub principal_vectors_bot: Vec<[[f64; 3]; 3]>,
    pub von_mises_top: Vec<f64>,
    pub von_mises_bot: Vec<f64>,
}

/// Full shell SPR pipeline: recover nodal stresses for membrane, top, and
/// bottom fiber fields, then decompose into principals and von Mises.
pub fn spr_shell_full(
    node_positions: &Array2<f64>,
    elements: &[Vec<usize>],
    membrane_stresses: &[[f64; 6]],
    top_stresses: &[[f64; 6]],
    bottom_stresses: &[[f64; 6]],
) -> ShellSprResult {
    let (nodal_membrane, element_errors) =
        spr_shell_recover_nodal_stresses(node_positions, elements, membrane_stresses);
    let (nodal_top, _) =
        spr_shell_recover_nodal_stresses(node_positions, elements, top_stresses);
    let (nodal_bottom, _) =
        spr_shell_recover_nodal_stresses(node_positions, elements, bottom_stresses);

    let nn = node_positions.nrows();

    // Principal decomposition + von Mises in parallel
    let top_decomp: Vec<([f64; 3], [[f64; 3]; 3], f64)> = nodal_top
        .par_iter()
        .map(|s| {
            let (vals, vecs) = principal_stresses(s);
            let vm = von_mises_from_principals(vals[0], vals[1], vals[2]);
            (vals, vecs, vm)
        })
        .collect();

    let bot_decomp: Vec<([f64; 3], [[f64; 3]; 3], f64)> = nodal_bottom
        .par_iter()
        .map(|s| {
            let (vals, vecs) = principal_stresses(s);
            let vm = von_mises_from_principals(vals[0], vals[1], vals[2]);
            (vals, vecs, vm)
        })
        .collect();

    let mut principal_values_top = Vec::with_capacity(nn);
    let mut principal_vectors_top = Vec::with_capacity(nn);
    let mut von_mises_top = Vec::with_capacity(nn);
    for (vals, vecs, vm) in &top_decomp {
        principal_values_top.push(*vals);
        principal_vectors_top.push(*vecs);
        von_mises_top.push(*vm);
    }

    let mut principal_values_bot = Vec::with_capacity(nn);
    let mut principal_vectors_bot = Vec::with_capacity(nn);
    let mut von_mises_bot = Vec::with_capacity(nn);
    for (vals, vecs, vm) in &bot_decomp {
        principal_values_bot.push(*vals);
        principal_vectors_bot.push(*vecs);
        von_mises_bot.push(*vm);
    }

    ShellSprResult {
        nodal_membrane,
        nodal_top,
        nodal_bottom,
        element_errors,
        principal_values_top,
        principal_values_bot,
        principal_vectors_top,
        principal_vectors_bot,
        von_mises_top,
        von_mises_bot,
    }
}

// ─────────────────────────────────────────────────────────────
//  Tests
// ─────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_spr_single_triangle() {
        let positions = Array2::from_shape_vec(
            (3, 3),
            vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0],
        )
        .unwrap();
        let elements = vec![vec![0, 1, 2]];
        let stresses = vec![[100.0, 50.0, 25.0, 10.0, 0.0, 0.0]];

        let (nodal, errors) =
            spr_shell_recover_nodal_stresses(&positions, &elements, &stresses);

        assert_eq!(nodal.len(), 3);
        assert_eq!(errors.len(), 1);

        for ni in 0..3 {
            for c in 0..6 {
                assert!(
                    (nodal[ni][c] - stresses[0][c]).abs() < 1e-6,
                    "node {ni}, component {c}: expected {}, got {}",
                    stresses[0][c],
                    nodal[ni][c]
                );
            }
        }
    }

    #[test]
    fn test_spr_fan_patch_no_blowup() {
        // 6 triangles in a fan sharing node 0, lying in the z=0 plane
        let positions = Array2::from_shape_vec(
            (7, 3),
            vec![
                0.0, 0.0, 0.0,
                2.0, 0.0, 0.0,
                1.0, 1.73, 0.0,
                -1.0, 1.73, 0.0,
                -2.0, 0.0, 0.0,
                -1.0, -1.73, 0.0,
                1.0, -1.73, 0.0,
            ],
        )
        .unwrap();

        let elements = vec![
            vec![0, 1, 2],
            vec![0, 2, 3],
            vec![0, 3, 4],
            vec![0, 4, 5],
            vec![0, 5, 6],
            vec![0, 6, 1],
        ];

        let stresses = vec![
            [100.0, 50.0, 25.0, 10.0, 0.0, 0.0],
            [110.0, 55.0, 28.0, 12.0, 0.0, 0.0],
            [90.0, 45.0, 22.0, 8.0, 0.0, 0.0],
            [105.0, 52.0, 26.0, 11.0, 0.0, 0.0],
            [95.0, 48.0, 24.0, 9.0, 0.0, 0.0],
            [102.0, 51.0, 25.5, 10.5, 0.0, 0.0],
        ];

        let (nodal, _) =
            spr_shell_recover_nodal_stresses(&positions, &elements, &stresses);

        for c in 0..6 {
            let elem_max = stresses.iter().map(|s| s[c].abs()).fold(0.0f64, f64::max);
            assert!(
                nodal[0][c].abs() < 5.0 * elem_max.max(1.0),
                "node 0, component {c}: recovered {:.1} unreasonably large (max elem: {:.1})",
                nodal[0][c],
                elem_max,
            );
        }
    }

    #[test]
    fn test_spr_boundary_node_fallback() {
        // Node 0 is on the boundary with only 2 adjacent elements → fallback
        let positions = Array2::from_shape_vec(
            (4, 3),
            vec![
                0.0, 0.0, 0.0,
                1.0, 0.0, 0.0,
                0.5, 1.0, 0.0,
                -0.5, 1.0, 0.0,
            ],
        )
        .unwrap();

        let elements = vec![vec![0, 1, 2], vec![0, 2, 3]];
        let stresses = vec![
            [100.0, 50.0, 25.0, 10.0, 0.0, 0.0],
            [80.0, 40.0, 20.0, 8.0, 0.0, 0.0],
        ];

        let (nodal, _) =
            spr_shell_recover_nodal_stresses(&positions, &elements, &stresses);

        // Node 0 should get area-weighted average (2 elements < 4 threshold)
        for c in 0..6 {
            let elem_max = stresses.iter().map(|s| s[c].abs()).fold(0.0f64, f64::max);
            assert!(
                nodal[0][c].abs() <= elem_max * 1.5,
                "boundary node 0, component {c}: {:.1} exceeds reasonable bound",
                nodal[0][c],
            );
        }
    }
}
