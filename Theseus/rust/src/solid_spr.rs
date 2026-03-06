//! Zienkiewicz–Zhu Superconvergent Patch Recovery (SPR) for tet4 elements
//! and closed-form principal stress decomposition via Cardano's method.

use ndarray::Array2;
use rayon::prelude::*;

// ─────────────────────────────────────────────────────────────
//  SPR nodal stress recovery
// ─────────────────────────────────────────────────────────────

/// Recover smooth nodal stresses from constant-strain tet4 element stresses
/// using the Zienkiewicz–Zhu Superconvergent Patch Recovery.
///
/// For each node, collects all adjacent element stresses and their centroid
/// positions, fits a **centered** linear polynomial
/// σ*(x,y,z) = a₀ + a₁(x−xₙ) + a₂(y−yₙ) + a₃(z−zₙ)
/// via least squares. Centering ensures a₀ directly estimates the stress at
/// the node, eliminating extrapolation amplification of gradient errors.
///
/// Falls back to inverse-distance-weighted averaging when the patch has fewer
/// than 6 elements, the normal matrix is near-singular, or the fit produces
/// outlier values.
///
/// Returns `(nodal_stresses, element_errors)` where:
/// - `nodal_stresses[i]` is the 6-component recovered stress tensor at node i
/// - `element_errors[e]` is the SPR error estimator η_e for element e
pub fn spr_recover_nodal_stresses(
    node_positions: &Array2<f64>,
    elements: &[[usize; 4]],
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

    // 2. Compute element centroids
    let centroids: Vec<[f64; 3]> = elements
        .iter()
        .map(|elem| {
            let mut c = [0.0; 3];
            for &ni in elem {
                for d in 0..3 {
                    c[d] += node_positions[[ni, d]];
                }
            }
            for d in 0..3 {
                c[d] *= 0.25;
            }
            c
        })
        .collect();

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

            if n_adj < 6 {
                return idw_average(&centroids, adj, &node_pos, element_stresses);
            }

            // Centered polynomial basis: P = [1, x−xₙ, y−yₙ, z−zₙ]
            let mut ptp = [[0.0f64; 4]; 4];
            let mut pts = [[0.0f64; 4]; 6];
            let mut max_stress_sq = 0.0f64;

            for &e in adj {
                let p = [
                    1.0,
                    centroids[e][0] - node_pos[0],
                    centroids[e][1] - node_pos[1],
                    centroids[e][2] - node_pos[2],
                ];
                for i in 0..4 {
                    for j in 0..4 {
                        ptp[i][j] += p[i] * p[j];
                    }
                    for c in 0..6 {
                        pts[c][i] += p[i] * element_stresses[e][c];
                    }
                }
                let s_sq: f64 = element_stresses[e].iter().map(|v| v * v).sum();
                if s_sq > max_stress_sq { max_stress_sq = s_sq; }
            }

            let max_stress = max_stress_sq.sqrt();

            let coeffs = match solve_4x4_multi(&ptp, &pts) {
                Some(c) => c,
                None => {
                    return idw_average(&centroids, adj, &node_pos, element_stresses);
                }
            };

            // With centered basis, P_node = [1, 0, 0, 0] → result is simply a₀
            let mut result = [0.0; 6];
            for c in 0..6 {
                result[c] = coeffs[c][0];
            }

            // Sanity: if any component exceeds 3× the local stress magnitude,
            // the fit is unreliable — fall back to IDW.
            let threshold = 3.0 * max_stress;
            if threshold > 1e-30 {
                for c in 0..6 {
                    if result[c].abs() > threshold {
                        return idw_average(&centroids, adj, &node_pos, element_stresses);
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
            // Average the recovered nodal stresses over the element
            let mut spr_avg = [0.0; 6];
            for &ni in elem {
                for c in 0..6 {
                    spr_avg[c] += nodal_stresses[ni][c];
                }
            }
            for c in 0..6 {
                spr_avg[c] *= 0.25;
            }

            // ||sigma_SPR - sigma_h||
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

/// Inverse-distance-weighted average of element stresses around a node.
fn idw_average(
    centroids: &[[f64; 3]],
    adj: &[usize],
    node_pos: &[f64; 3],
    element_stresses: &[[f64; 6]],
) -> [f64; 6] {
    let mut result = [0.0; 6];
    let mut total_w = 0.0;
    for &e in adj {
        let dx = centroids[e][0] - node_pos[0];
        let dy = centroids[e][1] - node_pos[1];
        let dz = centroids[e][2] - node_pos[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-30);
        let w = 1.0 / dist;
        total_w += w;
        for c in 0..6 {
            result[c] += w * element_stresses[e][c];
        }
    }
    let inv_w = 1.0 / total_w;
    for c in 0..6 {
        result[c] *= inv_w;
    }
    result
}

/// Solve 4×4 linear system A x = b for 6 right-hand sides simultaneously
/// using Gaussian elimination with partial pivoting.
/// Returns `None` if the system is near-singular.
fn solve_4x4_multi(a: &[[f64; 4]; 4], rhs: &[[f64; 4]; 6]) -> Option<[[f64; 4]; 6]> {
    let mut aug_a = *a;
    let mut aug_rhs = *rhs;
    let mut max_pivot = 0.0f64;

    for col in 0..4 {
        // Partial pivoting
        let mut max_val = aug_a[col][col].abs();
        let mut max_row = col;
        for row in (col + 1)..4 {
            let v = aug_a[row][col].abs();
            if v > max_val {
                max_val = v;
                max_row = row;
            }
        }

        if max_val > max_pivot {
            max_pivot = max_val;
        }

        // Scale-relative near-singularity check
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

        // Eliminate below
        let pivot = aug_a[col][col];
        for row in (col + 1)..4 {
            let factor = aug_a[row][col] / pivot;
            for j in col..4 {
                aug_a[row][j] -= factor * aug_a[col][j];
            }
            for c in 0..6 {
                aug_rhs[c][row] -= factor * aug_rhs[c][col];
            }
        }
    }

    // Back-substitution
    let mut result = [[0.0f64; 4]; 6];
    for c in 0..6 {
        for i in (0..4).rev() {
            let mut s = aug_rhs[c][i];
            for j in (i + 1)..4 {
                s -= aug_a[i][j] * result[c][j];
            }
            result[c][i] = s / aug_a[i][i];
        }
    }
    Some(result)
}

// ─────────────────────────────────────────────────────────────
//  Principal stress decomposition (Cardano's method)
// ─────────────────────────────────────────────────────────────

/// Compute the three principal stresses and their directions from a symmetric
/// 3×3 stress tensor using the trigonometric solution of the characteristic cubic.
///
/// Input: `stress = [xx, yy, zz, xy, yz, xz]`
///
/// Returns `(eigenvalues, eigenvectors)` where eigenvalues are sorted descending
/// (σ₁ ≥ σ₂ ≥ σ₃) and eigenvectors[i] is the unit direction for eigenvalue i.
pub fn principal_stresses(stress: &[f64; 6]) -> ([f64; 3], [[f64; 3]; 3]) {
    let sxx = stress[0];
    let syy = stress[1];
    let szz = stress[2];
    let sxy = stress[3];
    let syz = stress[4];
    let sxz = stress[5];

    // Stress matrix:
    // [[sxx, sxy, sxz],
    //  [sxy, syy, syz],
    //  [sxz, syz, szz]]

    // Invariants of the stress tensor
    let i1 = sxx + syy + szz;
    let i2 = sxx * syy + syy * szz + szz * sxx - sxy * sxy - syz * syz - sxz * sxz;
    let i3 = sxx * (syy * szz - syz * syz)
           - sxy * (sxy * szz - syz * sxz)
           + sxz * (sxy * syz - syy * sxz);

    // Solve λ³ - I₁λ² + I₂λ - I₃ = 0
    // Substitute λ = t + I₁/3 to get depressed cubic t³ + pt + q = 0
    let i1_3 = i1 / 3.0;
    let p = i2 - i1 * i1_3;
    let q = -2.0 * i1_3 * i1_3 * i1_3 + i1_3 * i2 - i3;

    let eigenvalues = if p.abs() < 1e-30 {
        // All eigenvalues equal
        [i1_3, i1_3, i1_3]
    } else {
        let r = ((-p / 3.0).max(0.0)).sqrt();
        let cos_arg = (-q / (2.0 * r * r * r)).clamp(-1.0, 1.0);
        let phi = cos_arg.acos();

        let two_r = 2.0 * r;
        let mut vals = [
            two_r * (phi / 3.0).cos() + i1_3,
            two_r * ((phi + 2.0 * std::f64::consts::PI) / 3.0).cos() + i1_3,
            two_r * ((phi + 4.0 * std::f64::consts::PI) / 3.0).cos() + i1_3,
        ];

        // Sort descending
        vals.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
        vals
    };

    // Compute eigenvectors using the cross-product method
    let s = [
        [sxx, sxy, sxz],
        [sxy, syy, syz],
        [sxz, syz, szz],
    ];

    let mut eigenvectors = [[0.0f64; 3]; 3];
    for k in 0..3 {
        eigenvectors[k] = eigenvector_cross_product(&s, eigenvalues[k]);
    }

    // Ensure orthogonality: if v2 is not orthogonal to v1, re-derive via cross product
    let dot01 = dot3(&eigenvectors[0], &eigenvectors[1]);
    if dot01.abs() > 1e-6 {
        eigenvectors[2] = cross3(&eigenvectors[0], &eigenvectors[1]);
        let len = norm3(&eigenvectors[2]);
        if len > 1e-30 {
            for d in 0..3 {
                eigenvectors[2][d] /= len;
            }
        }
        // Re-derive v2 = v3 × v1
        eigenvectors[1] = cross3(&eigenvectors[2], &eigenvectors[0]);
        let len = norm3(&eigenvectors[1]);
        if len > 1e-30 {
            for d in 0..3 {
                eigenvectors[1][d] /= len;
            }
        }
    }

    (eigenvalues, eigenvectors)
}

/// Compute eigenvector for a symmetric 3×3 matrix via cross-product of rows of (S - λI).
fn eigenvector_cross_product(s: &[[f64; 3]; 3], lambda: f64) -> [f64; 3] {
    let m = [
        [s[0][0] - lambda, s[0][1], s[0][2]],
        [s[1][0], s[1][1] - lambda, s[1][2]],
        [s[2][0], s[2][1], s[2][2] - lambda],
    ];

    // Try all three pairs of rows, pick the one with the largest cross product
    let pairs = [(0, 1), (0, 2), (1, 2)];
    let mut best = [0.0f64; 3];
    let mut best_len = 0.0f64;

    for &(i, j) in &pairs {
        let c = cross3(&m[i], &m[j]);
        let len = norm3(&c);
        if len > best_len {
            best_len = len;
            best = c;
        }
    }

    if best_len > 1e-30 {
        let inv = 1.0 / best_len;
        [best[0] * inv, best[1] * inv, best[2] * inv]
    } else {
        // Degenerate: return a canonical basis vector not parallel to any row
        for d in 0..3 {
            let mut e = [0.0; 3];
            e[d] = 1.0;
            let c = cross3(&m[0], &e);
            let len = norm3(&c);
            if len > 1e-15 {
                let inv = 1.0 / len;
                return [c[0] * inv, c[1] * inv, c[2] * inv];
            }
        }
        [1.0, 0.0, 0.0]
    }
}

// ─────────────────────────────────────────────────────────────
//  Von Mises from principal stresses
// ─────────────────────────────────────────────────────────────

/// Compute von Mises stress from three principal stresses.
#[inline]
pub fn von_mises_from_principals(s1: f64, s2: f64, s3: f64) -> f64 {
    (0.5 * ((s1 - s2).powi(2) + (s2 - s3).powi(2) + (s3 - s1).powi(2))).sqrt()
}

// ─────────────────────────────────────────────────────────────
//  Vector helpers
// ─────────────────────────────────────────────────────────────

#[inline]
fn cross3(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[inline]
fn dot3(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

#[inline]
fn norm3(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

// ─────────────────────────────────────────────────────────────
//  Tests
// ─────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_principal_stresses_diagonal() {
        // Pure diagonal stress: eigenvalues should be the diagonal entries
        let stress = [100.0, 50.0, 25.0, 0.0, 0.0, 0.0];
        let (vals, vecs) = principal_stresses(&stress);
        assert!((vals[0] - 100.0).abs() < 1e-10);
        assert!((vals[1] - 50.0).abs() < 1e-10);
        assert!((vals[2] - 25.0).abs() < 1e-10);

        // Eigenvectors should be axis-aligned
        for k in 0..3 {
            let len = norm3(&vecs[k]);
            assert!((len - 1.0).abs() < 1e-10, "eigenvector {k} not unit length");
        }
    }

    #[test]
    fn test_principal_stresses_hydrostatic() {
        // Hydrostatic: all eigenvalues equal
        let p = 42.0;
        let stress = [p, p, p, 0.0, 0.0, 0.0];
        let (vals, _) = principal_stresses(&stress);
        for v in &vals {
            assert!((v - p).abs() < 1e-10);
        }
    }

    #[test]
    fn test_von_mises_uniaxial() {
        // Uniaxial tension: σ₁ = 100, σ₂ = σ₃ = 0 → VM = 100
        let vm = von_mises_from_principals(100.0, 0.0, 0.0);
        assert!((vm - 100.0).abs() < 1e-10);
    }

    #[test]
    fn test_von_mises_hydrostatic() {
        // Pure hydrostatic: VM = 0
        let vm = von_mises_from_principals(50.0, 50.0, 50.0);
        assert!(vm.abs() < 1e-10);
    }

    #[test]
    fn test_spr_single_tet() {
        // Single tet: SPR should return the element stress at all nodes
        let positions = Array2::from_shape_vec(
            (4, 3),
            vec![0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        )
        .unwrap();
        let elements = vec![[0, 1, 2, 3]];
        let stresses = vec![[100.0, 50.0, 25.0, 10.0, 5.0, 2.0]];

        let (nodal, errors) = spr_recover_nodal_stresses(&positions, &elements, &stresses);

        assert_eq!(nodal.len(), 4);
        assert_eq!(errors.len(), 1);

        // With a single element, all nodes share the same patch and the recovered
        // stress should be close to the element stress
        for ni in 0..4 {
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
    fn test_spr_coplanar_centroids_no_blowup() {
        // 7 tets in a fan sharing node 0 at the origin.
        // All centroids lie nearly in the z=0 plane, simulating a boundary
        // surface node. This would cause catastrophic extrapolation with an
        // uncentered polynomial basis.
        let positions = Array2::from_shape_vec(
            (10, 3),
            vec![
                0.0,  0.0,  0.0,   // 0: shared corner
                2.0,  0.0,  0.0,   // 1
                1.0,  2.0,  0.0,   // 2
               -1.0,  2.0,  0.0,   // 3
               -2.0,  0.0,  0.0,   // 4
               -1.0, -2.0,  0.0,   // 5
                1.0, -2.0,  0.0,   // 6
                2.0,  1.0,  0.0,   // 7
                0.0,  0.0,  0.2,   // 8: slightly above z=0
                0.0,  0.0, -0.2,   // 9: slightly below z=0
            ],
        )
        .unwrap();

        let elements = vec![
            [0, 1, 2, 8],
            [0, 2, 3, 8],
            [0, 3, 4, 9],
            [0, 4, 5, 9],
            [0, 5, 6, 8],
            [0, 6, 7, 9],
            [0, 7, 1, 8],
        ];

        let stresses = vec![
            [100.0, 50.0, 25.0, 10.0, 5.0, 2.0],
            [110.0, 55.0, 28.0, 12.0, 6.0, 3.0],
            [ 90.0, 45.0, 22.0,  8.0, 4.0, 1.0],
            [105.0, 52.0, 26.0, 11.0, 5.5, 2.5],
            [ 95.0, 48.0, 24.0,  9.0, 4.5, 1.5],
            [102.0, 51.0, 25.5, 10.5, 5.2, 2.2],
            [ 98.0, 49.0, 24.5,  9.5, 4.8, 1.8],
        ];

        let (nodal, _) = spr_recover_nodal_stresses(&positions, &elements, &stresses);

        // Node 0 has 7 adjacent elements (>= 6), so the linear fit path is
        // attempted. With nearly coplanar centroids the normal matrix may be
        // ill-conditioned; centering + near-singularity detection + sanity
        // check should prevent blowup.
        for c in 0..6 {
            let elem_max = stresses.iter().map(|s| s[c].abs()).fold(0.0f64, f64::max);
            assert!(
                nodal[0][c].abs() < 5.0 * elem_max.max(1.0),
                "node 0, component {c}: recovered {:.1} is unreasonably large (max element: {:.1})",
                nodal[0][c], elem_max,
            );
        }
    }

    #[test]
    fn test_spr_large_coordinates_no_blowup() {
        // Same geometry as coplanar test but translated far from the origin
        // (simulating a model in millimeters). Without centering, the large
        // coordinates would amplify gradient errors catastrophically.
        let offset = 500.0;
        let positions = Array2::from_shape_vec(
            (10, 3),
            vec![
                offset + 0.0,  offset + 0.0,   0.0,
                offset + 2.0,  offset + 0.0,   0.0,
                offset + 1.0,  offset + 2.0,   0.0,
                offset - 1.0,  offset + 2.0,   0.0,
                offset - 2.0,  offset + 0.0,   0.0,
                offset - 1.0,  offset - 2.0,   0.0,
                offset + 1.0,  offset - 2.0,   0.0,
                offset + 2.0,  offset + 1.0,   0.0,
                offset + 0.0,  offset + 0.0,   0.2,
                offset + 0.0,  offset + 0.0,  -0.2,
            ],
        )
        .unwrap();

        let elements = vec![
            [0, 1, 2, 8],
            [0, 2, 3, 8],
            [0, 3, 4, 9],
            [0, 4, 5, 9],
            [0, 5, 6, 8],
            [0, 6, 7, 9],
            [0, 7, 1, 8],
        ];

        let stresses = vec![
            [100.0, 50.0, 25.0, 10.0, 5.0, 2.0],
            [110.0, 55.0, 28.0, 12.0, 6.0, 3.0],
            [ 90.0, 45.0, 22.0,  8.0, 4.0, 1.0],
            [105.0, 52.0, 26.0, 11.0, 5.5, 2.5],
            [ 95.0, 48.0, 24.0,  9.0, 4.5, 1.5],
            [102.0, 51.0, 25.5, 10.5, 5.2, 2.2],
            [ 98.0, 49.0, 24.5,  9.5, 4.8, 1.8],
        ];

        let (nodal, _) = spr_recover_nodal_stresses(&positions, &elements, &stresses);

        for c in 0..6 {
            let elem_max = stresses.iter().map(|s| s[c].abs()).fold(0.0f64, f64::max);
            assert!(
                nodal[0][c].abs() < 5.0 * elem_max.max(1.0),
                "node 0 (offset={offset}), component {c}: recovered {:.1} is unreasonably large (max element: {:.1})",
                nodal[0][c], elem_max,
            );
        }
    }
}
