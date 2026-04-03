//! Forward FEA solver for flat shell triangle elements with 6 DOFs per node.
//! Pipeline: assemble K → build RHS → factorise → solve K u = f → post-process.

use crate::shell_assembly::{assemble_shell_k, assemble_shell_loads, compute_tri_shell_ke_global, tri_local_frame, dkt_bending_b_matrix, shell_element_selfweight};
use crate::shell_types::{ShellCache, ShellProblem, ShellResult};
use crate::types::{Factorization, FactorizationStrategy, TheseusError};
use rayon::prelude::*;

// ─────────────────────────────────────────────────────────────
//  Factor and solve
// ─────────────────────────────────────────────────────────────

fn factor_and_solve_shell(cache: &mut ShellCache) -> Result<(), TheseusError> {
    let mut need_ldl_fallback = false;

    match &mut cache.factorization {
        Some(fac) => {
            if let Err(_e) = fac.update(&cache.k_matrix) {
                if fac.strategy() == FactorizationStrategy::Cholesky {
                    need_ldl_fallback = true;
                } else {
                    return Err(_e.into());
                }
            }
        }
        None => {
            match Factorization::new(&cache.k_matrix, FactorizationStrategy::Cholesky) {
                Ok(fac) => {
                    cache.factorization = Some(fac);
                }
                Err(_e) => {
                    need_ldl_fallback = true;
                }
            }
        }
    }

    if need_ldl_fallback {
        cache.factorization = None;
        cache.factorization =
            Some(Factorization::new(&cache.k_matrix, FactorizationStrategy::LDL)?);
    }

    let fac = cache
        .factorization
        .as_ref()
        .ok_or(TheseusError::MissingFactorization)?;
    let solution = fac.solve(&cache.rhs);

    if solution.iter().any(|v| !v.is_finite()) {
        return Err(TheseusError::Solver(
            "Shell FEA linear solve produced non-finite displacements (singular or ill-conditioned \
             stiffness matrix). Check supports and element connectivity."
                .into(),
        ));
    }

    cache.displacements.copy_from_slice(&solution);
    Ok(())
}

// ─────────────────────────────────────────────────────────────
//  Post-processing: reactions (R = Ke*u - f at fixed DOFs)
// ─────────────────────────────────────────────────────────────

fn post_process(
    cache: &mut ShellCache,
    problem: &ShellProblem,
    u_full: &[f64],
) {
    let ne = problem.num_elements;
    let n_total = cache.reactions.len();

    // Parallel fold/reduce for reaction assembly (same pattern as assemble_shell_k)
    let merged = (0..ne)
        .into_par_iter()
        .fold(
            || vec![0.0f64; n_total],
            |mut buf, e| {
                let nodes = &problem.elements[e];
                if nodes.len() != 3 {
                    return buf;
                }

                let mut u_e = [0.0f64; 18];
                for (i, &ni) in nodes.iter().enumerate() {
                    for d in 0..6 {
                        u_e[i * 6 + d] = u_full[ni * 6 + d];
                    }
                }

                let mut ke_flat = [0.0f64; 18 * 18];
                compute_tri_shell_ke_global(e, problem, &mut ke_flat);

                for li in 0..18 {
                    let ni = nodes[li / 6];
                    let gi = ni * 6 + (li % 6);
                    let mut s = 0.0;
                    for lj in 0..18 {
                        s += ke_flat[li * 18 + lj] * u_e[lj];
                    }
                    buf[gi] += s;
                }
                buf
            },
        )
        .reduce(
            || vec![0.0f64; n_total],
            |mut a, b| {
                for (ai, bi) in a.iter_mut().zip(b.iter()) {
                    *ai += *bi;
                }
                a
            },
        );

    cache.reactions.copy_from_slice(&merged);

    for load in &problem.loads {
        for d in 0..3 {
            cache.reactions[load.node_idx * 6 + d] -= load.force[d];
        }
        for d in 0..3 {
            cache.reactions[load.node_idx * 6 + 3 + d] -= load.moment[d];
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
                    cache.reactions[ni * 6 + d] -= forces[node_local][d];
                }
            }
        }
    }

    for &gi in &cache.dof_map.free_dofs {
        cache.reactions[gi] = 0.0;
    }
}

// ─────────────────────────────────────────────────────────────
//  Membrane stress recovery (CST principal stresses)
// ─────────────────────────────────────────────────────────────

fn compute_principal_stresses(problem: &ShellProblem, u_full: &[f64]) -> Vec<[f64; 2]> {
    let ne = problem.num_elements;

    (0..ne)
        .into_par_iter()
        .map(|e| {
            let nodes = &problem.elements[e];
            if nodes.len() != 3 {
                return [0.0, 0.0];
            }

            let p1 = problem.node_positions.row(nodes[0]);
            let p2 = problem.node_positions.row(nodes[1]);
            let p3 = problem.node_positions.row(nodes[2]);

            let (ex, ey, _ez, lc, area) = tri_local_frame(
                p1.as_slice().unwrap(),
                p2.as_slice().unwrap(),
                p3.as_slice().unwrap(),
            );

            let (_x1, _y1, x2, _y2, x3, y3) = (lc[0], lc[1], lc[2], lc[3], lc[4], lc[5]);

            let b_mem = [
                [-y3, 0.0, y3, 0.0, 0.0, 0.0],
                [0.0, x3 - x2, 0.0, -x3, 0.0, x2],
                [x3 - x2, -y3, -x3, y3, x2, 0.0],
            ];
            let inv_2a = 1.0 / (2.0 * area);

            let mut u_loc = [0.0f64; 6];
            for n in 0..3 {
                let ug = &u_full[nodes[n] * 6..nodes[n] * 6 + 3];
                u_loc[n * 2] = ug[0] * ex[0] + ug[1] * ex[1] + ug[2] * ex[2];
                u_loc[n * 2 + 1] = ug[0] * ey[0] + ug[1] * ey[1] + ug[2] * ey[2];
            }

            let mut eps = [0.0f64; 3];
            for i in 0..3 {
                for j in 0..6 {
                    eps[i] += b_mem[i][j] * u_loc[j] * inv_2a;
                }
            }

            let props = &problem.element_props[e];
            let mat = &problem.materials[props.material_idx];
            let e_mod = mat.e;
            let nu = mat.nu;
            let c = e_mod / (1.0 - nu * nu);
            let s_xx = c * (eps[0] + nu * eps[1]);
            let s_yy = c * (nu * eps[0] + eps[1]);
            let s_xy = (e_mod / (2.0 * (1.0 + nu))) * eps[2];

            let center = (s_xx + s_yy) / 2.0;
            let radius = (((s_xx - s_yy) / 2.0).powi(2) + s_xy.powi(2)).sqrt();
            [center + radius, center - radius]
        })
        .collect()
}

// ─────────────────────────────────────────────────────────────
//  Von Mises from 6-component Voigt stress tensor
// ─────────────────────────────────────────────────────────────

#[inline]
fn von_mises_6(s: &[f64; 6]) -> f64 {
    let (sxx, syy, szz, sxy, syz, sxz) = (s[0], s[1], s[2], s[3], s[4], s[5]);
    (0.5 * ((sxx - syy).powi(2) + (syy - szz).powi(2) + (szz - sxx).powi(2))
        + 3.0 * (sxy * sxy + syz * syz + sxz * sxz))
    .sqrt()
}

// ─────────────────────────────────────────────────────────────
//  Rotate local in-plane stress [sxx, syy, sxy] to global 6-component Voigt
// ─────────────────────────────────────────────────────────────

/// Transform a 2D in-plane stress tensor from element-local frame (ex,ey) to
/// a 6-component symmetric 3D Voigt tensor [XX, YY, ZZ, XY, YZ, XZ] in global coordinates.
#[inline]
fn local_stress_to_global(
    sxx: f64, syy: f64, sxy: f64,
    ex: &[f64; 3], ey: &[f64; 3],
) -> [f64; 6] {
    [
        ex[0]*ex[0]*sxx + ey[0]*ey[0]*syy + 2.0*ex[0]*ey[0]*sxy,
        ex[1]*ex[1]*sxx + ey[1]*ey[1]*syy + 2.0*ex[1]*ey[1]*sxy,
        ex[2]*ex[2]*sxx + ey[2]*ey[2]*syy + 2.0*ex[2]*ey[2]*sxy,
        ex[0]*ex[1]*sxx + ey[0]*ey[1]*syy + (ex[0]*ey[1]+ey[0]*ex[1])*sxy,
        ex[1]*ex[2]*sxx + ey[1]*ey[2]*syy + (ex[1]*ey[2]+ey[1]*ex[2])*sxy,
        ex[0]*ex[2]*sxx + ey[0]*ey[2]*syy + (ex[0]*ey[2]+ey[0]*ex[2])*sxy,
    ]
}

// ─────────────────────────────────────────────────────────────
//  Combined membrane + bending stress recovery in global frame
// ─────────────────────────────────────────────────────────────

/// Per-element stress output: membrane, top-fiber, and bottom-fiber in global Voigt notation.
pub struct ShellElementStresses {
    pub membrane: Vec<[f64; 6]>,
    pub top: Vec<[f64; 6]>,
    pub bottom: Vec<[f64; 6]>,
}

/// Compute per-element membrane, top-fiber, and bottom-fiber stress tensors
/// in global 6-component Voigt notation [XX, YY, ZZ, XY, YZ, XZ].
///
/// Membrane stresses come from the CST formulation.
/// Bending stresses are evaluated at the element centroid using DKT curvatures:
///   sigma_bend(z) = C_plane * z * kappa
/// Top fiber: z = +t/2 + offset, bottom fiber: z = -t/2 + offset.
pub fn compute_shell_element_stresses_global(
    problem: &ShellProblem, u_full: &[f64],
) -> ShellElementStresses {
    let ne = problem.num_elements;

    let results: Vec<([f64; 6], [f64; 6], [f64; 6])> = (0..ne)
        .into_par_iter()
        .map(|e| {
            let nodes = &problem.elements[e];
            if nodes.len() != 3 {
                return ([0.0; 6], [0.0; 6], [0.0; 6]);
            }

            let p1 = problem.node_positions.row(nodes[0]);
            let p2 = problem.node_positions.row(nodes[1]);
            let p3 = problem.node_positions.row(nodes[2]);

            let (ex, ey, ez, lc, area) = tri_local_frame(
                p1.as_slice().unwrap(),
                p2.as_slice().unwrap(),
                p3.as_slice().unwrap(),
            );

            let props = &problem.element_props[e];
            let mat = &problem.materials[props.material_idx];
            let e_mod = mat.e;
            let nu = mat.nu;

            let t_avg = (problem.node_thicknesses[nodes[0]]
                + problem.node_thicknesses[nodes[1]]
                + problem.node_thicknesses[nodes[2]])
                / 3.0;

            let section_offset = if props.section_idx < problem.sections.len() {
                problem.sections[props.section_idx].offset
            } else {
                0.0
            };

            // ── Plane-stress constitutive matrix C_plane ──
            let c11 = e_mod / (1.0 - nu * nu);
            let c12 = nu * c11;
            let c33 = e_mod / (2.0 * (1.0 + nu));

            // ── Membrane stress (CST) ──
            let (_x1, _y1, x2, _y2, x3, y3) = (lc[0], lc[1], lc[2], lc[3], lc[4], lc[5]);
            let b_mem = [
                [-y3, 0.0, y3, 0.0, 0.0, 0.0],
                [0.0, x3 - x2, 0.0, -x3, 0.0, x2],
                [x3 - x2, -y3, -x3, y3, x2, 0.0],
            ];
            let inv_2a = 1.0 / (2.0 * area);

            let mut u_mem_loc = [0.0f64; 6];
            for n in 0..3 {
                let ug = &u_full[nodes[n] * 6..nodes[n] * 6 + 3];
                u_mem_loc[n * 2] = ug[0] * ex[0] + ug[1] * ex[1] + ug[2] * ex[2];
                u_mem_loc[n * 2 + 1] = ug[0] * ey[0] + ug[1] * ey[1] + ug[2] * ey[2];
            }

            let mut eps_mem = [0.0f64; 3];
            for i in 0..3 {
                for j in 0..6 {
                    eps_mem[i] += b_mem[i][j] * u_mem_loc[j] * inv_2a;
                }
            }

            let s_mem_xx = c11 * eps_mem[0] + c12 * eps_mem[1];
            let s_mem_yy = c12 * eps_mem[0] + c11 * eps_mem[1];
            let s_mem_xy = c33 * eps_mem[2];

            // ── Bending curvatures (DKT at centroid) ──
            let two_a = 2.0 * area;
            let bb = dkt_bending_b_matrix(&lc, two_a, 1.0 / 3.0, 1.0 / 3.0);

            let mut u_bend_loc = [0.0f64; 9];
            for n in 0..3 {
                let ut = &u_full[nodes[n] * 6..nodes[n] * 6 + 3];
                let ur = &u_full[nodes[n] * 6 + 3..nodes[n] * 6 + 6];
                // w = translation along element normal
                u_bend_loc[n * 3] = ut[0] * ez[0] + ut[1] * ez[1] + ut[2] * ez[2];
                // thx = rotation about local x
                u_bend_loc[n * 3 + 1] = ur[0] * ex[0] + ur[1] * ex[1] + ur[2] * ex[2];
                // thy = rotation about local y
                u_bend_loc[n * 3 + 2] = ur[0] * ey[0] + ur[1] * ey[1] + ur[2] * ey[2];
            }

            let mut kappa = [0.0f64; 3];
            for i in 0..3 {
                for j in 0..9 {
                    kappa[i] += bb[i][j] * u_bend_loc[j];
                }
            }

            // sigma_bend(z) = C_plane * z * kappa
            let z_top = t_avg / 2.0 + section_offset;
            let z_bot = -t_avg / 2.0 + section_offset;

            let s_bend_top = [
                (c11 * kappa[0] + c12 * kappa[1]) * z_top,
                (c12 * kappa[0] + c11 * kappa[1]) * z_top,
                c33 * kappa[2] * z_top,
            ];
            let s_bend_bot = [
                (c11 * kappa[0] + c12 * kappa[1]) * z_bot,
                (c12 * kappa[0] + c11 * kappa[1]) * z_bot,
                c33 * kappa[2] * z_bot,
            ];

            // ── Combined top/bottom stresses in local frame ──
            let s_top_xx = s_mem_xx + s_bend_top[0];
            let s_top_yy = s_mem_yy + s_bend_top[1];
            let s_top_xy = s_mem_xy + s_bend_top[2];

            let s_bot_xx = s_mem_xx + s_bend_bot[0];
            let s_bot_yy = s_mem_yy + s_bend_bot[1];
            let s_bot_xy = s_mem_xy + s_bend_bot[2];

            // ── Rotate to global Voigt notation ──
            let g_mem = local_stress_to_global(s_mem_xx, s_mem_yy, s_mem_xy, &ex, &ey);
            let g_top = local_stress_to_global(s_top_xx, s_top_yy, s_top_xy, &ex, &ey);
            let g_bot = local_stress_to_global(s_bot_xx, s_bot_yy, s_bot_xy, &ex, &ey);

            (g_mem, g_top, g_bot)
        })
        .collect();

    let mut membrane = Vec::with_capacity(ne);
    let mut top = Vec::with_capacity(ne);
    let mut bottom = Vec::with_capacity(ne);
    for (m, t, b) in results {
        membrane.push(m);
        top.push(t);
        bottom.push(b);
    }

    ShellElementStresses { membrane, top, bottom }
}

// ─────────────────────────────────────────────────────────────
//  Top-level forward solve
// ─────────────────────────────────────────────────────────────

/// Full forward shell FEA solve. Returns a ShellResult.
pub fn solve_shell(
    cache: &mut ShellCache,
    problem: &ShellProblem,
) -> Result<ShellResult, TheseusError> {
    assemble_shell_k(cache, problem)?;
    assemble_shell_loads(cache, problem);
    factor_and_solve_shell(cache)?;

    let n_total = cache.dof_map.num_total_dofs;
    let mut u_full = vec![0.0; n_total];
    for (fi, &gi) in cache.dof_map.free_dofs.iter().enumerate() {
        u_full[gi] = cache.displacements[fi];
    }

    post_process(cache, problem, &u_full);

    let principal_stresses = compute_principal_stresses(problem, &u_full);
    let elem_stresses = compute_shell_element_stresses_global(problem, &u_full);

    let ne = problem.num_elements;
    let von_mises: Vec<f64> = (0..ne)
        .into_par_iter()
        .map(|e| {
            let vm_top = von_mises_6(&elem_stresses.top[e]);
            let vm_bot = von_mises_6(&elem_stresses.bottom[e]);
            vm_top.max(vm_bot)
        })
        .collect();

    Ok(ShellResult {
        displacements: u_full,
        reactions: cache.reactions.clone(),
        principal_stresses,
        membrane_stresses_global: elem_stresses.membrane,
        top_stresses_global: elem_stresses.top,
        bottom_stresses_global: elem_stresses.bottom,
        von_mises,
    })
}
