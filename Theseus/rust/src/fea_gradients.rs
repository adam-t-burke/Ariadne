//! Hand-coded gradients for FEA objectives + adjoint solve.
//!
//! The adjoint method for FEA:
//!   1. Forward solve: K(θ) u = f(θ)
//!   2. Accumulate dJ/du from each objective
//!   3. Adjoint solve: K λ = dJ/du  (reuse K factorization, K symmetric)
//!   4. Implicit gradient: dJ/dθ = -λ^T (dK/dθ) u + (df/dθ)^T λ + dJ/dθ_explicit
//!
//! Design variables θ can include node positions, cross-section areas,
//! and support positions.

use crate::fea_assembly::bar_element_stiffness;
use crate::fea_objectives::{fea_bounds_penalty, fea_bounds_penalty_grad, fea_total_loss};
use crate::fea_solve::solve_fea;
use crate::fea_types::{FeaCache, FeaProblem, FeaSnapshot, FeaVariableMask};
use crate::types::TheseusError;
use rayon::prelude::*;

// ─────────────────────────────────────────────────────────────
//  Adjoint solve
// ─────────────────────────────────────────────────────────────

/// Solve K λ = dJ/du using the factorization from the forward solve.
fn solve_adjoint(cache: &mut FeaCache) -> Result<(), TheseusError> {
    let fac = cache.factorization.as_ref()
        .ok_or(TheseusError::MissingFactorization)?;
    let solution = fac.solve(&cache.grad_u);
    cache.lambda.copy_from_slice(&solution);
    Ok(())
}

// ─────────────────────────────────────────────────────────────
//  Owned snapshot data (avoids borrow conflicts with cache)
// ─────────────────────────────────────────────────────────────

struct CacheSnapshot {
    axial_forces: Vec<f64>,
    stresses: Vec<f64>,
    strains: Vec<f64>,
    reactions: Vec<f64>,
    utilization: Vec<f64>,
    elem_lengths: Vec<f64>,
    elem_cos: Vec<[f64; 3]>,
}

impl CacheSnapshot {
    fn from_cache(cache: &FeaCache) -> Self {
        Self {
            axial_forces: cache.axial_forces.clone(),
            stresses: cache.stresses.clone(),
            strains: cache.strains.clone(),
            reactions: cache.reactions.clone(),
            utilization: cache.utilization.clone(),
            elem_lengths: cache.elem_lengths.clone(),
            elem_cos: cache.elem_cos.clone(),
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Implicit gradient: -λ^T (dK/d(node_pos)) u
// ─────────────────────────────────────────────────────────────

/// For bar element e connecting nodes i,j:
///   ke = (EA/L) * [T -T; -T T] where T = c*c^T, c = (xj-xi)/L
///
/// Derivatives w.r.t. xi_d (node i, direction d):
///   dL/d(xi_d) = -c_d
///   d(EA/L)/d(xi_d) = EA*c_d/L^2
///   dc_a/d(xi_d) = -(δ_{ad} - c_a*c_d) / L
///
/// For xj_d: all signs flip relative to xi_d.
fn accumulate_dk_dpos(
    grad_pos: &mut [f64],
    cache: &FeaCache,
    problem: &FeaProblem,
    areas: &[f64],
    u_full: &[f64],
    lambda_full: &[f64],
) {
    let ne = problem.num_elements;
    let n_dofs = grad_pos.len();

    let partial = (0..ne)
        .into_par_iter()
        .fold(
            || vec![0.0; n_dofs],
            |mut local_grad, e| {
                let (ni, nj) = problem.edge_nodes[e];
                let props = &problem.element_props[e];
                let mat = &problem.materials[props.material_idx];
                let a = areas[props.section_idx];
                let l = cache.elem_lengths[e];
                if l < f64::EPSILON { return local_grad; }
                let c = cache.elem_cos[e];
                let ea = mat.e * a;
                let ea_over_l = ea / l;

                let u_e: [f64; 6] = [
                    u_full[ni * 3], u_full[ni * 3 + 1], u_full[ni * 3 + 2],
                    u_full[nj * 3], u_full[nj * 3 + 1], u_full[nj * 3 + 2],
                ];
                let lam_e: [f64; 6] = [
                    lambda_full[ni * 3], lambda_full[ni * 3 + 1], lambda_full[ni * 3 + 2],
                    lambda_full[nj * 3], lambda_full[nj * 3 + 1], lambda_full[nj * 3 + 2],
                ];

                for d in 0..3 {
                    let d_ea_over_l_dxi = ea * c[d] / (l * l);
                    let d_ea_over_l_dxj = -ea * c[d] / (l * l);

                    let mut lam_dke_u_xi = 0.0;
                    let mut lam_dke_u_xj = 0.0;

                    for a_idx in 0..6 {
                        for b_idx in 0..6 {
                            let sa = if a_idx < 3 { 1.0 } else { -1.0 };
                            let sb = if b_idx < 3 { 1.0 } else { -1.0 };
                            let block_sign = sa * sb;

                            let ca = c[a_idx % 3];
                            let cb = c[b_idx % 3];

                            let t1_xi = d_ea_over_l_dxi * block_sign * ca * cb;
                            let t1_xj = d_ea_over_l_dxj * block_sign * ca * cb;

                            let ad = a_idx % 3;
                            let bd = b_idx % 3;
                            let delta_ad = if ad == d { 1.0 } else { 0.0 };
                            let delta_bd = if bd == d { 1.0 } else { 0.0 };
                            let dca_dxi = -(delta_ad - ca * c[d]) / l;
                            let dcb_dxi = -(delta_bd - cb * c[d]) / l;
                            let dt_ab_dxi = dca_dxi * cb + ca * dcb_dxi;
                            let dt_ab_dxj = -dt_ab_dxi;

                            let t2_xi = ea_over_l * block_sign * dt_ab_dxi;
                            let t2_xj = ea_over_l * block_sign * dt_ab_dxj;

                            lam_dke_u_xi += lam_e[a_idx] * (t1_xi + t2_xi) * u_e[b_idx];
                            lam_dke_u_xj += lam_e[a_idx] * (t1_xj + t2_xj) * u_e[b_idx];
                        }
                    }

                    local_grad[ni * 3 + d] -= lam_dke_u_xi;
                    local_grad[nj * 3 + d] -= lam_dke_u_xj;
                }
                local_grad
            },
        )
        .reduce(
            || vec![0.0; n_dofs],
            |mut a, b| {
                for (ai, bi) in a.iter_mut().zip(b.iter()) {
                    *ai += *bi;
                }
                a
            },
        );

    for (g, p) in grad_pos.iter_mut().zip(partial.iter()) {
        *g += *p;
    }
}

// ─────────────────────────────────────────────────────────────
//  Implicit gradient: -λ^T (dK/dA) u
// ─────────────────────────────────────────────────────────────

/// For bar: ke = (EA/L) [T -T; -T T], so dke/dA = (E/L) [T -T; -T T].
fn accumulate_dk_da(
    grad_area: &mut [f64],
    cache: &FeaCache,
    problem: &FeaProblem,
    areas: &[f64],
    u_full: &[f64],
    lambda_full: &[f64],
) {
    let ne = problem.num_elements;
    for e in 0..ne {
        let (ni, nj) = problem.edge_nodes[e];
        let props = &problem.element_props[e];
        let mat = &problem.materials[props.material_idx];
        let a = areas[props.section_idx];
        let l = cache.elem_lengths[e];
        if l < f64::EPSILON || a < f64::EPSILON { continue; }
        let c = cache.elem_cos[e];
        let e_over_l = mat.e / l;

        let dke = bar_element_stiffness(e_over_l, &c);

        let u_e: [f64; 6] = [
            u_full[ni * 3], u_full[ni * 3 + 1], u_full[ni * 3 + 2],
            u_full[nj * 3], u_full[nj * 3 + 1], u_full[nj * 3 + 2],
        ];
        let lam_e: [f64; 6] = [
            lambda_full[ni * 3], lambda_full[ni * 3 + 1], lambda_full[ni * 3 + 2],
            lambda_full[nj * 3], lambda_full[nj * 3 + 1], lambda_full[nj * 3 + 2],
        ];

        let mut lam_dke_u = 0.0;
        for i in 0..6 {
            for j in 0..6 {
                lam_dke_u += lam_e[i] * dke[i][j] * u_e[j];
            }
        }

        grad_area[props.section_idx] -= lam_dke_u;
    }
}

// ─────────────────────────────────────────────────────────────
//  Self-weight gradient: (df/dpos)^T λ
// ─────────────────────────────────────────────────────────────

/// Self-weight load at each node of element e:
///   f_{ni,d_grav} += 0.5 * ρ * A * L * g[d_grav]
///   f_{nj,d_grav} += 0.5 * ρ * A * L * g[d_grav]
///
/// Differentiating w.r.t. xi_{d_pos}:
///   dL/d(xi_{d_pos}) = -c[d_pos]
///   df_{ni,d_grav}/d(xi_{d_pos}) = 0.5 * ρ * A * (-c[d_pos]) * g[d_grav]
///   df_{nj,d_grav}/d(xi_{d_pos}) = 0.5 * ρ * A * (-c[d_pos]) * g[d_grav]
///   (both nodes get half the element weight, both change by the same amount)
///
/// The (df/dpos)^T λ contribution:
///   grad_pos[ni*3+d_pos] += df_dxi * (λ[ni*3+d_grav] + λ[nj*3+d_grav])
///   grad_pos[nj*3+d_pos] += df_dxj * (λ[ni*3+d_grav] + λ[nj*3+d_grav])
fn accumulate_df_dpos(
    grad_pos: &mut [f64],
    cache: &FeaCache,
    problem: &FeaProblem,
    areas: &[f64],
    lambda_full: &[f64],
) {
    if !problem.include_self_weight { return; }

    let ne = problem.num_elements;
    for e in 0..ne {
        let (ni, nj) = problem.edge_nodes[e];
        let props = &problem.element_props[e];
        let mat = &problem.materials[props.material_idx];
        let a = areas[props.section_idx];
        let l = cache.elem_lengths[e];
        if l < f64::EPSILON { continue; }
        let c = cache.elem_cos[e];
        let rho_a_half = mat.density * a * 0.5;

        for d_pos in 0..3 {
            // dL/d(xi_{d_pos}) = -c[d_pos], dL/d(xj_{d_pos}) = +c[d_pos]
            let dl_dxi = -c[d_pos];
            let dl_dxj = c[d_pos];

            for d_grav in 0..3 {
                let g = problem.gravity[d_grav];
                let df_dxi = rho_a_half * dl_dxi * g;
                let df_dxj = rho_a_half * dl_dxj * g;

                // Both f_ni and f_nj change when we move xi or xj
                let lam_sum = lambda_full[ni * 3 + d_grav] + lambda_full[nj * 3 + d_grav];
                grad_pos[ni * 3 + d_pos] += df_dxi * lam_sum;
                grad_pos[nj * 3 + d_pos] += df_dxj * lam_sum;
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Combined value_and_gradient
// ─────────────────────────────────────────────────────────────

/// Evaluate total loss and gradient for the FEA optimization problem.
///
/// theta layout: [node_pos (3*n_var_nodes), areas (n_var_sections), support_pos (3*n_var_supports)]
pub fn fea_value_and_gradient(
    cache: &mut FeaCache,
    problem: &FeaProblem,
    theta: &[f64],
    grad: &mut [f64],
    mask: &FeaVariableMask,
    lb: &[f64],
    ub: &[f64],
    lb_idx: &[usize],
    ub_idx: &[usize],
) -> Result<f64, TheseusError> {
    // 1. Unpack theta into positions and areas
    let mut positions = problem.node_positions.clone();
    let mut areas: Vec<f64> = problem.sections.iter().map(|s| s.area).collect();
    let mut offset = 0;

    if mask.node_positions {
        for &ni in &mask.node_indices {
            for d in 0..3 {
                positions[[ni, d]] = theta[offset];
                offset += 1;
            }
        }
    }
    if mask.cross_section_areas {
        for &ei in &mask.area_element_indices {
            areas[ei] = theta[offset];
            offset += 1;
        }
    }
    if mask.support_positions {
        for &ni in &mask.support_node_indices {
            for d in 0..3 {
                positions[[ni, d]] = theta[offset];
                offset += 1;
            }
        }
    }

    // 2. Forward solve
    let result = solve_fea(cache, problem, &positions, &areas)?;

    // 3. Build full displacement vector
    let n_total = cache.dof_map.num_total_dofs;
    let mut u_full = vec![0.0; n_total];
    for (fi, &gi) in cache.dof_map.free_dofs.iter().enumerate() {
        u_full[gi] = cache.displacements[fi];
    }

    // 4. Build per-element auxiliary arrays for snapshot
    let ne = problem.num_elements;
    let densities: Vec<f64> = (0..ne)
        .into_par_iter()
        .map(|e| problem.materials[problem.element_props[e].material_idx].density)
        .collect();
    let section_indices: Vec<usize> = (0..ne)
        .into_par_iter()
        .map(|e| problem.element_props[e].section_idx)
        .collect();
    let materials_e: Vec<f64> = (0..ne)
        .into_par_iter()
        .map(|e| problem.materials[problem.element_props[e].material_idx].e)
        .collect();

    // 5. Extract owned data from cache to avoid borrow conflicts
    let snap_data = CacheSnapshot::from_cache(cache);
    let deformed = result.deformed_positions;

    // 6. Build snapshot, evaluate loss, accumulate dJ/du — all before mutating cache
    let snap = FeaSnapshot {
        displacements: &u_full,
        deformed_xyz: &deformed,
        axial_forces: &snap_data.axial_forces,
        stresses: &snap_data.stresses,
        strains: &snap_data.strains,
        reactions: &snap_data.reactions,
        utilization: &snap_data.utilization,
        elem_lengths: &snap_data.elem_lengths,
        node_positions: &positions,
        areas: &areas,
        densities: &densities,
        section_indices: &section_indices,
        edge_nodes: &problem.edge_nodes,
        elem_cos: &snap_data.elem_cos,
        materials_e: &materials_e,
    };

    let obj_loss = fea_total_loss(&problem.objectives, &snap);
    let barrier_loss = fea_bounds_penalty(
        theta, lb, ub, lb_idx, ub_idx, problem.solver.barrier_sharpness,
    );
    let total_loss = obj_loss + problem.solver.barrier_weight * barrier_loss;

    // Accumulate dJ/du from objectives into a temporary buffer (full-DOF space)
    let mut grad_u_full = vec![0.0; n_total];
    for obj in &problem.objectives {
        obj.accumulate_grad_u(&mut grad_u_full, &snap);
    }

    // Compliance special case: J = w * u^T K u.
    //
    // The correct dJ/du to feed the adjoint for compliance is w*f (NOT 2w*f).
    //
    // Derivation: dJ/dθ = 2w u^T K (du/dθ) + w u^T (dK/dθ) u.
    // The adjoint method gives dJ/dθ = -λ^T (dK/dθ) u + λ^T (df/dθ).
    // Setting K λ = w*f gives λ = w*u, so -λ^T (dK/dθ) u = -w u^T (dK/dθ) u,
    // and λ^T (df/dθ) = w u^T (df/dθ). The total matches the true derivative
    // when f doesn't depend on θ: dJ/dθ = -w u^T (dK/dθ) u.
    // When f depends on θ (self-weight): the (df/dθ)^T λ term handles it.
    for obj in &problem.objectives {
        if obj.is_compliance() {
            let w = obj.weight();
            for (fi, &gi) in cache.dof_map.free_dofs.iter().enumerate() {
                grad_u_full[gi] += w * cache.rhs[fi];
            }
        }
    }

    // Accumulate explicit objective gradients into temporary buffers
    let mut grad_pos_explicit = vec![0.0; n_total];
    let mut grad_area_explicit = vec![0.0; problem.sections.len()];
    for obj in &problem.objectives {
        if obj.has_position_gradient() {
            obj.accumulate_grad_pos(&mut grad_pos_explicit, &snap);
        }
        if obj.has_area_gradient() {
            obj.accumulate_grad_area(&mut grad_area_explicit, &snap);
        }
    }

    // Done with snapshot — now we can mutate cache
    drop(snap);

    // 7. Project grad_u to free DOFs and write into cache
    for v in cache.grad_u.iter_mut() { *v = 0.0; }
    for (fi, &gi) in cache.dof_map.free_dofs.iter().enumerate() {
        cache.grad_u[fi] = grad_u_full[gi];
    }

    // 8. Adjoint solve: K λ = dJ/du
    solve_adjoint(cache)?;

    // 9. Build full lambda vector
    let mut lambda_full = vec![0.0; n_total];
    for (fi, &gi) in cache.dof_map.free_dofs.iter().enumerate() {
        lambda_full[gi] = cache.lambda[fi];
    }

    // 10. Accumulate implicit + explicit gradients
    for v in grad.iter_mut() { *v = 0.0; }

    let mut grad_pos = grad_pos_explicit;
    let mut grad_area = grad_area_explicit;

    if mask.node_positions || mask.support_positions {
        accumulate_dk_dpos(&mut grad_pos, cache, problem, &areas, &u_full, &lambda_full);
        accumulate_df_dpos(&mut grad_pos, cache, problem, &areas, &lambda_full);
    }
    if mask.cross_section_areas {
        accumulate_dk_da(&mut grad_area, cache, problem, &areas, &u_full, &lambda_full);
    }

    // 11. Pack gradients back into theta-space
    offset = 0;
    if mask.node_positions {
        for &ni in &mask.node_indices {
            for d in 0..3 {
                grad[offset] = grad_pos[ni * 3 + d];
                offset += 1;
            }
        }
    }
    if mask.cross_section_areas {
        for &ei in &mask.area_element_indices {
            grad[offset] = grad_area[ei];
            offset += 1;
        }
    }
    if mask.support_positions {
        for &ni in &mask.support_node_indices {
            for d in 0..3 {
                grad[offset] = grad_pos[ni * 3 + d];
                offset += 1;
            }
        }
    }

    // 12. Barrier gradient
    fea_bounds_penalty_grad(
        grad, theta, lb, ub, lb_idx, ub_idx,
        problem.solver.barrier_sharpness, problem.solver.barrier_weight,
    );

    Ok(total_loss)
}
