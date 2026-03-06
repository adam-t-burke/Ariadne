//! Forward FEA solver: assemble K, build RHS, factorise, solve K u = f.
//!
//! Post-processing: axial forces, stresses, strains, reactions, utilization.

use crate::fea_assembly::{assemble_global_k, assemble_loads, bar_element_stiffness, update_element_geometry};
use crate::fea_types::{FeaCache, FeaProblem, FeaResult};
use crate::types::{Factorization, FactorizationStrategy, TheseusError};
use ndarray::Array2;
use rayon::prelude::*;

// ─────────────────────────────────────────────────────────────
//  Factor and solve
// ─────────────────────────────────────────────────────────────

/// Factor K and solve K u = f.  Uses Cholesky (K is SPD for well-constrained trusses).
/// Falls back to LDL if Cholesky fails.
fn factor_and_solve(cache: &mut FeaCache) -> Result<(), TheseusError> {
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
        cache.factorization = Some(Factorization::new(&cache.k_matrix, FactorizationStrategy::LDL)?);
    }

    let fac = cache.factorization.as_ref()
        .ok_or(TheseusError::MissingFactorization)?;
    let solution = fac.solve(&cache.rhs);

    if solution.iter().any(|v| !v.is_finite()) {
        return Err(TheseusError::Solver(
            "FEA linear solve produced non-finite displacements (singular or ill-conditioned stiffness matrix). \
             Check supports and element connectivity.".into(),
        ));
    }

    cache.displacements.copy_from_slice(&solution);
    Ok(())
}

// ─────────────────────────────────────────────────────────────
//  Post-processing
// ─────────────────────────────────────────────────────────────

/// Compute axial forces, stresses, strains, reactions, and utilization
/// from the displacement solution.
fn post_process(
    cache: &mut FeaCache,
    problem: &FeaProblem,
    areas: &[f64],
    _positions: &Array2<f64>,
) -> Result<(), TheseusError> {
    let ne = problem.num_elements;
    let n_total = cache.dof_map.num_total_dofs;

    // Expand displacements to full DOF vector (constrained DOFs = 0)
    let mut u_full = vec![0.0; n_total];
    for (fi, &gi) in cache.dof_map.free_dofs.iter().enumerate() {
        u_full[gi] = cache.displacements[fi];
    }

    // Per-element: elongation, axial force, stress, strain, utilization (parallel)
    cache.axial_forces[..ne]
        .par_iter_mut()
        .zip(cache.stresses[..ne].par_iter_mut())
        .zip(cache.strains[..ne].par_iter_mut())
        .zip(cache.utilization[..ne].par_iter_mut())
        .enumerate()
        .for_each(|(e, (((af, st), sn), ut))| {
            let (ni, nj) = problem.edge_nodes[e];
            let props = &problem.element_props[e];
            let mat = &problem.materials[props.material_idx];
            let a = areas[props.section_idx];
            let l = cache.elem_lengths[e];
            let c = &cache.elem_cos[e];

            if l < f64::EPSILON {
                *af = 0.0;
                *st = 0.0;
                *sn = 0.0;
                *ut = 0.0;
                return;
            }

            let mut delta = 0.0;
            for d in 0..3 {
                delta += c[d] * (u_full[nj * 3 + d] - u_full[ni * 3 + d]);
            }

            let strain = delta / l;
            let stress = mat.e * strain;
            let force = stress * a;

            *sn = strain;
            *st = stress;
            *af = force;
            *ut = if mat.yield_stress > 0.0 {
                stress.abs() / mat.yield_stress
            } else {
                0.0
            };
        });

    // Reactions: R = K_full * u - f
    // Parallel fold/reduce with per-thread buffers to avoid data races.
    let n_reactions = cache.reactions.len();

    let reactions = (0..ne)
        .into_par_iter()
        .fold(
            || vec![0.0; n_reactions],
            |mut buf, e| {
                let (ni, nj) = problem.edge_nodes[e];
                let props = &problem.element_props[e];
                let mat = &problem.materials[props.material_idx];
                let a = areas[props.section_idx];
                let l = cache.elem_lengths[e];
                if l < f64::EPSILON { return buf; }
                let c = &cache.elem_cos[e];
                let ea_over_l = mat.e * a / l;
                let ke = bar_element_stiffness(ea_over_l, c);

                let global_dofs = [ni * 3, ni * 3 + 1, ni * 3 + 2,
                                   nj * 3, nj * 3 + 1, nj * 3 + 2];

                for (li, &gi) in global_dofs.iter().enumerate() {
                    let mut ke_u = 0.0;
                    for (lj, &gj) in global_dofs.iter().enumerate() {
                        ke_u += ke[li][lj] * u_full[gj];
                    }
                    buf[gi] += ke_u;
                }
                buf
            },
        )
        .reduce(
            || vec![0.0; n_reactions],
            |mut a_buf, b_buf| {
                for (a, b) in a_buf.iter_mut().zip(b_buf.iter()) {
                    *a += *b;
                }
                a_buf
            },
        );

    cache.reactions.copy_from_slice(&reactions);

    // Subtract applied loads from reactions
    for load in &problem.loads {
        for d in 0..3 {
            cache.reactions[load.node_idx * 3 + d] -= load.force[d];
        }
    }
    if problem.include_self_weight {
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
                cache.reactions[ni * 3 + d] -= f;
                cache.reactions[nj * 3 + d] -= f;
            }
        }
    }

    // Zero out reactions at free DOFs (only constrained DOFs have meaningful reactions)
    for &gi in &cache.dof_map.free_dofs {
        cache.reactions[gi] = 0.0;
    }

    Ok(())
}

// ─────────────────────────────────────────────────────────────
//  Top-level forward solve
// ─────────────────────────────────────────────────────────────

/// Full forward FEA solve. Returns an FeaResult.
pub fn solve_fea(
    cache: &mut FeaCache,
    problem: &FeaProblem,
    positions: &Array2<f64>,
    areas: &[f64],
) -> Result<FeaResult, TheseusError> {
    // 1. Update element geometry
    update_element_geometry(cache, problem, positions);

    // 2. Assemble K
    assemble_global_k(cache, problem, areas)?;

    // 3. Assemble loads
    assemble_loads(cache, problem, areas, positions);

    // 4. Factor and solve
    factor_and_solve(cache)?;

    // 5. Post-process
    post_process(cache, problem, areas, positions)?;

    // 6. Build deformed positions
    let nn = problem.num_nodes;
    let n_total = cache.dof_map.num_total_dofs;
    let mut u_full = vec![0.0; n_total];
    for (fi, &gi) in cache.dof_map.free_dofs.iter().enumerate() {
        u_full[gi] = cache.displacements[fi];
    }

    let mut deformed = positions.clone();
    for i in 0..nn {
        for d in 0..3 {
            deformed[[i, d]] += u_full[i * 3 + d];
        }
    }

    // 7. Validate
    for (e, &f) in cache.axial_forces.iter().enumerate() {
        if !f.is_finite() {
            return Err(TheseusError::Solver(format!(
                "FEA produced non-finite axial force at element {e}.",
            )));
        }
    }

    Ok(FeaResult {
        displacements: u_full,
        deformed_positions: deformed,
        reactions: cache.reactions.clone(),
        axial_forces: cache.axial_forces.clone(),
        stresses: cache.stresses.clone(),
        strains: cache.strains.clone(),
        utilization: cache.utilization.clone(),
    })
}
