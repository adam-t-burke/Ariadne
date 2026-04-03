//! Forward FEA solver: assemble K, build RHS, factorise, solve K u = f.
//!
//! Post-processing: axial forces, stresses, strains, reactions, utilization.
//! For beams: internal forces (N, Vy, Vz, T, My, Mz) per element end.

use crate::fea_assembly::{assemble_global_k, assemble_loads, bar_element_stiffness, beam_element_stiffness, update_element_geometry};
use crate::fea_types::{BeamFormulation, FeaCache, FeaProblem, FeaResult};
use crate::types::{Factorization, FactorizationStrategy, TheseusError};
use ndarray::Array2;
use rayon::prelude::*;

// ─────────────────────────────────────────────────────────────
//  Factor and solve
// ─────────────────────────────────────────────────────────────

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

fn post_process(
    cache: &mut FeaCache,
    problem: &FeaProblem,
    areas: &[f64],
    _positions: &Array2<f64>,
) -> Result<(), TheseusError> {
    let ne = problem.num_elements;
    let n_total = cache.dof_map.num_total_dofs;
    let dpn = problem.dofs_per_node;
    let is_beam = problem.beam_formulation != BeamFormulation::Truss;

    let mut u_full = vec![0.0; n_total];
    for (fi, &gi) in cache.dof_map.free_dofs.iter().enumerate() {
        u_full[gi] = cache.displacements[fi];
    }

    // Per-element truss results (axial force, stress, strain, utilization)
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
                *af = 0.0; *st = 0.0; *sn = 0.0; *ut = 0.0;
                return;
            }

            let mut delta = 0.0;
            for d in 0..3 {
                delta += c[d] * (u_full[nj * dpn + d] - u_full[ni * dpn + d]);
            }

            let strain = delta / l;
            let stress = mat.e * strain;
            let force = stress * a;

            *sn = strain;
            *st = stress;
            *af = force;
            *ut = if mat.yield_stress > 0.0 { stress.abs() / mat.yield_stress } else { 0.0 };
        });

    // Beam internal forces (if beam mode)
    if is_beam && !cache.internal_forces.is_empty() {
        for e in 0..ne {
            let (ni, nj) = problem.edge_nodes[e];
            let props = &problem.element_props[e];
            let mat = &problem.materials[props.material_idx];
            let sec = &problem.sections[props.section_idx];
            let l = cache.elem_lengths[e];
            let c = &cache.elem_cos[e];

            if l < f64::EPSILON { continue; }

            let kg = beam_element_stiffness(mat.e, 0.3, sec, l, c, problem.beam_formulation);

            // Extract element DOFs from global displacement
            let mut u_e = [0.0f64; 12];
            for d in 0..6 { u_e[d] = u_full[ni * 6 + d]; }
            for d in 0..6 { u_e[6 + d] = u_full[nj * 6 + d]; }

            // f = K * u (in global coords, which gives element end forces)
            for i in 0..12 {
                let mut s = 0.0;
                for j in 0..12 {
                    s += kg[i * 12 + j] * u_e[j];
                }
                cache.internal_forces[e * 12 + i] = s;
            }
        }
    }

    // Reactions: R = K_full * u - f
    let n_reactions = cache.reactions.len();
    let elem_dofs = dpn * 2;

    let reactions = (0..ne)
        .into_par_iter()
        .fold(
            || vec![0.0; n_reactions],
            |mut buf, e| {
                let (ni, nj) = problem.edge_nodes[e];
                let props = &problem.element_props[e];
                let mat = &problem.materials[props.material_idx];
                let l = cache.elem_lengths[e];
                if l < f64::EPSILON { return buf; }
                let c = &cache.elem_cos[e];

                if is_beam {
                    let sec = &problem.sections[props.section_idx];
                    let kg = beam_element_stiffness(mat.e, 0.3, sec, l, c, problem.beam_formulation);

                    let mut global_dofs = Vec::with_capacity(elem_dofs);
                    for d in 0..dpn { global_dofs.push(ni * dpn + d); }
                    for d in 0..dpn { global_dofs.push(nj * dpn + d); }

                    for (li, &gi) in global_dofs.iter().enumerate() {
                        let mut ke_u = 0.0;
                        for (lj, &gj) in global_dofs.iter().enumerate() {
                            ke_u += kg[li * 12 + lj] * u_full[gj];
                        }
                        buf[gi] += ke_u;
                    }
                } else {
                    let a = areas[props.section_idx];
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

    for load in &problem.loads {
        for d in 0..3 {
            cache.reactions[load.node_idx * dpn + d] -= load.force[d];
        }
        if dpn == 6 {
            for d in 0..3 {
                cache.reactions[load.node_idx * dpn + 3 + d] -= load.moment[d];
            }
        }
    }
    if problem.include_self_weight {
        for e in 0..ne {
            let props = &problem.element_props[e];
            let mat = &problem.materials[props.material_idx];
            let a = areas[props.section_idx];
            let l = cache.elem_lengths[e];
            let half_weight = mat.density * a * l * 0.5;
            let (ni, nj) = problem.edge_nodes[e];
            for d in 0..3 {
                let f = half_weight * problem.gravity[d];
                cache.reactions[ni * dpn + d] -= f;
                cache.reactions[nj * dpn + d] -= f;
            }
        }
    }

    for &gi in &cache.dof_map.free_dofs {
        cache.reactions[gi] = 0.0;
    }

    Ok(())
}

// ─────────────────────────────────────────────────────────────
//  Top-level forward solve
// ─────────────────────────────────────────────────────────────

pub fn solve_fea(
    cache: &mut FeaCache,
    problem: &FeaProblem,
    positions: &Array2<f64>,
    areas: &[f64],
) -> Result<FeaResult, TheseusError> {
    update_element_geometry(cache, problem, positions);
    assemble_global_k(cache, problem, areas)?;
    assemble_loads(cache, problem, areas, positions);
    factor_and_solve(cache)?;
    post_process(cache, problem, areas, positions)?;

    let nn = problem.num_nodes;
    let n_total = cache.dof_map.num_total_dofs;
    let dpn = problem.dofs_per_node;

    let mut u_full = vec![0.0; n_total];
    for (fi, &gi) in cache.dof_map.free_dofs.iter().enumerate() {
        u_full[gi] = cache.displacements[fi];
    }

    let mut deformed = positions.clone();
    for i in 0..nn {
        for d in 0..3 {
            deformed[[i, d]] += u_full[i * dpn + d];
        }
    }

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
        internal_forces: cache.internal_forces.clone(),
    })
}
