//! Forward FEA solver for tet4 solid elements.
//!
//! Pipeline: assemble K → build RHS → factorise → solve K u = f → post-process.

use crate::solid_assembly::{assemble_global_k, assemble_loads, tet4_b_matrix, constitutive_matrix};
use crate::solid_types::{SolidCache, SolidProblem, SolidResult};
use crate::types::{Factorization, FactorizationStrategy, TheseusError};
use rayon::prelude::*;

// ─────────────────────────────────────────────────────────────
//  Factor and solve
// ─────────────────────────────────────────────────────────────

fn factor_and_solve(cache: &mut SolidCache) -> Result<(), TheseusError> {
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
            "Solid FEA linear solve produced non-finite displacements (singular or ill-conditioned \
             stiffness matrix). Check supports and element connectivity.".into(),
        ));
    }

    cache.displacements.copy_from_slice(&solution);
    Ok(())
}

// ─────────────────────────────────────────────────────────────
//  Post-processing: stress, strain, von Mises, reactions
// ─────────────────────────────────────────────────────────────

fn post_process(
    cache: &mut SolidCache,
    problem: &SolidProblem,
) -> Result<(), TheseusError> {
    let ne = problem.num_elements;
    let n_total = cache.dof_map.num_total_dofs;

    let mut u_full = vec![0.0; n_total];
    for (fi, &gi) in cache.dof_map.free_dofs.iter().enumerate() {
        u_full[gi] = cache.displacements[fi];
    }

    let positions = &problem.node_positions;
    let materials = &problem.materials;
    let element_props = &problem.element_props;
    let elements = &problem.elements;

    // Per-element stress/strain/von Mises (parallel)
    let elem_results: Vec<([f64; 6], [f64; 6], f64)> = (0..ne)
        .into_par_iter()
        .map(|e| {
            let nodes = &elements[e];
            let mat = &materials[element_props[e].material_idx];

            let coords = [
                [positions[[nodes[0], 0]], positions[[nodes[0], 1]], positions[[nodes[0], 2]]],
                [positions[[nodes[1], 0]], positions[[nodes[1], 1]], positions[[nodes[1], 2]]],
                [positions[[nodes[2], 0]], positions[[nodes[2], 1]], positions[[nodes[2], 2]]],
                [positions[[nodes[3], 0]], positions[[nodes[3], 1]], positions[[nodes[3], 2]]],
            ];

            let (b, _vol) = match tet4_b_matrix(&coords) {
                Some(bv) => bv,
                None => return ([0.0; 6], [0.0; 6], 0.0),
            };

            let d = constitutive_matrix(mat.e, mat.nu);

            // u_e: element displacement vector (12 entries)
            let mut u_e = [0.0f64; 12];
            for (i, &ni) in nodes.iter().enumerate() {
                for dd in 0..3 {
                    u_e[i * 3 + dd] = u_full[ni * 3 + dd];
                }
            }

            // strain = B * u_e  (6 entries)
            let mut strain = [0.0f64; 6];
            for i in 0..6 {
                let mut s = 0.0;
                for j in 0..12 {
                    s += b[i][j] * u_e[j];
                }
                strain[i] = s;
            }

            // stress = D * strain  (6 entries)
            let mut stress = [0.0f64; 6];
            for i in 0..6 {
                let mut s = 0.0;
                for j in 0..6 {
                    s += d[i][j] * strain[j];
                }
                stress[i] = s;
            }

            // von Mises
            let sxx = stress[0]; let syy = stress[1]; let szz = stress[2];
            let sxy = stress[3]; let syz = stress[4]; let sxz = stress[5];
            let vm = (0.5 * ((sxx - syy).powi(2) + (syy - szz).powi(2) + (szz - sxx).powi(2)
                + 6.0 * (sxy * sxy + syz * syz + sxz * sxz))).sqrt();

            (stress, strain, vm)
        })
        .collect();

    for (e, (stress, strain, vm)) in elem_results.into_iter().enumerate() {
        cache.stresses[e] = stress;
        cache.strains[e] = strain;
        cache.von_mises[e] = vm;
    }

    // Reactions: R = K_full * u - f  (parallel fold/reduce)
    let n_reactions = cache.reactions.len();

    let reactions = (0..ne)
        .into_par_iter()
        .fold(
            || vec![0.0; n_reactions],
            |mut buf, e| {
                let nodes = &elements[e];
                let mat = &materials[element_props[e].material_idx];

                let coords = [
                    [positions[[nodes[0], 0]], positions[[nodes[0], 1]], positions[[nodes[0], 2]]],
                    [positions[[nodes[1], 0]], positions[[nodes[1], 1]], positions[[nodes[1], 2]]],
                    [positions[[nodes[2], 0]], positions[[nodes[2], 1]], positions[[nodes[2], 2]]],
                    [positions[[nodes[3], 0]], positions[[nodes[3], 1]], positions[[nodes[3], 2]]],
                ];

                let (b, vol) = match tet4_b_matrix(&coords) {
                    Some(bv) => bv,
                    None => return buf,
                };

                let d = constitutive_matrix(mat.e, mat.nu);
                let ke = crate::solid_assembly::tet4_element_stiffness(&b, &d, vol);

                let mut global_dofs = [0usize; 12];
                for (i, &ni) in nodes.iter().enumerate() {
                    for dd in 0..3 {
                        global_dofs[i * 3 + dd] = ni * 3 + dd;
                    }
                }

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

    // Subtract applied loads
    for load in &problem.loads {
        for d in 0..3 {
            cache.reactions[load.node_idx * 3 + d] -= load.force[d];
        }
    }

    if problem.include_self_weight {
        for e in 0..ne {
            let nodes = &elements[e];
            let mat = &materials[element_props[e].material_idx];

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
                        cache.reactions[ni * 3 + d] -= quarter_weight * problem.gravity[d];
                    }
                }
            }
        }
    }

    // Zero out reactions at free DOFs
    for &gi in &cache.dof_map.free_dofs {
        cache.reactions[gi] = 0.0;
    }

    Ok(())
}

// ─────────────────────────────────────────────────────────────
//  Top-level forward solve
// ─────────────────────────────────────────────────────────────

/// Full forward solid FEA solve for tet4 elements. Returns a SolidResult.
pub fn solve_solid(
    cache: &mut SolidCache,
    problem: &SolidProblem,
) -> Result<SolidResult, TheseusError> {
    // 1. Assemble K (parallel)
    assemble_global_k(cache, problem)?;

    // 2. Assemble loads
    assemble_loads(cache, problem);

    // 3. Factor and solve
    factor_and_solve(cache)?;

    // 4. Post-process
    post_process(cache, problem)?;

    // 5. Build full displacement vector and deformed positions
    let nn = problem.num_nodes;
    let n_total = cache.dof_map.num_total_dofs;
    let mut u_full = vec![0.0; n_total];
    for (fi, &gi) in cache.dof_map.free_dofs.iter().enumerate() {
        u_full[gi] = cache.displacements[fi];
    }

    let mut deformed = problem.node_positions.clone();
    for i in 0..nn {
        for d in 0..3 {
            deformed[[i, d]] += u_full[i * 3 + d];
        }
    }

    // 6. Validate
    for (e, &vm) in cache.von_mises.iter().enumerate() {
        if !vm.is_finite() {
            return Err(TheseusError::Solver(format!(
                "Solid FEA produced non-finite von Mises stress at element {e}.",
            )));
        }
    }

    Ok(SolidResult {
        displacements: u_full,
        deformed_positions: deformed,
        reactions: cache.reactions.clone(),
        stresses: cache.stresses.clone(),
        strains: cache.strains.clone(),
        von_mises: cache.von_mises.clone(),
    })
}
