//! Forward FEA solver for tet4 solid elements with 4-point quadrature.
//!
//! Pipeline: assemble K → build RHS → factorise → solve K u = f → post-process.

use crate::solid_assembly::{
    assemble_global_k, assemble_loads, precompute_d_matrices,
    TET4_GAUSS_POINTS, NUM_GP,
};
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
//  Post-processing: merged stress/strain/VM + reactions
// ─────────────────────────────────────────────────────────────

/// Per-element result produced by the merged parallel pass.
struct ElemPostResult {
    stresses: [[f64; 6]; NUM_GP],
    strains: [[f64; 6]; NUM_GP],
    von_mises: [f64; NUM_GP],
    reaction_contribs: [f64; 30],
    global_dofs: [usize; 30],
    num_dofs: usize,
}

fn post_process(
    cache: &mut SolidCache,
    problem: &SolidProblem,
    u_full: &[f64],
) -> Result<(), TheseusError> {
    let ne = problem.num_elements;

    let positions = &problem.node_positions;
    let element_props = &problem.element_props;
    let elements = &problem.elements;

    let d_matrices = precompute_d_matrices(problem);

    let elem_results: Vec<ElemPostResult> = (0..ne)
        .into_par_iter()
        .map(|e| {
            let nodes = &elements[e];
            let d = &d_matrices[element_props[e].material_idx];
            let num_nodes = nodes.len();
            let num_dofs = num_nodes * 3;

            let mut u_e = [0.0f64; 30];
            let mut global_dofs = [0usize; 30];
            for (i, &ni) in nodes.iter().enumerate() {
                for dd in 0..3 {
                    u_e[i * 3 + dd] = u_full[ni * 3 + dd];
                    global_dofs[i * 3 + dd] = ni * 3 + dd;
                }
            }

            let mut gp_stress = [[0.0f64; 6]; NUM_GP];
            let mut gp_strain = [[0.0f64; 6]; NUM_GP];
            let mut gp_vm = [0.0f64; NUM_GP];
            let mut r_local = [0.0f64; 30];

            if num_nodes == 10 {
                let coords = crate::solid_assembly::extract_tet10_coords(positions, nodes);
                for (gp, &(pt, _w)) in TET4_GAUSS_POINTS.iter().enumerate() {
                    let (b, _) = match crate::solid_assembly::tet10_b_matrix_at(&coords, pt[0], pt[1], pt[2]) {
                        Some(bv) => bv,
                        None => {
                            return ElemPostResult {
                                stresses: [[0.0; 6]; NUM_GP],
                                strains: [[0.0; 6]; NUM_GP],
                                von_mises: [0.0; NUM_GP],
                                reaction_contribs: [0.0; 30],
                                global_dofs: [0; 30],
                                num_dofs,
                            };
                        }
                    };

                    let mut strain = [0.0f64; 6];
                    for i in 0..6 {
                        let mut s = 0.0;
                        for j in 0..30 {
                            s += b[i][j] * u_e[j];
                        }
                        strain[i] = s;
                    }

                    let mut stress = [0.0f64; 6];
                    for i in 0..6 {
                        let mut s = 0.0;
                        for j in 0..6 {
                            s += d[i][j] * strain[j];
                        }
                        stress[i] = s;
                    }

                    let sxx = stress[0]; let syy = stress[1]; let szz = stress[2];
                    let sxy = stress[3]; let syz = stress[4]; let sxz = stress[5];
                    let vm = (0.5 * ((sxx - syy).powi(2) + (syy - szz).powi(2) + (szz - sxx).powi(2)
                        + 6.0 * (sxy * sxy + syz * syz + sxz * sxz))).sqrt();

                    gp_stress[gp] = stress;
                    gp_strain[gp] = strain;
                    gp_vm[gp] = vm;
                }

                let ke = crate::solid_assembly::tet10_element_stiffness_quad(&coords, d)
                    .unwrap_or([[0.0; 30]; 30]);

                for li in 0..30 {
                    let mut s = 0.0;
                    for lj in 0..30 {
                        s += ke[li][lj] * u_e[lj];
                    }
                    r_local[li] = s;
                }
            } else {
                let coords = crate::solid_assembly::extract_tet4_coords(positions, nodes);
                for (gp, &(pt, _w)) in TET4_GAUSS_POINTS.iter().enumerate() {
                    let (b, _) = match crate::solid_assembly::tet4_b_matrix_at(&coords, pt[0], pt[1], pt[2]) {
                        Some(bv) => bv,
                        None => {
                            return ElemPostResult {
                                stresses: [[0.0; 6]; NUM_GP],
                                strains: [[0.0; 6]; NUM_GP],
                                von_mises: [0.0; NUM_GP],
                                reaction_contribs: [0.0; 30],
                                global_dofs: [0; 30],
                                num_dofs,
                            };
                        }
                    };

                    let mut strain = [0.0f64; 6];
                    for i in 0..6 {
                        let mut s = 0.0;
                        for j in 0..12 {
                            s += b[i][j] * u_e[j];
                        }
                        strain[i] = s;
                    }

                    let mut stress = [0.0f64; 6];
                    for i in 0..6 {
                        let mut s = 0.0;
                        for j in 0..6 {
                            s += d[i][j] * strain[j];
                        }
                        stress[i] = s;
                    }

                    let sxx = stress[0]; let syy = stress[1]; let szz = stress[2];
                    let sxy = stress[3]; let syz = stress[4]; let sxz = stress[5];
                    let vm = (0.5 * ((sxx - syy).powi(2) + (syy - szz).powi(2) + (szz - sxx).powi(2)
                        + 6.0 * (sxy * sxy + syz * syz + sxz * sxz))).sqrt();

                    gp_stress[gp] = stress;
                    gp_strain[gp] = strain;
                    gp_vm[gp] = vm;
                }

                let ke = crate::solid_assembly::tet4_element_stiffness_quad(&coords, d)
                    .unwrap_or([[0.0; 12]; 12]);

                for li in 0..12 {
                    let mut s = 0.0;
                    for lj in 0..12 {
                        s += ke[li][lj] * u_e[lj];
                    }
                    r_local[li] = s;
                }
            }

            ElemPostResult {
                stresses: gp_stress,
                strains: gp_strain,
                von_mises: gp_vm,
                reaction_contribs: r_local,
                global_dofs,
                num_dofs,
            }
        })
        .collect();

    // Scatter results into cache
    let n_reactions = cache.reactions.len();
    for v in cache.reactions.iter_mut() { *v = 0.0; }

    for (e, res) in elem_results.into_iter().enumerate() {
        cache.stresses[e] = res.stresses;
        cache.strains[e] = res.strains;
        cache.von_mises[e] = res.von_mises;

        for li in 0..res.num_dofs {
            let gi = res.global_dofs[li];
            if gi < n_reactions {
                cache.reactions[gi] += res.reaction_contribs[li];
            }
        }
    }

    // Subtract applied loads
    for load in &problem.loads {
        for d in 0..3 {
            cache.reactions[load.node_idx * 3 + d] -= load.force[d];
        }
    }

    if problem.include_self_weight {
        for e in 0..ne {
            let nodes = &problem.elements[e];
            let mat = &problem.materials[problem.element_props[e].material_idx];
            
            if nodes.len() == 10 {
                let coords = crate::solid_assembly::extract_tet10_coords(&problem.node_positions, nodes);
                for &(pt, w) in &TET4_GAUSS_POINTS {
                    if let Some((_, abs_det_j)) = crate::solid_assembly::tet10_b_matrix_at(&coords, pt[0], pt[1], pt[2]) {
                        let n_vals = crate::solid_assembly::tet10_shape_functions(pt[0], pt[1], pt[2]);
                        let factor = w * abs_det_j * mat.density;
                        for node_local in 0..10 {
                            let ni = nodes[node_local];
                            for d in 0..3 {
                                cache.reactions[ni * 3 + d] -= factor * n_vals[node_local] * problem.gravity[d];
                            }
                        }
                    }
                }
            } else {
                let coords = crate::solid_assembly::extract_tet4_coords(&problem.node_positions, nodes);
                for &(pt, w) in &TET4_GAUSS_POINTS {
                    if let Some((_, abs_det_j)) = crate::solid_assembly::tet4_b_matrix_at(&coords, pt[0], pt[1], pt[2]) {
                        let n_vals = crate::solid_assembly::tet4_shape_functions(pt[0], pt[1], pt[2]);
                        let factor = w * abs_det_j * mat.density;
                        for node_local in 0..4 {
                            let ni = nodes[node_local];
                            for d in 0..3 {
                                cache.reactions[ni * 3 + d] -= factor * n_vals[node_local] * problem.gravity[d];
                            }
                        }
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

    // 4. Build full displacement vector (once, shared with post_process)
    let nn = problem.num_nodes;
    let n_total = cache.dof_map.num_total_dofs;
    let mut u_full = vec![0.0; n_total];
    for (fi, &gi) in cache.dof_map.free_dofs.iter().enumerate() {
        u_full[gi] = cache.displacements[fi];
    }

    // 5. Post-process (stress, strain, VM, reactions -- single merged pass)
    post_process(cache, problem, &u_full)?;

    // 6. Deformed positions
    let mut deformed = problem.node_positions.clone();
    for i in 0..nn {
        for d in 0..3 {
            deformed[[i, d]] += u_full[i * 3 + d];
        }
    }

    // 7. Validate
    for (e, gp_vm) in cache.von_mises.iter().enumerate() {
        for (gp, &vm) in gp_vm.iter().enumerate() {
            if !vm.is_finite() {
                return Err(TheseusError::Solver(format!(
                    "Solid FEA produced non-finite von Mises stress at element {e}, GP {gp}.",
                )));
            }
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
