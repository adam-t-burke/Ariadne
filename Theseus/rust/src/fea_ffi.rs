//! C-compatible FFI for the FEA solver (C# P/Invoke).
//!
//! Same architecture as the FDM FFI: thin extern "C" functions that marshal
//! raw pointers into safe Rust types and delegate to the core library.

use crate::fea_objectives::*;
use crate::fea_optimizer::{self, FeaProgressCallback};
use crate::fea_solve::solve_fea;
use crate::fea_types::*;
use crate::ffi::set_last_error;
use crate::types::TheseusError;
use ndarray::Array2;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::slice;

/// Wrap an `extern "C"` body with error handling.
unsafe fn fea_ffi_guard<F>(f: F) -> i32
where
    F: FnOnce() -> Result<(), TheseusError> + std::panic::UnwindSafe,
{
    match catch_unwind(f) {
        Ok(Ok(())) => 0,
        Ok(Err(e)) => {
            set_last_error(&e.to_string());
            -1
        }
        Err(_panic) => {
            set_last_error("internal panic in FEA solver (this is a bug — please report it)");
            -2
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Opaque handle
// ─────────────────────────────────────────────────────────────

pub struct FeaHandle {
    pub problem: FeaProblem,
    pub state: FeaOptimizationState,
    pub mask: FeaVariableMask,
    pub progress_callback: Option<FeaProgressCallback>,
    pub report_frequency: usize,
}

// ─────────────────────────────────────────────────────────────
//  Problem construction
// ─────────────────────────────────────────────────────────────

/// Create a new FEA problem from raw arrays.
///
/// # Safety
/// All pointers must be valid for the given lengths.
#[no_mangle]
pub unsafe extern "C" fn theseus_fea_create(
    num_nodes: usize,
    num_elements: usize,
    node_positions: *const f64,
    edge_nodes: *const usize,
    num_materials: usize,
    materials: *const f64,
    num_sections: usize,
    sections: *const f64,
    element_props: *const usize,
    num_supports: usize,
    supports: *const usize,
    num_loads: usize,
    loads_data: *const f64,
    load_moments: *const f64,
    load_nodes: *const usize,
    include_self_weight: i32,
    gravity: *const f64,
    beam_formulation: u8,
) -> *mut FeaHandle {
    let result = catch_unwind(AssertUnwindSafe(|| {
        fea_create_inner(
            num_nodes, num_elements, node_positions, edge_nodes,
            num_materials, materials, num_sections, sections, element_props,
            num_supports, supports, num_loads, loads_data, load_moments, load_nodes,
            include_self_weight, gravity, beam_formulation,
        )
    }));

    match result {
        Ok(Ok(ptr)) => ptr,
        Ok(Err(e)) => {
            set_last_error(&e.to_string());
            std::ptr::null_mut()
        }
        Err(_panic) => {
            set_last_error("internal panic in theseus_fea_create (this is a bug)");
            std::ptr::null_mut()
        }
    }
}

unsafe fn fea_create_inner(
    num_nodes: usize,
    num_elements: usize,
    node_positions_ptr: *const f64,
    edge_nodes_ptr: *const usize,
    num_materials: usize,
    materials_ptr: *const f64,
    num_sections: usize,
    sections_ptr: *const f64,
    element_props_ptr: *const usize,
    num_supports: usize,
    supports_ptr: *const usize,
    num_loads: usize,
    loads_data_ptr: *const f64,
    load_moments_ptr: *const f64,
    load_nodes_ptr: *const usize,
    include_self_weight: i32,
    gravity_ptr: *const f64,
    beam_formulation_raw: u8,
) -> Result<*mut FeaHandle, TheseusError> {
    let beam_formulation = match beam_formulation_raw {
        1 => BeamFormulation::EulerBernoulli,
        2 => BeamFormulation::Timoshenko,
        _ => BeamFormulation::Truss,
    };
    let dofs_per_node: usize = if beam_formulation == BeamFormulation::Truss { 3 } else { 6 };

    let pos_slice = slice::from_raw_parts(node_positions_ptr, num_nodes * 3);
    let node_positions = Array2::from_shape_vec((num_nodes, 3), pos_slice.to_vec())
        .map_err(|e| TheseusError::Shape(format!("node_positions: {e}")))?;

    let edge_slice = slice::from_raw_parts(edge_nodes_ptr, num_elements * 2);
    let edge_nodes: Vec<(usize, usize)> = (0..num_elements)
        .map(|e| (edge_slice[e * 2], edge_slice[e * 2 + 1]))
        .collect();

    let mat_slice = slice::from_raw_parts(materials_ptr, num_materials * 3);
    let materials: Vec<FeaMaterial> = (0..num_materials)
        .map(|i| FeaMaterial {
            e: mat_slice[i * 3],
            density: mat_slice[i * 3 + 1],
            yield_stress: mat_slice[i * 3 + 2],
        })
        .collect();

    // Sections: 6 doubles per section [A, Iy, Iz, J, Asy, Asz]
    let sec_slice = slice::from_raw_parts(sections_ptr, num_sections * 6);
    let sections: Vec<FeaSection> = (0..num_sections)
        .map(|i| FeaSection {
            area: sec_slice[i * 6],
            iy: sec_slice[i * 6 + 1],
            iz: sec_slice[i * 6 + 2],
            j: sec_slice[i * 6 + 3],
            asy: sec_slice[i * 6 + 4],
            asz: sec_slice[i * 6 + 5],
        })
        .collect();

    let ep_slice = slice::from_raw_parts(element_props_ptr, num_elements * 2);
    let element_props: Vec<FeaElementProps> = (0..num_elements)
        .map(|e| FeaElementProps {
            material_idx: ep_slice[e * 2],
            section_idx: ep_slice[e * 2 + 1],
        })
        .collect();

    // Supports: 7 usize per support [node, fix_x, fix_y, fix_z, fix_rx, fix_ry, fix_rz]
    let sup_stride = 7;
    let sup_slice = slice::from_raw_parts(supports_ptr, num_supports * sup_stride);
    let supports: Vec<FeaSupport> = (0..num_supports)
        .map(|i| {
            let base = i * sup_stride;
            let mut fixed = Vec::with_capacity(dofs_per_node);
            for d in 0..dofs_per_node {
                fixed.push(sup_slice[base + 1 + d] != 0);
            }
            FeaSupport {
                node_idx: sup_slice[base],
                fixed_dofs: fixed,
            }
        })
        .collect();

    // Loads: forces (3) + moments (3) per load
    let load_forces = slice::from_raw_parts(loads_data_ptr, num_loads * 3);
    let load_moms = slice::from_raw_parts(load_moments_ptr, num_loads * 3);
    let load_node_indices = slice::from_raw_parts(load_nodes_ptr, num_loads);
    let loads: Vec<FeaLoad> = (0..num_loads)
        .map(|i| FeaLoad {
            node_idx: load_node_indices[i],
            force: [load_forces[i * 3], load_forces[i * 3 + 1], load_forces[i * 3 + 2]],
            moment: [load_moms[i * 3], load_moms[i * 3 + 1], load_moms[i * 3 + 2]],
        })
        .collect();

    let grav = slice::from_raw_parts(gravity_ptr, 3);
    let gravity = [grav[0], grav[1], grav[2]];

    let dof_map = DofMap::from_supports(num_nodes, dofs_per_node, &supports);

    let areas: Vec<f64> = sections.iter().map(|s| s.area).collect();

    let problem = FeaProblem {
        num_nodes,
        num_elements,
        beam_formulation,
        dofs_per_node,
        materials,
        sections,
        element_props,
        supports,
        loads,
        node_positions: node_positions.clone(),
        edge_nodes,
        dof_map,
        include_self_weight: include_self_weight != 0,
        gravity,
        objectives: Vec::new(),
        solver: FeaSolverOptions::default(),
    };

    let state = FeaOptimizationState::new(node_positions, areas);

    Ok(Box::into_raw(Box::new(FeaHandle {
        problem,
        state,
        mask: FeaVariableMask::default(),
        progress_callback: None,
        report_frequency: 1,
    })))
}

/// Free a handle.
#[no_mangle]
pub unsafe extern "C" fn theseus_fea_free(handle: *mut FeaHandle) {
    if handle.is_null() { return; }
    let _ = catch_unwind(AssertUnwindSafe(|| {
        drop(Box::from_raw(handle));
    }));
}

// ─────────────────────────────────────────────────────────────
//  Objective registration
// ─────────────────────────────────────────────────────────────

#[no_mangle]
pub unsafe extern "C" fn theseus_fea_add_compliance(
    handle: *mut FeaHandle,
    weight: f64,
) -> i32 {
    fea_ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        h.problem.objectives.push(Box::new(FeaCompliance { weight }));
        Ok(())
    }))
}

#[no_mangle]
pub unsafe extern "C" fn theseus_fea_add_max_displacement(
    handle: *mut FeaHandle,
    weight: f64,
    node_indices: *const usize,
    num_nodes: usize,
) -> i32 {
    fea_ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let idx = slice::from_raw_parts(node_indices, num_nodes).to_vec();
        h.problem.objectives.push(Box::new(FeaMaxDisplacement { weight, node_indices: idx }));
        Ok(())
    }))
}

#[no_mangle]
pub unsafe extern "C" fn theseus_fea_add_target_displacement(
    handle: *mut FeaHandle,
    weight: f64,
    node_indices: *const usize,
    num_nodes: usize,
    targets: *const f64,
) -> i32 {
    fea_ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let idx = slice::from_raw_parts(node_indices, num_nodes).to_vec();
        let tgt_slice = slice::from_raw_parts(targets, num_nodes * 3);
        let tgts: Vec<[f64; 3]> = (0..num_nodes)
            .map(|i| [tgt_slice[i * 3], tgt_slice[i * 3 + 1], tgt_slice[i * 3 + 2]])
            .collect();
        h.problem.objectives.push(Box::new(FeaTargetDisplacement {
            weight, node_indices: idx, targets: tgts,
        }));
        Ok(())
    }))
}

#[no_mangle]
pub unsafe extern "C" fn theseus_fea_add_min_weight(
    handle: *mut FeaHandle,
    weight: f64,
) -> i32 {
    fea_ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        h.problem.objectives.push(Box::new(FeaMinWeight { weight }));
        Ok(())
    }))
}

#[no_mangle]
pub unsafe extern "C" fn theseus_fea_add_max_stress(
    handle: *mut FeaHandle,
    weight: f64,
    edge_indices: *const usize,
    num_edges: usize,
    thresholds: *const f64,
    sharpness: f64,
) -> i32 {
    fea_ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let idx = slice::from_raw_parts(edge_indices, num_edges).to_vec();
        let thr = slice::from_raw_parts(thresholds, num_edges).to_vec();
        h.problem.objectives.push(Box::new(FeaMaxStress {
            weight, edge_indices: idx, threshold: thr, sharpness,
        }));
        Ok(())
    }))
}

#[no_mangle]
pub unsafe extern "C" fn theseus_fea_add_target_geometry(
    handle: *mut FeaHandle,
    weight: f64,
    node_indices: *const usize,
    num_nodes: usize,
    targets: *const f64,
) -> i32 {
    fea_ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let idx = slice::from_raw_parts(node_indices, num_nodes).to_vec();
        let tgt_slice = slice::from_raw_parts(targets, num_nodes * 3);
        let tgts: Vec<[f64; 3]> = (0..num_nodes)
            .map(|i| [tgt_slice[i * 3], tgt_slice[i * 3 + 1], tgt_slice[i * 3 + 2]])
            .collect();
        h.problem.objectives.push(Box::new(FeaTargetGeometry {
            weight, node_indices: idx, targets: tgts,
        }));
        Ok(())
    }))
}

// ─────────────────────────────────────────────────────────────
//  Solver options
// ─────────────────────────────────────────────────────────────

#[no_mangle]
pub unsafe extern "C" fn theseus_fea_set_solver_options(
    handle: *mut FeaHandle,
    max_iterations: usize,
    abs_tol: f64,
    rel_tol: f64,
    barrier_weight: f64,
    barrier_sharpness: f64,
) -> i32 {
    fea_ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        h.problem.solver = FeaSolverOptions {
            max_iterations,
            absolute_tolerance: abs_tol,
            relative_tolerance: rel_tol,
            report_frequency: 1,
            barrier_weight,
            barrier_sharpness,
        };
        Ok(())
    }))
}

// ─────────────────────────────────────────────────────────────
//  Variable mask
// ─────────────────────────────────────────────────────────────

#[no_mangle]
pub unsafe extern "C" fn theseus_fea_set_variable_mask(
    handle: *mut FeaHandle,
    opt_node_positions: i32,
    node_indices: *const usize,
    num_node_indices: usize,
    opt_areas: i32,
    area_indices: *const usize,
    num_area_indices: usize,
    opt_support_positions: i32,
    support_node_indices: *const usize,
    num_support_indices: usize,
) -> i32 {
    fea_ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        h.mask = FeaVariableMask {
            node_positions: opt_node_positions != 0,
            node_indices: slice::from_raw_parts(node_indices, num_node_indices).to_vec(),
            cross_section_areas: opt_areas != 0,
            area_element_indices: slice::from_raw_parts(area_indices, num_area_indices).to_vec(),
            support_positions: opt_support_positions != 0,
            support_node_indices: slice::from_raw_parts(support_node_indices, num_support_indices).to_vec(),
        };
        Ok(())
    }))
}

// ─────────────────────────────────────────────────────────────
//  Progress callback
// ─────────────────────────────────────────────────────────────

#[no_mangle]
pub unsafe extern "C" fn theseus_fea_set_progress_callback(
    handle: *mut FeaHandle,
    callback: Option<FeaProgressCallback>,
    frequency: usize,
) -> i32 {
    fea_ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        h.progress_callback = callback;
        h.report_frequency = if frequency == 0 { 1 } else { frequency };
        Ok(())
    }))
}

// ─────────────────────────────────────────────────────────────
//  Forward solve
// ─────────────────────────────────────────────────────────────

#[no_mangle]
pub unsafe extern "C" fn theseus_fea_solve_forward(
    handle: *mut FeaHandle,
    out_displacements: *mut f64,
    out_reactions: *mut f64,
    out_axial_forces: *mut f64,
    out_stresses: *mut f64,
    out_strains: *mut f64,
    out_deformed_xyz: *mut f64,
    out_utilization: *mut f64,
    out_internal_forces: *mut f64,
) -> i32 {
    fea_ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let areas: Vec<f64> = h.problem.sections.iter().map(|s| s.area).collect();
        let mut cache = FeaCache::new(&h.problem)?;
        let result = solve_fea(&mut cache, &h.problem, &h.problem.node_positions, &areas)?;

        let nn = h.problem.num_nodes;
        let ne = h.problem.num_elements;
        let dpn = h.problem.dofs_per_node;
        let n_total = nn * dpn;

        slice::from_raw_parts_mut(out_displacements, n_total).copy_from_slice(&result.displacements);
        slice::from_raw_parts_mut(out_reactions, n_total).copy_from_slice(&result.reactions);
        slice::from_raw_parts_mut(out_axial_forces, ne).copy_from_slice(&result.axial_forces);
        slice::from_raw_parts_mut(out_stresses, ne).copy_from_slice(&result.stresses);
        slice::from_raw_parts_mut(out_strains, ne).copy_from_slice(&result.strains);
        slice::from_raw_parts_mut(out_utilization, ne).copy_from_slice(&result.utilization);

        let xyz_out = slice::from_raw_parts_mut(out_deformed_xyz, nn * 3);
        for i in 0..nn {
            for d in 0..3 {
                xyz_out[i * 3 + d] = result.deformed_positions[[i, d]];
            }
        }

        if dpn == 6 && !result.internal_forces.is_empty() {
            let if_out = slice::from_raw_parts_mut(out_internal_forces, ne * 12);
            if_out.copy_from_slice(&result.internal_forces);
        }

        Ok(())
    }))
}

// ─────────────────────────────────────────────────────────────
//  Optimization
// ─────────────────────────────────────────────────────────────

#[no_mangle]
pub unsafe extern "C" fn theseus_fea_optimize(
    handle: *mut FeaHandle,
    out_displacements: *mut f64,
    out_reactions: *mut f64,
    out_axial_forces: *mut f64,
    out_stresses: *mut f64,
    out_strains: *mut f64,
    out_deformed_xyz: *mut f64,
    out_utilization: *mut f64,
    out_node_positions: *mut f64,
    out_areas: *mut f64,
    out_iterations: *mut usize,
    out_converged: *mut bool,
) -> i32 {
    fea_ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let cb = h.progress_callback;
        let freq = h.report_frequency;
        let mask = h.mask.clone();
        let result = fea_optimizer::fea_optimize(
            &h.problem, &mut h.state, &mask, cb, freq,
        )?;

        let nn = h.problem.num_nodes;
        let ne = h.problem.num_elements;
        let dpn = h.problem.dofs_per_node;
        let n_total = nn * dpn;
        let ns = h.problem.sections.len();

        let r = &result.fea_result;
        slice::from_raw_parts_mut(out_displacements, n_total).copy_from_slice(&r.displacements);
        slice::from_raw_parts_mut(out_reactions, n_total).copy_from_slice(&r.reactions);
        slice::from_raw_parts_mut(out_axial_forces, ne).copy_from_slice(&r.axial_forces);
        slice::from_raw_parts_mut(out_stresses, ne).copy_from_slice(&r.stresses);
        slice::from_raw_parts_mut(out_strains, ne).copy_from_slice(&r.strains);
        slice::from_raw_parts_mut(out_utilization, ne).copy_from_slice(&r.utilization);

        let xyz_out = slice::from_raw_parts_mut(out_deformed_xyz, nn * 3);
        for i in 0..nn {
            for d in 0..3 {
                xyz_out[i * 3 + d] = r.deformed_positions[[i, d]];
            }
        }

        let pos_out = slice::from_raw_parts_mut(out_node_positions, nn * 3);
        for i in 0..nn {
            for d in 0..3 {
                pos_out[i * 3 + d] = result.node_positions[[i, d]];
            }
        }

        slice::from_raw_parts_mut(out_areas, ns).copy_from_slice(&result.areas);
        *out_iterations = result.iterations;
        *out_converged = result.converged;

        Ok(())
    }))
}
