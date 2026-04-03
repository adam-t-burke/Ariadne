//! C-compatible FFI for the shell FEA solver (C# P/Invoke).
//!
//! Exposes 6-DOF shell elements with support for triangles and quadrilaterals.

use crate::ffi::set_last_error;
use crate::shell_solve::solve_shell;
use crate::shell_types::*;
use crate::types::TheseusError;
use ndarray::Array2;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::slice;

unsafe fn shell_ffi_guard<F>(f: F) -> i32
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
            set_last_error("internal panic in shell FEA solver (this is a bug — please report it)");
            -2
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Opaque handle
// ─────────────────────────────────────────────────────────────

pub struct ShellHandle {
    pub problem: ShellProblem,
}

// ─────────────────────────────────────────────────────────────
//  Problem construction
// ─────────────────────────────────────────────────────────────

/// Create a new shell FEA problem from raw arrays.
///
/// # Safety
/// All pointers must be valid for the given lengths.
#[no_mangle]
pub unsafe extern "C" fn theseus_shell_create(
    num_nodes: usize,
    num_elements: usize,
    node_positions: *const f64,
    node_thicknesses: *const f64,
    elements: *const i32,
    num_nodes_per_element: *const i32,
    num_materials: usize,
    materials: *const f64,
    num_sections: usize,
    sections: *const f64,
    element_props: *const i32,
    num_supports: usize,
    supports: *const i32,
    num_loads: usize,
    load_forces: *const f64,
    load_moments: *const f64,
    load_nodes: *const i32,
    include_self_weight: u8,
    gravity: *const f64,
) -> *mut ShellHandle {
    let result = catch_unwind(AssertUnwindSafe(|| {
        shell_create_inner(
            num_nodes, num_elements, node_positions, node_thicknesses, elements, num_nodes_per_element,
            num_materials, materials, num_sections, sections, element_props,
            num_supports, supports,
            num_loads, load_forces, load_moments, load_nodes,
            include_self_weight, gravity,
        )
    }));

    match result {
        Ok(Ok(ptr)) => ptr,
        Ok(Err(e)) => {
            set_last_error(&e.to_string());
            std::ptr::null_mut()
        }
        Err(_panic) => {
            set_last_error("internal panic in theseus_shell_create (this is a bug)");
            std::ptr::null_mut()
        }
    }
}

unsafe fn shell_create_inner(
    num_nodes: usize,
    num_elements: usize,
    node_positions_ptr: *const f64,
    node_thicknesses_ptr: *const f64,
    elements_ptr: *const i32,
    num_nodes_per_element_ptr: *const i32,
    num_materials: usize,
    materials_ptr: *const f64,
    num_sections: usize,
    sections_ptr: *const f64,
    element_props_ptr: *const i32,
    num_supports: usize,
    supports_ptr: *const i32,
    num_loads: usize,
    load_forces_ptr: *const f64,
    load_moments_ptr: *const f64,
    load_nodes_ptr: *const i32,
    include_self_weight: u8,
    gravity_ptr: *const f64,
) -> Result<*mut ShellHandle, TheseusError> {
    // Node positions: [x, y, z, ...]
    let pos_slice = slice::from_raw_parts(node_positions_ptr, num_nodes * 3);
    let node_positions = Array2::from_shape_vec((num_nodes, 3), pos_slice.to_vec())
        .map_err(|e| TheseusError::Shape(format!("node_positions: {e}")))?;

    // Node thicknesses: [t0, t1, t2, ...]
    let node_thicknesses = slice::from_raw_parts(node_thicknesses_ptr, num_nodes).to_vec();

    // Elements: varied number of nodes per element (tri=3, quad=4)
    let num_nodes_per_elem = slice::from_raw_parts(num_nodes_per_element_ptr, num_elements);
    let total_indices: usize = num_nodes_per_elem.iter().map(|&n| n as usize).sum();
    let elem_slice = slice::from_raw_parts(elements_ptr, total_indices);
    
    let mut elements = Vec::with_capacity(num_elements);
    let mut offset = 0;
    for &n in num_nodes_per_elem {
        let n = n as usize;
        let mut nodes = Vec::with_capacity(n);
        for i in 0..n {
            nodes.push(elem_slice[offset + i] as usize);
        }
        elements.push(nodes);
        offset += n;
    }

    // Materials: [E, nu, density, yield_stress, ...]
    let mat_slice = slice::from_raw_parts(materials_ptr, num_materials * 4);
    let materials: Vec<ShellMaterial> = (0..num_materials)
        .map(|i| ShellMaterial {
            e: mat_slice[i * 4],
            nu: mat_slice[i * 4 + 1],
            density: mat_slice[i * 4 + 2],
            yield_stress: mat_slice[i * 4 + 3],
        })
        .collect();

    let sec_slice = if num_sections > 0 {
        slice::from_raw_parts(sections_ptr, num_sections)
    } else {
        &[]
    };
    let sections: Vec<ShellSection> = (0..num_sections)
        .map(|i| ShellSection { offset: sec_slice[i] })
        .collect();

    // Element properties: [material_idx, section_idx, ...]
    let ep_slice = slice::from_raw_parts(element_props_ptr, num_elements * 2);
    let element_props: Vec<ShellElementProps> = (0..num_elements)
        .map(|i| ShellElementProps {
            material_idx: ep_slice[i * 2] as usize,
            section_idx: ep_slice[i * 2 + 1] as usize,
        })
        .collect();

    // Supports: [node_idx, fix_x, fix_y, fix_z, fix_rx, fix_ry, fix_rz, ...]
    let sup_slice = slice::from_raw_parts(supports_ptr, num_supports * 7);
    let supports: Vec<ShellSupport> = (0..num_supports)
        .map(|i| ShellSupport {
            node_idx: sup_slice[i * 7] as usize,
            fixed_dofs: [
                sup_slice[i * 7 + 1] != 0,
                sup_slice[i * 7 + 2] != 0,
                sup_slice[i * 7 + 3] != 0,
                sup_slice[i * 7 + 4] != 0,
                sup_slice[i * 7 + 5] != 0,
                sup_slice[i * 7 + 6] != 0,
            ],
        })
        .collect();

    // Loads: force + moment per loaded node
    let loads: Vec<ShellLoad> = if num_loads == 0 {
        Vec::new()
    } else {
        let force_slice = slice::from_raw_parts(load_forces_ptr, num_loads * 3);
        let moment_slice = slice::from_raw_parts(load_moments_ptr, num_loads * 3);
        let node_slice = slice::from_raw_parts(load_nodes_ptr, num_loads);
        (0..num_loads)
            .map(|i| ShellLoad {
                node_idx: node_slice[i] as usize,
                force: [force_slice[i * 3], force_slice[i * 3 + 1], force_slice[i * 3 + 2]],
                moment: [moment_slice[i * 3], moment_slice[i * 3 + 1], moment_slice[i * 3 + 2]],
            })
            .collect()
    };

    let grav = slice::from_raw_parts(gravity_ptr, 3);
    let gravity = [grav[0], grav[1], grav[2]];

    let dof_map = ShellDofMap::from_supports(num_nodes, &supports);

    let problem = ShellProblem {
        num_nodes,
        num_elements,
        materials,
        sections,
        element_props,
        supports,
        loads,
        node_positions,
        node_thicknesses,
        elements,
        dof_map,
        include_self_weight: include_self_weight != 0,
        gravity,
    };

    Ok(Box::into_raw(Box::new(ShellHandle { problem })))
}

// ─────────────────────────────────────────────────────────────
//  Forward solve
// ─────────────────────────────────────────────────────────────

/// Run a forward shell FEA solve.
///
/// # Safety
/// All output pointers must have the correct sizes:
/// - `out_displacements`: num_nodes * 6
/// - `out_reactions`: num_nodes * 6
#[no_mangle]
pub unsafe extern "C" fn theseus_shell_solve_forward(
    handle: *mut ShellHandle,
    out_displacements: *mut f64,
    out_reactions: *mut f64,
    out_principal_stresses: *mut f64,
    out_top_stresses: *mut f64,
    out_bottom_stresses: *mut f64,
    out_von_mises: *mut f64,
) -> i32 {
    shell_ffi_guard(AssertUnwindSafe(|| {
        if handle.is_null() {
            return Err(TheseusError::Solver("null solver handle".into()));
        }
        let h = &mut *handle;
        let mut cache = ShellCache::new(&h.problem)?;
        let result = solve_shell(&mut cache, &h.problem)?;

        let n_total = h.problem.dof_map.num_total_dofs;
        let ne = h.problem.num_elements;

        slice::from_raw_parts_mut(out_displacements, n_total)
            .copy_from_slice(&result.displacements);

        slice::from_raw_parts_mut(out_reactions, n_total)
            .copy_from_slice(&result.reactions);

        let stress_out = slice::from_raw_parts_mut(out_principal_stresses, ne * 2);
        for e in 0..ne {
            stress_out[e * 2] = result.principal_stresses[e][0];
            stress_out[e * 2 + 1] = result.principal_stresses[e][1];
        }

        let top_out = slice::from_raw_parts_mut(out_top_stresses, ne * 6);
        let bot_out = slice::from_raw_parts_mut(out_bottom_stresses, ne * 6);
        for e in 0..ne {
            for c in 0..6 {
                top_out[e * 6 + c] = result.top_stresses_global[e][c];
                bot_out[e * 6 + c] = result.bottom_stresses_global[e][c];
            }
        }

        let vm_out = slice::from_raw_parts_mut(out_von_mises, ne);
        vm_out.copy_from_slice(&result.von_mises);

        Ok(())
    }))
}

// ─────────────────────────────────────────────────────────────
//  Destruction
// ─────────────────────────────────────────────────────────────

/// Free a shell solver handle.
#[no_mangle]
pub unsafe extern "C" fn theseus_shell_destroy(handle: *mut ShellHandle) {
    if handle.is_null() {
        return;
    }
    let _ = catch_unwind(AssertUnwindSafe(|| {
        drop(Box::from_raw(handle));
    }));
}

// ─────────────────────────────────────────────────────────────
//  Error retrieval
// ─────────────────────────────────────────────────────────────

/// Retrieve the last error message.
///
/// # Safety
/// `buf` must point to at least `buf_len` writable bytes.
#[no_mangle]
pub unsafe extern "C" fn theseus_shell_last_error(buf: *mut u8, buf_len: usize) -> i32 {
    crate::ffi::theseus_last_error(buf, buf_len)
}
