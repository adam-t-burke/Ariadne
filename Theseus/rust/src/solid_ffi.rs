//! C-compatible FFI for the solid FEA solver (C# P/Invoke).
//!
//! Same architecture as `fea_ffi.rs`: thin extern "C" functions that marshal
//! raw pointers into safe Rust types and delegate to the core library.

use crate::ffi::set_last_error;
use crate::solid_solve::solve_solid;
use crate::solid_types::*;
use crate::types::TheseusError;
use ndarray::Array2;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::slice;

unsafe fn solid_ffi_guard<F>(f: F) -> i32
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
            set_last_error("internal panic in solid FEA solver (this is a bug — please report it)");
            -2
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Opaque handle
// ─────────────────────────────────────────────────────────────

pub struct SolidHandle {
    pub problem: SolidProblem,
}

// ─────────────────────────────────────────────────────────────
//  Problem construction
// ─────────────────────────────────────────────────────────────

/// Create a new solid FEA problem from raw arrays.
///
/// # Safety
/// All pointers must be valid for the given lengths.
#[no_mangle]
pub unsafe extern "C" fn theseus_solid_create(
    num_nodes: usize,
    num_elements: usize,
    node_positions: *const f64,
    elements: *const i32,
    num_materials: usize,
    materials: *const f64,
    element_props: *const i32,
    num_supports: usize,
    supports: *const i32,
    num_loads: usize,
    load_forces: *const f64,
    load_nodes: *const i32,
    include_self_weight: u8,
    gravity: *const f64,
) -> *mut SolidHandle {
    let result = catch_unwind(AssertUnwindSafe(|| {
        solid_create_inner(
            num_nodes, num_elements, node_positions, elements,
            num_materials, materials, element_props,
            num_supports, supports,
            num_loads, load_forces, load_nodes,
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
            set_last_error("internal panic in theseus_solid_create (this is a bug)");
            std::ptr::null_mut()
        }
    }
}

unsafe fn solid_create_inner(
    num_nodes: usize,
    num_elements: usize,
    node_positions_ptr: *const f64,
    elements_ptr: *const i32,
    num_materials: usize,
    materials_ptr: *const f64,
    element_props_ptr: *const i32,
    num_supports: usize,
    supports_ptr: *const i32,
    num_loads: usize,
    load_forces_ptr: *const f64,
    load_nodes_ptr: *const i32,
    include_self_weight: u8,
    gravity_ptr: *const f64,
) -> Result<*mut SolidHandle, TheseusError> {
    // Node positions: flat [x0,y0,z0, x1,y1,z1, ...]
    let pos_slice = slice::from_raw_parts(node_positions_ptr, num_nodes * 3);
    let node_positions = Array2::from_shape_vec((num_nodes, 3), pos_slice.to_vec())
        .map_err(|e| TheseusError::Shape(format!("node_positions: {e}")))?;

    // Elements: flat [n0,n1,n2,n3, ...] groups of 4 (tet4)
    let elem_slice = slice::from_raw_parts(elements_ptr, num_elements * 4);
    let elements: Vec<[usize; 4]> = (0..num_elements)
        .map(|e| [
            elem_slice[e * 4] as usize,
            elem_slice[e * 4 + 1] as usize,
            elem_slice[e * 4 + 2] as usize,
            elem_slice[e * 4 + 3] as usize,
        ])
        .collect();

    // Materials: flat [E, nu, density, yield_stress, ...] groups of 4
    let mat_slice = slice::from_raw_parts(materials_ptr, num_materials * 4);
    let materials: Vec<SolidMaterial> = (0..num_materials)
        .map(|i| SolidMaterial {
            e: mat_slice[i * 4],
            nu: mat_slice[i * 4 + 1],
            density: mat_slice[i * 4 + 2],
            yield_stress: mat_slice[i * 4 + 3],
        })
        .collect();

    // Element properties: material index per element
    let ep_slice = slice::from_raw_parts(element_props_ptr, num_elements);
    let element_props: Vec<SolidElementProps> = ep_slice
        .iter()
        .map(|&idx| SolidElementProps {
            material_idx: idx as usize,
        })
        .collect();

    // Supports: flat [node_idx, fix_x, fix_y, fix_z, ...] groups of 4
    let sup_slice = slice::from_raw_parts(supports_ptr, num_supports * 4);
    let supports: Vec<SolidSupport> = (0..num_supports)
        .map(|i| SolidSupport {
            node_idx: sup_slice[i * 4] as usize,
            fixed_dofs: [
                sup_slice[i * 4 + 1] != 0,
                sup_slice[i * 4 + 2] != 0,
                sup_slice[i * 4 + 3] != 0,
            ],
        })
        .collect();

    // Loads
    let loads: Vec<SolidLoad> = if num_loads == 0 {
        Vec::new()
    } else {
        let force_slice = slice::from_raw_parts(load_forces_ptr, num_loads * 3);
        let node_slice = slice::from_raw_parts(load_nodes_ptr, num_loads);
        (0..num_loads)
            .map(|i| SolidLoad {
                node_idx: node_slice[i] as usize,
                force: [
                    force_slice[i * 3],
                    force_slice[i * 3 + 1],
                    force_slice[i * 3 + 2],
                ],
            })
            .collect()
    };

    // Gravity
    let grav = slice::from_raw_parts(gravity_ptr, 3);
    let gravity = [grav[0], grav[1], grav[2]];

    let dof_map = SolidDofMap::from_supports(num_nodes, &supports);

    let problem = SolidProblem {
        num_nodes,
        num_elements,
        materials,
        element_props,
        supports,
        loads,
        node_positions,
        elements,
        dof_map,
        include_self_weight: include_self_weight != 0,
        gravity,
    };

    Ok(Box::into_raw(Box::new(SolidHandle { problem })))
}

// ─────────────────────────────────────────────────────────────
//  Forward solve
// ─────────────────────────────────────────────────────────────

/// Run a forward solid FEA solve, writing results into caller-provided buffers.
///
/// Returns 0 on success, -1 on error, -2 on internal panic.
///
/// # Safety
/// All output pointers must have the correct sizes:
/// - `out_displacements`: num_nodes * 3
/// - `out_deformed_xyz`: num_nodes * 3
/// - `out_reactions`: num_nodes * 3
/// - `out_stresses`: num_elements * 4 * 6  (4 GPs × 6 components per element)
/// - `out_von_mises`: num_elements * 4     (4 GPs per element)
#[no_mangle]
pub unsafe extern "C" fn theseus_solid_solve_forward(
    handle: *mut SolidHandle,
    out_displacements: *mut f64,
    out_deformed_xyz: *mut f64,
    out_reactions: *mut f64,
    out_stresses: *mut f64,
    out_von_mises: *mut f64,
) -> i32 {
    solid_ffi_guard(AssertUnwindSafe(|| {
        if handle.is_null() {
            return Err(TheseusError::Solver("null solver handle".into()));
        }
        let h = &mut *handle;
        let mut cache = SolidCache::new(&h.problem)?;
        let result = solve_solid(&mut cache, &h.problem)?;

        let nn = h.problem.num_nodes;
        let ne = h.problem.num_elements;
        let ngp = crate::solid_assembly::NUM_GP;

        // Displacements (flat, num_nodes * 3)
        slice::from_raw_parts_mut(out_displacements, nn * 3)
            .copy_from_slice(&result.displacements);

        // Deformed positions (flat, num_nodes * 3)
        let xyz_out = slice::from_raw_parts_mut(out_deformed_xyz, nn * 3);
        for i in 0..nn {
            for d in 0..3 {
                xyz_out[i * 3 + d] = result.deformed_positions[[i, d]];
            }
        }

        // Reactions (flat, num_nodes * 3)
        slice::from_raw_parts_mut(out_reactions, nn * 3)
            .copy_from_slice(&result.reactions);

        // Stresses: [elem0_gp0(6), elem0_gp1(6), ..., elem0_gp3(6), elem1_gp0(6), ...]
        let stress_out = slice::from_raw_parts_mut(out_stresses, ne * ngp * 6);
        for e in 0..ne {
            for gp in 0..ngp {
                for c in 0..6 {
                    stress_out[e * ngp * 6 + gp * 6 + c] = result.stresses[e][gp][c];
                }
            }
        }

        // Von Mises: [elem0_gp0, elem0_gp1, ..., elem0_gp3, elem1_gp0, ...]
        let vm_out = slice::from_raw_parts_mut(out_von_mises, ne * ngp);
        for e in 0..ne {
            for gp in 0..ngp {
                vm_out[e * ngp + gp] = result.von_mises[e][gp];
            }
        }

        Ok(())
    }))
}

// ─────────────────────────────────────────────────────────────
//  Destruction
// ─────────────────────────────────────────────────────────────

/// Free a solid solver handle.
#[no_mangle]
pub unsafe extern "C" fn theseus_solid_destroy(handle: *mut SolidHandle) {
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

/// Retrieve the last error message (delegates to the shared thread-local in ffi.rs).
///
/// # Safety
/// `buf` must point to at least `buf_len` writable bytes.
#[no_mangle]
pub unsafe extern "C" fn theseus_solid_last_error(buf: *mut u8, buf_len: usize) -> i32 {
    crate::ffi::theseus_last_error(buf, buf_len)
}
