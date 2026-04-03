//! C-compatible FFI for shell SPR recovery and principal stress decomposition.

use crate::ffi::set_last_error;
use crate::shell_spr::spr_shell_full;
use crate::types::TheseusError;
use ndarray::Array2;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::slice;

unsafe fn shell_spr_guard<F>(f: F) -> i32
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
            set_last_error(
                "internal panic in shell SPR recovery (this is a bug — please report it)",
            );
            -2
        }
    }
}

/// Perform shell SPR nodal stress recovery with principal stress decomposition.
///
/// Accepts per-element stress tensors (6-component Voigt, global frame) for
/// membrane, top-fiber, and bottom-fiber fields. Runs SPR on each field,
/// then decomposes recovered nodal tensors into principals + von Mises.
///
/// # Safety
/// All pointers must be valid for the given lengths.
#[no_mangle]
pub unsafe extern "C" fn theseus_shell_spr_recover(
    num_nodes: usize,
    num_elements: usize,
    node_positions: *const f64,
    elements: *const i32,
    num_nodes_per_element: *const i32,
    membrane_stresses: *const f64,
    top_stresses: *const f64,
    bottom_stresses: *const f64,
    // outputs
    out_nodal_membrane: *mut f64,
    out_nodal_top: *mut f64,
    out_nodal_bottom: *mut f64,
    out_element_errors: *mut f64,
    out_principal_values_top: *mut f64,
    out_principal_values_bot: *mut f64,
    out_principal_vectors_top: *mut f64,
    out_principal_vectors_bot: *mut f64,
    out_von_mises_top: *mut f64,
    out_von_mises_bot: *mut f64,
) -> i32 {
    shell_spr_guard(AssertUnwindSafe(|| {
        let pos_slice = slice::from_raw_parts(node_positions, num_nodes * 3);
        let positions = Array2::from_shape_vec((num_nodes, 3), pos_slice.to_vec())
            .map_err(|e| TheseusError::Shape(format!("node_positions: {e}")))?;

        let num_nodes_per_elem = slice::from_raw_parts(num_nodes_per_element, num_elements);
        let total_indices: usize = num_nodes_per_elem.iter().map(|&n| n as usize).sum();
        let elem_slice = slice::from_raw_parts(elements, total_indices);

        let mut elems = Vec::with_capacity(num_elements);
        let mut offset = 0;
        for &n in num_nodes_per_elem {
            let n = n as usize;
            let mut nodes = Vec::with_capacity(n);
            for i in 0..n {
                nodes.push(elem_slice[offset + i] as usize);
            }
            elems.push(nodes);
            offset += n;
        }

        let mem_slice = slice::from_raw_parts(membrane_stresses, num_elements * 6);
        let top_slice = slice::from_raw_parts(top_stresses, num_elements * 6);
        let bot_slice = slice::from_raw_parts(bottom_stresses, num_elements * 6);

        let mem_s: Vec<[f64; 6]> = (0..num_elements)
            .map(|e| {
                let mut s = [0.0; 6];
                s.copy_from_slice(&mem_slice[e * 6..(e + 1) * 6]);
                s
            })
            .collect();
        let top_s: Vec<[f64; 6]> = (0..num_elements)
            .map(|e| {
                let mut s = [0.0; 6];
                s.copy_from_slice(&top_slice[e * 6..(e + 1) * 6]);
                s
            })
            .collect();
        let bot_s: Vec<[f64; 6]> = (0..num_elements)
            .map(|e| {
                let mut s = [0.0; 6];
                s.copy_from_slice(&bot_slice[e * 6..(e + 1) * 6]);
                s
            })
            .collect();

        let result = spr_shell_full(&positions, &elems, &mem_s, &top_s, &bot_s);

        // Copy outputs
        let out_nm = slice::from_raw_parts_mut(out_nodal_membrane, num_nodes * 6);
        let out_nt = slice::from_raw_parts_mut(out_nodal_top, num_nodes * 6);
        let out_nb = slice::from_raw_parts_mut(out_nodal_bottom, num_nodes * 6);
        for (i, s) in result.nodal_membrane.iter().enumerate() {
            for c in 0..6 {
                out_nm[i * 6 + c] = s[c];
            }
        }
        for (i, s) in result.nodal_top.iter().enumerate() {
            for c in 0..6 {
                out_nt[i * 6 + c] = s[c];
            }
        }
        for (i, s) in result.nodal_bottom.iter().enumerate() {
            for c in 0..6 {
                out_nb[i * 6 + c] = s[c];
            }
        }

        let out_err = slice::from_raw_parts_mut(out_element_errors, num_elements);
        out_err.copy_from_slice(&result.element_errors);

        let out_pv_top = slice::from_raw_parts_mut(out_principal_values_top, num_nodes * 3);
        let out_pv_bot = slice::from_raw_parts_mut(out_principal_values_bot, num_nodes * 3);
        let out_pvec_top = slice::from_raw_parts_mut(out_principal_vectors_top, num_nodes * 9);
        let out_pvec_bot = slice::from_raw_parts_mut(out_principal_vectors_bot, num_nodes * 9);
        let out_vm_top = slice::from_raw_parts_mut(out_von_mises_top, num_nodes);
        let out_vm_bot = slice::from_raw_parts_mut(out_von_mises_bot, num_nodes);

        for i in 0..num_nodes {
            for k in 0..3 {
                out_pv_top[i * 3 + k] = result.principal_values_top[i][k];
                out_pv_bot[i * 3 + k] = result.principal_values_bot[i][k];
                for d in 0..3 {
                    out_pvec_top[i * 9 + k * 3 + d] = result.principal_vectors_top[i][k][d];
                    out_pvec_bot[i * 9 + k * 3 + d] = result.principal_vectors_bot[i][k][d];
                }
            }
            out_vm_top[i] = result.von_mises_top[i];
            out_vm_bot[i] = result.von_mises_bot[i];
        }

        Ok(())
    }))
}
