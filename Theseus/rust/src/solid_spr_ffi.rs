//! C-compatible FFI for SPR recovery and principal stress decomposition.

use crate::ffi::set_last_error;
use crate::solid_spr::{principal_stresses, spr_recover_nodal_stresses, von_mises_from_principals};
use crate::types::TheseusError;
use ndarray::Array2;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::slice;

unsafe fn spr_ffi_guard<F>(f: F) -> i32
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
            set_last_error("internal panic in SPR recovery (this is a bug — please report it)");
            -2
        }
    }
}

/// Perform SPR nodal stress recovery and principal stress decomposition.
///
/// # Safety
/// All pointers must be valid for the given lengths:
/// - `node_positions`: `num_nodes * 3` f64
/// - `elements`: `num_elements * 4` i32
/// - `element_stresses`: `num_elements * 6` f64
/// - `nodal_stresses`: `num_nodes * 6` f64 (output)
/// - `element_errors`: `num_elements` f64 (output)
/// - `principal_values`: `num_nodes * 3` f64 (output)
/// - `principal_vectors`: `num_nodes * 9` f64 (output)
/// - `nodal_von_mises`: `num_nodes` f64 (output)
///
/// Returns 0 on success, -1 on error, -2 on internal panic.
#[no_mangle]
pub unsafe extern "C" fn theseus_solid_spr_recover(
    num_nodes: usize,
    num_elements: usize,
    node_positions: *const f64,
    elements: *const i32,
    element_stresses: *const f64,
    nodal_stresses: *mut f64,
    element_errors: *mut f64,
    principal_values: *mut f64,
    principal_vectors: *mut f64,
    nodal_von_mises: *mut f64,
) -> i32 {
    spr_ffi_guard(AssertUnwindSafe(|| {
        // Marshal node positions
        let pos_slice = slice::from_raw_parts(node_positions, num_nodes * 3);
        let positions = Array2::from_shape_vec((num_nodes, 3), pos_slice.to_vec())
            .map_err(|e| TheseusError::Shape(format!("node_positions: {e}")))?;

        // Marshal elements
        let elem_slice = slice::from_raw_parts(elements, num_elements * 4);
        let elems: Vec<[usize; 4]> = (0..num_elements)
            .map(|e| [
                elem_slice[e * 4] as usize,
                elem_slice[e * 4 + 1] as usize,
                elem_slice[e * 4 + 2] as usize,
                elem_slice[e * 4 + 3] as usize,
            ])
            .collect();

        // Marshal element stresses
        let stress_slice = slice::from_raw_parts(element_stresses, num_elements * 6);
        let elem_stresses: Vec<[f64; 6]> = (0..num_elements)
            .map(|e| {
                let mut s = [0.0f64; 6];
                s.copy_from_slice(&stress_slice[e * 6..(e + 1) * 6]);
                s
            })
            .collect();

        // Run SPR recovery
        let (nodal_s, elem_err) = spr_recover_nodal_stresses(&positions, &elems, &elem_stresses);

        // Write nodal stresses
        let out_nodal = slice::from_raw_parts_mut(nodal_stresses, num_nodes * 6);
        for (i, s) in nodal_s.iter().enumerate() {
            for c in 0..6 {
                out_nodal[i * 6 + c] = s[c];
            }
        }

        // Write element errors
        let out_errors = slice::from_raw_parts_mut(element_errors, num_elements);
        out_errors.copy_from_slice(&elem_err);

        // Compute principal stresses and von Mises for each node
        let out_pvals = slice::from_raw_parts_mut(principal_values, num_nodes * 3);
        let out_pvecs = slice::from_raw_parts_mut(principal_vectors, num_nodes * 9);
        let out_vm = slice::from_raw_parts_mut(nodal_von_mises, num_nodes);

        for (i, s) in nodal_s.iter().enumerate() {
            let (vals, vecs) = principal_stresses(s);

            out_pvals[i * 3] = vals[0];
            out_pvals[i * 3 + 1] = vals[1];
            out_pvals[i * 3 + 2] = vals[2];

            for k in 0..3 {
                for d in 0..3 {
                    out_pvecs[i * 9 + k * 3 + d] = vecs[k][d];
                }
            }

            out_vm[i] = von_mises_from_principals(vals[0], vals[1], vals[2]);
        }

        Ok(())
    }))
}
