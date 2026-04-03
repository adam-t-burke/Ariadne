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
/// Accepts per-GP stresses (4 Gauss points per element) and averages them
/// to one stress per element before running SPR recovery.
///
/// # Safety
/// All pointers must be valid for the given lengths:
/// - `node_positions`: `num_nodes * 3` f64
/// - `elements`: `num_elements * 4` i32
/// - `element_stresses`: `num_elements * NUM_GP * 6` f64 (4 GPs × 6 components)
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
    num_nodes_per_element: *const i32,
    element_stresses: *const f64,
    nodal_stresses: *mut f64,
    element_errors: *mut f64,
    principal_values: *mut f64,
    principal_vectors: *mut f64,
    nodal_von_mises: *mut f64,
) -> i32 {
    let ngp = crate::solid_assembly::NUM_GP;

    spr_ffi_guard(AssertUnwindSafe(|| {
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

        // Average 4 GP stresses per element down to 1 for SPR
        let stress_slice = slice::from_raw_parts(element_stresses, num_elements * ngp * 6);
        let elem_stresses: Vec<[f64; 6]> = (0..num_elements)
            .map(|e| {
                let mut avg = [0.0f64; 6];
                for gp in 0..ngp {
                    for c in 0..6 {
                        avg[c] += stress_slice[e * ngp * 6 + gp * 6 + c];
                    }
                }
                let inv = 1.0 / ngp as f64;
                for c in 0..6 {
                    avg[c] *= inv;
                }
                avg
            })
            .collect();

        let (nodal_s, elem_err) = spr_recover_nodal_stresses(&positions, &elems, &elem_stresses);

        let out_nodal = slice::from_raw_parts_mut(nodal_stresses, num_nodes * 6);
        for (i, s) in nodal_s.iter().enumerate() {
            for c in 0..6 {
                out_nodal[i * 6 + c] = s[c];
            }
        }

        let out_errors = slice::from_raw_parts_mut(element_errors, num_elements);
        out_errors.copy_from_slice(&elem_err);

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
