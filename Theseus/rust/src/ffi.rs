//! C-compatible FFI for Grasshopper (C# P/Invoke).
//!
//! # Architecture
//!
//! The FFI layer is deliberately thin — each `extern "C"` function does
//! only two things:
//!
//!   1. **Marshal** raw pointers / lengths into safe Rust types.
//!   2. **Delegate** to a safe, `Result`-returning inner function.
//!
//! All real logic lives in the core library (`types`, `fdm`, `gradients`,
//! `optimizer`), which never panics.  Errors propagate as
//! `Result<_, TheseusError>` and are translated at the FFI boundary to:
//!
//!   - `i32` return codes: 0 = success, negative = error.
//!   - Thread-local error message retrievable via `theseus_last_error`.
//!
//! `catch_unwind` wraps every `extern "C"` as a **safety net only** — if it
//! ever fires, that means we have a bug (an uncovered panic path in the
//! core).  It exists solely to prevent UB from stack unwinding across FFI.
//!
//! # Memory convention
//!
//!   - Caller allocates flat arrays and passes pointers + lengths.
//!   - Opaque handles (`*mut TheseusHandle`) are created by Rust and freed
//!     by Rust via `theseus_free`.
//!   - No JSON, no WebSocket — pure value types over the boundary.

use crate::types::*;
use crate::optimizer;
use ndarray::Array2;
use sprs::{CsMat, TriMat};
use std::cell::RefCell;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::slice;

// ─────────────────────────────────────────────────────────────
//  Thread-local error message  (the SQLite pattern)
// ─────────────────────────────────────────────────────────────

thread_local! {
    static LAST_ERROR: RefCell<String> = RefCell::new(String::new());
}

/// Store an error message for later retrieval by `theseus_last_error`.
fn set_last_error(msg: &str) {
    LAST_ERROR.with(|e| *e.borrow_mut() = msg.to_owned());
}

/// Wrap an `extern "C"` body: calls the closure, translates `Result` to
/// `i32`, stores error message, and uses `catch_unwind` as a final safety
/// net against bugs.
unsafe fn ffi_guard<F>(f: F) -> i32
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
            set_last_error("internal panic (this is a bug — please report it)");
            -2
        }
    }
}

/// Retrieve the last error message.
///
/// Copies the UTF-8 message into a caller-provided buffer.  Returns the
/// number of bytes written (excluding null terminator), or −1 if the
/// buffer is too small.  A return of 0 means no error has been recorded.
///
/// # Safety
/// `buf` must point to at least `buf_len` writable bytes.
#[no_mangle]
pub unsafe extern "C" fn theseus_last_error(buf: *mut u8, buf_len: usize) -> i32 {
    LAST_ERROR.with(|e| {
        let msg = e.borrow();
        if msg.is_empty() {
            return 0;
        }
        let bytes = msg.as_bytes();
        if buf_len < bytes.len() + 1 {
            return -1; // buffer too small
        }
        let out = slice::from_raw_parts_mut(buf, buf_len);
        out[..bytes.len()].copy_from_slice(bytes);
        out[bytes.len()] = 0; // null terminator
        bytes.len() as i32
    })
}

// ─────────────────────────────────────────────────────────────
//  Opaque handle
// ─────────────────────────────────────────────────────────────

/// C-callable progress callback.
///
/// Called every `report_frequency` evaluations during optimization with:
///   - `iteration`: evaluation count (0-based)
///   - `loss`: current objective value
///   - `xyz`: pointer to `num_nodes * 3` doubles (row-major node positions)
///   - `num_nodes`: total number of nodes
pub type ProgressCallback = unsafe extern "C" fn(
    iteration: usize,
    loss: f64,
    xyz: *const f64,
    num_nodes: usize,
);

/// Solver handle that owns the problem + state.
pub struct TheseusHandle {
    pub problem: Problem,
    pub state: OptimizationState,
    pub progress_callback: Option<ProgressCallback>,
    pub report_frequency: usize,
}

// ─────────────────────────────────────────────────────────────
//  Problem construction
// ─────────────────────────────────────────────────────────────

/// Create a new problem from raw arrays.
///
/// Returns a valid handle pointer on success, or null on failure.
/// On failure call `theseus_last_error` for details.
///
/// # Safety
/// All pointers must be valid for the given lengths.
#[no_mangle]
pub unsafe extern "C" fn theseus_create(
    // ── Incidence (COO triplets) ──
    num_edges: usize,
    num_nodes: usize,
    num_free: usize,
    coo_rows: *const usize,
    coo_cols: *const usize,
    coo_vals: *const f64,
    coo_nnz: usize,
    free_node_indices: *const usize,
    fixed_node_indices: *const usize,
    num_fixed: usize,
    // ── Loads & geometry ──
    loads: *const f64,
    fixed_positions: *const f64,
    // ── Initial q ──
    q_init: *const f64,
    // ── Bounds ──
    lower_bounds: *const f64,
    upper_bounds: *const f64,
) -> *mut TheseusHandle {
    let result = catch_unwind(AssertUnwindSafe(|| {
        create_inner(
            num_edges, num_nodes, num_free,
            coo_rows, coo_cols, coo_vals, coo_nnz,
            free_node_indices, fixed_node_indices, num_fixed,
            loads, fixed_positions,
            q_init, lower_bounds, upper_bounds,
        )
    }));

    match result {
        Ok(Ok(ptr)) => ptr,
        Ok(Err(e)) => {
            set_last_error(&e.to_string());
            std::ptr::null_mut()
        }
        Err(_panic) => {
            set_last_error("internal panic in theseus_create (this is a bug)");
            std::ptr::null_mut()
        }
    }
}

/// Safe inner function for `theseus_create`.  All pointer-to-slice
/// conversion happens here; everything downstream is pure safe Rust.
unsafe fn create_inner(
    num_edges: usize, num_nodes: usize, num_free: usize,
    coo_rows: *const usize, coo_cols: *const usize, coo_vals: *const f64, coo_nnz: usize,
    free_node_indices: *const usize, fixed_node_indices: *const usize, num_fixed: usize,
    loads: *const f64, fixed_positions: *const f64,
    q_init: *const f64, lower_bounds: *const f64, upper_bounds: *const f64,
) -> Result<*mut TheseusHandle, TheseusError> {
    let rows = slice::from_raw_parts(coo_rows, coo_nnz);
    let cols = slice::from_raw_parts(coo_cols, coo_nnz);
    let vals = slice::from_raw_parts(coo_vals, coo_nnz);
    let free_idx = slice::from_raw_parts(free_node_indices, num_free).to_vec();
    let fixed_idx = slice::from_raw_parts(fixed_node_indices, num_fixed).to_vec();
    let loads_slice = slice::from_raw_parts(loads, num_free * 3);
    let fixed_pos_slice = slice::from_raw_parts(fixed_positions, num_fixed * 3);
    let q_slice = slice::from_raw_parts(q_init, num_edges);
    let lb_slice = slice::from_raw_parts(lower_bounds, num_edges);
    let ub_slice = slice::from_raw_parts(upper_bounds, num_edges);

    // Build incidence matrix from COO
    let mut tri = TriMat::new((num_edges, num_nodes));
    for i in 0..coo_nnz {
        tri.add_triplet(rows[i], cols[i], vals[i]);
    }
    let incidence = tri.to_csc();

    let free_inc = extract_columns(&incidence, &free_idx);
    let fixed_inc = extract_columns(&incidence, &fixed_idx);

    let topology = NetworkTopology {
        incidence,
        free_incidence: free_inc,
        fixed_incidence: fixed_inc,
        num_edges,
        num_nodes,
        free_node_indices: free_idx,
        fixed_node_indices: fixed_idx.clone(),
    };

    let free_node_loads = Array2::from_shape_vec((num_free, 3), loads_slice.to_vec())
        .map_err(|e| TheseusError::Shape(format!("loads: {e}")))?;

    let fixed_node_positions = Array2::from_shape_vec((num_fixed, 3), fixed_pos_slice.to_vec())
        .map_err(|e| TheseusError::Shape(format!("fixed_positions: {e}")))?;

    let anchors = AnchorInfo::all_fixed(fixed_node_positions.clone());

    let bounds = Bounds {
        lower: lb_slice.to_vec(),
        upper: ub_slice.to_vec(),
    };

    let problem = Problem {
        topology,
        free_node_loads,
        fixed_node_positions,
        anchors,
        objectives: Vec::new(),
        bounds,
        solver: SolverOptions::default(),
    };

    let state = OptimizationState::new(q_slice.to_vec(), Array2::zeros((0, 3)));

    Ok(Box::into_raw(Box::new(TheseusHandle {
        problem,
        state,
        progress_callback: None,
        report_frequency: 1,
    })))
}

/// Free a handle.
///
/// # Safety
/// `handle` must be a pointer returned by `theseus_create`, or null.
#[no_mangle]
pub unsafe extern "C" fn theseus_free(handle: *mut TheseusHandle) {
    if handle.is_null() { return; }
    let _ = catch_unwind(AssertUnwindSafe(|| {
        drop(Box::from_raw(handle));
    }));
}

// ─────────────────────────────────────────────────────────────
//  Objective registration
// ─────────────────────────────────────────────────────────────

/// Add a TargetXYZ objective.  Returns 0 on success.
///
/// # Safety
/// Valid handle and arrays.
#[no_mangle]
pub unsafe extern "C" fn theseus_add_target_xyz(
    handle: *mut TheseusHandle,
    weight: f64,
    node_indices: *const usize,
    num_nodes: usize,
    target_xyz: *const f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let idx = slice::from_raw_parts(node_indices, num_nodes).to_vec();
        let target = Array2::from_shape_vec(
            (num_nodes, 3),
            slice::from_raw_parts(target_xyz, num_nodes * 3).to_vec(),
        ).map_err(|e| TheseusError::Shape(format!("target_xyz: {e}")))?;
        h.problem.objectives.push(Box::new(TargetXYZ { weight, node_indices: idx, target }));
        Ok(())
    }))
}

/// Add a TargetLength objective.  Returns 0 on success.
///
/// # Safety
/// Valid handle and arrays.
#[no_mangle]
pub unsafe extern "C" fn theseus_add_target_length(
    handle: *mut TheseusHandle,
    weight: f64,
    edge_indices: *const usize,
    num_edges: usize,
    targets: *const f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let idx = slice::from_raw_parts(edge_indices, num_edges).to_vec();
        let tgt = slice::from_raw_parts(targets, num_edges).to_vec();
        h.problem.objectives.push(Box::new(TargetLength { weight, edge_indices: idx, target: tgt }));
        Ok(())
    }))
}

/// Add a MinLength barrier objective.  Returns 0 on success.
///
/// # Safety
/// Valid handle and arrays.
#[no_mangle]
pub unsafe extern "C" fn theseus_add_min_length(
    handle: *mut TheseusHandle,
    weight: f64,
    edge_indices: *const usize,
    num_edges: usize,
    thresholds: *const f64,
    sharpness: f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let idx = slice::from_raw_parts(edge_indices, num_edges).to_vec();
        let thr = slice::from_raw_parts(thresholds, num_edges).to_vec();
        h.problem.objectives.push(Box::new(MinLength {
            weight, edge_indices: idx, threshold: thr, sharpness,
        }));
        Ok(())
    }))
}

/// Add a TargetXY objective (XY plane only).  Returns 0 on success.
///
/// `target_xy` is a flat row-major `num_nodes × 3` array (Z column is ignored
/// by the loss but must be present for uniform layout).
///
/// # Safety
/// Valid handle and arrays.
#[no_mangle]
pub unsafe extern "C" fn theseus_add_target_xy(
    handle: *mut TheseusHandle,
    weight: f64,
    node_indices: *const usize,
    num_nodes: usize,
    target_xy: *const f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let idx = slice::from_raw_parts(node_indices, num_nodes).to_vec();
        let target = Array2::from_shape_vec(
            (num_nodes, 3),
            slice::from_raw_parts(target_xy, num_nodes * 3).to_vec(),
        ).map_err(|e| TheseusError::Shape(format!("target_xy: {e}")))?;
        h.problem.objectives.push(Box::new(TargetXY { weight, node_indices: idx, target }));
        Ok(())
    }))
}

/// Add a LengthVariation objective (minimise range of edge lengths).
///
/// # Safety
/// Valid handle and arrays.
#[no_mangle]
pub unsafe extern "C" fn theseus_add_length_variation(
    handle: *mut TheseusHandle,
    weight: f64,
    edge_indices: *const usize,
    num_edges: usize,
    sharpness: f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let idx = slice::from_raw_parts(edge_indices, num_edges).to_vec();
        h.problem.objectives.push(Box::new(LengthVariation {
            weight, edge_indices: idx, sharpness,
        }));
        Ok(())
    }))
}

/// Add a ForceVariation objective (minimise range of member forces).
///
/// # Safety
/// Valid handle and arrays.
#[no_mangle]
pub unsafe extern "C" fn theseus_add_force_variation(
    handle: *mut TheseusHandle,
    weight: f64,
    edge_indices: *const usize,
    num_edges: usize,
    sharpness: f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let idx = slice::from_raw_parts(edge_indices, num_edges).to_vec();
        h.problem.objectives.push(Box::new(ForceVariation {
            weight, edge_indices: idx, sharpness,
        }));
        Ok(())
    }))
}

/// Add a SumForceLength objective (minimise Σ |f_k| × ℓ_k).
///
/// # Safety
/// Valid handle and arrays.
#[no_mangle]
pub unsafe extern "C" fn theseus_add_sum_force_length(
    handle: *mut TheseusHandle,
    weight: f64,
    edge_indices: *const usize,
    num_edges: usize,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let idx = slice::from_raw_parts(edge_indices, num_edges).to_vec();
        h.problem.objectives.push(Box::new(SumForceLength {
            weight, edge_indices: idx,
        }));
        Ok(())
    }))
}

/// Add a MaxLength barrier objective (penalty for edges exceeding threshold).
///
/// # Safety
/// Valid handle and arrays.
#[no_mangle]
pub unsafe extern "C" fn theseus_add_max_length(
    handle: *mut TheseusHandle,
    weight: f64,
    edge_indices: *const usize,
    num_edges: usize,
    thresholds: *const f64,
    sharpness: f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let idx = slice::from_raw_parts(edge_indices, num_edges).to_vec();
        let thr = slice::from_raw_parts(thresholds, num_edges).to_vec();
        h.problem.objectives.push(Box::new(MaxLength {
            weight, edge_indices: idx, threshold: thr, sharpness,
        }));
        Ok(())
    }))
}

/// Add a MinForce barrier objective (penalty for forces below threshold).
///
/// # Safety
/// Valid handle and arrays.
#[no_mangle]
pub unsafe extern "C" fn theseus_add_min_force(
    handle: *mut TheseusHandle,
    weight: f64,
    edge_indices: *const usize,
    num_edges: usize,
    thresholds: *const f64,
    sharpness: f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let idx = slice::from_raw_parts(edge_indices, num_edges).to_vec();
        let thr = slice::from_raw_parts(thresholds, num_edges).to_vec();
        h.problem.objectives.push(Box::new(MinForce {
            weight, edge_indices: idx, threshold: thr, sharpness,
        }));
        Ok(())
    }))
}

/// Add a MaxForce barrier objective (penalty for forces exceeding threshold).
///
/// # Safety
/// Valid handle and arrays.
#[no_mangle]
pub unsafe extern "C" fn theseus_add_max_force(
    handle: *mut TheseusHandle,
    weight: f64,
    edge_indices: *const usize,
    num_edges: usize,
    thresholds: *const f64,
    sharpness: f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let idx = slice::from_raw_parts(edge_indices, num_edges).to_vec();
        let thr = slice::from_raw_parts(thresholds, num_edges).to_vec();
        h.problem.objectives.push(Box::new(MaxForce {
            weight, edge_indices: idx, threshold: thr, sharpness,
        }));
        Ok(())
    }))
}

/// Add a RigidSetCompare objective (compare pairwise distances of a node set
/// against target positions).
///
/// `target_xyz` is a flat row-major `num_nodes × 3` array.
///
/// # Safety
/// Valid handle and arrays.
#[no_mangle]
pub unsafe extern "C" fn theseus_add_rigid_set_compare(
    handle: *mut TheseusHandle,
    weight: f64,
    node_indices: *const usize,
    num_nodes: usize,
    target_xyz: *const f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let idx = slice::from_raw_parts(node_indices, num_nodes).to_vec();
        let target = Array2::from_shape_vec(
            (num_nodes, 3),
            slice::from_raw_parts(target_xyz, num_nodes * 3).to_vec(),
        ).map_err(|e| TheseusError::Shape(format!("rigid_set target: {e}")))?;
        h.problem.objectives.push(Box::new(RigidSetCompare { weight, node_indices: idx, target }));
        Ok(())
    }))
}

/// Add a ReactionDirection objective (align anchor reaction directions).
///
/// `target_dirs` is a flat row-major `num_anchors × 3` array of unit vectors.
///
/// # Safety
/// Valid handle and arrays.
#[no_mangle]
pub unsafe extern "C" fn theseus_add_reaction_direction(
    handle: *mut TheseusHandle,
    weight: f64,
    anchor_indices: *const usize,
    num_anchors: usize,
    target_dirs: *const f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let idx = slice::from_raw_parts(anchor_indices, num_anchors).to_vec();
        let dirs = Array2::from_shape_vec(
            (num_anchors, 3),
            slice::from_raw_parts(target_dirs, num_anchors * 3).to_vec(),
        ).map_err(|e| TheseusError::Shape(format!("reaction_direction targets: {e}")))?;
        h.problem.objectives.push(Box::new(ReactionDirection {
            weight, anchor_indices: idx, target_directions: dirs,
        }));
        Ok(())
    }))
}

/// Add a ReactionDirectionMagnitude objective (align anchor reactions in both
/// direction and magnitude).
///
/// `target_dirs` is a flat row-major `num_anchors × 3` array of unit vectors.
/// `target_mags` is a flat `num_anchors`-element array of target magnitudes.
///
/// # Safety
/// Valid handle and arrays.
#[no_mangle]
pub unsafe extern "C" fn theseus_add_reaction_direction_magnitude(
    handle: *mut TheseusHandle,
    weight: f64,
    anchor_indices: *const usize,
    num_anchors: usize,
    target_dirs: *const f64,
    target_mags: *const f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let idx = slice::from_raw_parts(anchor_indices, num_anchors).to_vec();
        let dirs = Array2::from_shape_vec(
            (num_anchors, 3),
            slice::from_raw_parts(target_dirs, num_anchors * 3).to_vec(),
        ).map_err(|e| TheseusError::Shape(format!("reaction_dir_mag targets: {e}")))?;
        let mags = slice::from_raw_parts(target_mags, num_anchors).to_vec();
        h.problem.objectives.push(Box::new(ReactionDirectionMagnitude {
            weight, anchor_indices: idx, target_directions: dirs, target_magnitudes: mags,
        }));
        Ok(())
    }))
}

/// Configure solver options.  Returns 0 on success.
///
/// # Safety
/// Valid handle.
#[no_mangle]
pub unsafe extern "C" fn theseus_set_solver_options(
    handle: *mut TheseusHandle,
    max_iterations: usize,
    abs_tol: f64,
    rel_tol: f64,
    barrier_weight: f64,
    barrier_sharpness: f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        h.problem.solver = SolverOptions {
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
//  Progress callback
// ─────────────────────────────────────────────────────────────

/// Register a progress callback invoked every `frequency` evaluations.
///
/// Pass a null function pointer to clear the callback.
///
/// # Safety
/// Valid handle.  The callback pointer must remain valid for the
/// lifetime of any subsequent `theseus_optimize` call.
#[no_mangle]
pub unsafe extern "C" fn theseus_set_progress_callback(
    handle: *mut TheseusHandle,
    callback: Option<ProgressCallback>,
    frequency: usize,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        h.progress_callback = callback;
        h.report_frequency = if frequency == 0 { 1 } else { frequency };
        Ok(())
    }))
}

// ─────────────────────────────────────────────────────────────
//  Run optimisation
// ─────────────────────────────────────────────────────────────

/// Run L-BFGS optimisation.  Results are written into caller-provided buffers.
///
/// Returns 0 on success, -1 on error (call `theseus_last_error` for details),
/// -2 on internal panic (a bug).
///
/// # Safety
/// All output buffers must have the correct sizes.
#[no_mangle]
pub unsafe extern "C" fn theseus_optimize(
    handle: *mut TheseusHandle,
    out_xyz: *mut f64,
    out_lengths: *mut f64,
    out_forces: *mut f64,
    out_q: *mut f64,
    out_reactions: *mut f64,
    out_iterations: *mut usize,
    out_converged: *mut bool,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let cb = h.progress_callback;
        let freq = h.report_frequency;
        let result = optimizer::optimize(&h.problem, &mut h.state, cb, freq)?;

        let nn = h.problem.topology.num_nodes;
        let ne = h.problem.topology.num_edges;

        // Copy xyz
        let xyz_out = slice::from_raw_parts_mut(out_xyz, nn * 3);
        for i in 0..nn {
            for d in 0..3 {
                xyz_out[i * 3 + d] = result.xyz[[i, d]];
            }
        }

        // Copy lengths, forces, q
        slice::from_raw_parts_mut(out_lengths, ne).copy_from_slice(&result.member_lengths);
        slice::from_raw_parts_mut(out_forces, ne).copy_from_slice(&result.member_forces);
        slice::from_raw_parts_mut(out_q, ne).copy_from_slice(&result.q);

        // Copy reactions
        let r_out = slice::from_raw_parts_mut(out_reactions, nn * 3);
        for i in 0..nn {
            for d in 0..3 {
                r_out[i * 3 + d] = result.reactions[[i, d]];
            }
        }

        *out_iterations = result.iterations;
        *out_converged = result.converged;

        Ok(())
    }))
}

// ─────────────────────────────────────────────────────────────
//  Forward solve only  (no optimisation)
// ─────────────────────────────────────────────────────────────

/// Single forward FDM solve — useful for previewing geometry without optimising.
///
/// Returns 0 on success, -1 on error, -2 on internal panic.
///
/// # Safety
/// Valid handle and output buffers.
#[no_mangle]
pub unsafe extern "C" fn theseus_solve_forward(
    handle: *mut TheseusHandle,
    out_xyz: *mut f64,
    out_lengths: *mut f64,
    out_forces: *mut f64,
    out_q: *mut f64,
    out_reactions: *mut f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = &mut *handle;
        let mut cache = FdmCache::new(&h.problem)?;
        let anchors = h.state.variable_anchor_positions.clone();

        crate::fdm::solve_fdm(&mut cache, &h.state.force_densities, &h.problem, &anchors, 1e-12)?;

        let nn = h.problem.topology.num_nodes;
        let ne = h.problem.topology.num_edges;

        let xyz_out = slice::from_raw_parts_mut(out_xyz, nn * 3);
        for i in 0..nn {
            for d in 0..3 {
                xyz_out[i * 3 + d] = cache.nf[[i, d]];
            }
        }
        slice::from_raw_parts_mut(out_lengths, ne).copy_from_slice(&cache.member_lengths);
        slice::from_raw_parts_mut(out_forces, ne).copy_from_slice(&cache.member_forces);
        slice::from_raw_parts_mut(out_q, ne).copy_from_slice(&h.state.force_densities);

        let r_out = slice::from_raw_parts_mut(out_reactions, nn * 3);
        for i in 0..nn {
            for d in 0..3 {
                r_out[i * 3 + d] = cache.reactions[[i, d]];
            }
        }

        Ok(())
    }))
}


// ─────────────────────────────────────────────────────────────
//  Helpers
// ─────────────────────────────────────────────────────────────

/// Extract columns from a CSC matrix by index.
fn extract_columns(mat: &CsMat<f64>, cols: &[usize]) -> CsMat<f64> {
    let nrows = mat.rows();
    let ncols = cols.len();
    let mut tri = TriMat::new((nrows, ncols));

    let mat_csc = mat.to_csc();
    for (new_col, &old_col) in cols.iter().enumerate() {
        let start = mat_csc.indptr().raw_storage()[old_col];
        let end_ = mat_csc.indptr().raw_storage()[old_col + 1];
        for nz in start..end_ {
            tri.add_triplet(mat_csc.indices()[nz], new_col, mat_csc.data()[nz]);
        }
    }

    tri.to_csc()
}
