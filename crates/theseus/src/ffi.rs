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
//!   - No JSON — pure value types over the boundary.

use crate::optimizer;
use crate::sparse::SparseColMatOwned;
use crate::types::*;
use ndarray::Array2;
use std::cell::RefCell;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::slice;
use std::sync::atomic::{AtomicBool, Ordering};

// ─────────────────────────────────────────────────────────────
//  Thread-local error message  (the SQLite pattern)
// ─────────────────────────────────────────────────────────────

thread_local! {
    static LAST_ERROR: RefCell<String> = RefCell::new(String::new());
}

/// Store an error message for later retrieval by `theseus_last_error`.
pub(crate) fn set_last_error(msg: &str) {
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

/// Require a non-null mutable handle pointer.
unsafe fn require_handle(
    handle: *mut TheseusHandle,
) -> Result<&'static mut TheseusHandle, TheseusError> {
    if handle.is_null() {
        Err(TheseusError::Shape("null TheseusHandle".into()))
    } else {
        Ok(&mut *handle)
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
/// Called every `report_frequency` accepted L-BFGS iterations during optimization with:
///   - `iteration`: completed major iteration count (1-based)
///   - `loss`: current objective value
///   - `xyz`: pointer to `num_nodes * 3` doubles (row-major node positions)
///   - `num_nodes`: total number of nodes
///
/// Returns `1` to continue optimization, `0` to cancel.
pub type ProgressCallback = unsafe extern "C" fn(
    iteration: usize,
    loss: f64,
    xyz: *const f64,
    num_nodes: usize,
    q: *const f64,
    num_edges: usize,
) -> u8;

/// Solver handle that owns the problem + state.
pub struct TheseusHandle {
    pub problem: Problem,
    pub state: OptimizationState,
    pub progress_callback: Option<ProgressCallback>,
    pub report_frequency: usize,
    pub last_termination_reason: String,
    pub cancel_flag: AtomicBool,
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
            num_edges,
            num_nodes,
            num_free,
            coo_rows,
            coo_cols,
            coo_vals,
            coo_nnz,
            free_node_indices,
            fixed_node_indices,
            num_fixed,
            loads,
            fixed_positions,
            q_init,
            lower_bounds,
            upper_bounds,
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

/// Create a new problem with variable support definitions.
///
/// `support_kinds`: 0 = Sphere, 1 = Roller, 2 = Rail, 3 = NURBS curve, 4 = NURBS surface.
///
/// For each support i:
/// - `variable_node_indices[i]` is the global node index.
/// - Sphere uses `sphere_radii[i]`.
/// - Roller uses `roller_enabled[i*3..i*3+3]`, `roller_lower[...]`, `roller_upper[...]`.
/// - Rail uses `rail_start[i*3..]` and `rail_end[i*3..]`.
#[no_mangle]
pub unsafe extern "C" fn theseus_create_with_variable_supports(
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
    // ── Variable supports ──
    num_variable_supports: usize,
    variable_node_indices: *const usize,
    support_kinds: *const i32,
    sphere_radii: *const f64,
    roller_enabled: *const u8,
    roller_lower: *const f64,
    roller_upper: *const f64,
    rail_start: *const f64,
    rail_end: *const f64,
    nurbs_offsets: *const usize,
    nurbs_lengths: *const usize,
    nurbs_data: *const f64,
    nurbs_data_len: usize,
) -> *mut TheseusHandle {
    let result = catch_unwind(AssertUnwindSafe(|| {
        create_inner_with_variable_supports(
            num_edges,
            num_nodes,
            num_free,
            coo_rows,
            coo_cols,
            coo_vals,
            coo_nnz,
            free_node_indices,
            fixed_node_indices,
            num_fixed,
            loads,
            fixed_positions,
            q_init,
            lower_bounds,
            upper_bounds,
            num_variable_supports,
            variable_node_indices,
            support_kinds,
            sphere_radii,
            roller_enabled,
            roller_lower,
            roller_upper,
            rail_start,
            rail_end,
            nurbs_offsets,
            nurbs_lengths,
            nurbs_data,
            nurbs_data_len,
        )
    }));

    match result {
        Ok(Ok(ptr)) => ptr,
        Ok(Err(e)) => {
            set_last_error(&e.to_string());
            std::ptr::null_mut()
        }
        Err(_panic) => {
            set_last_error(
                "internal panic in theseus_create_with_variable_supports (this is a bug)",
            );
            std::ptr::null_mut()
        }
    }
}

/// Safe inner function for `theseus_create`.  All pointer-to-slice
/// conversion happens here; everything downstream is pure safe Rust.
unsafe fn create_inner(
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
    loads: *const f64,
    fixed_positions: *const f64,
    q_init: *const f64,
    lower_bounds: *const f64,
    upper_bounds: *const f64,
) -> Result<*mut TheseusHandle, TheseusError> {
    create_inner_with_variable_supports(
        num_edges,
        num_nodes,
        num_free,
        coo_rows,
        coo_cols,
        coo_vals,
        coo_nnz,
        free_node_indices,
        fixed_node_indices,
        num_fixed,
        loads,
        fixed_positions,
        q_init,
        lower_bounds,
        upper_bounds,
        0,
        std::ptr::null(),
        std::ptr::null(),
        std::ptr::null(),
        std::ptr::null(),
        std::ptr::null(),
        std::ptr::null(),
        std::ptr::null(),
        std::ptr::null(),
        std::ptr::null(),
        std::ptr::null(),
        std::ptr::null(),
        0,
    )
}

fn nurbs_support_data<'a>(
    index: usize,
    offsets: &[usize],
    lengths: &[usize],
    data: &'a [f64],
) -> Result<&'a [f64], TheseusError> {
    if offsets.len() <= index || lengths.len() <= index {
        return Err(TheseusError::Shape(format!(
            "NURBS support {index} is missing geometry offset/length"
        )));
    }
    let offset = offsets[index];
    let length = lengths[index];
    let end = offset
        .checked_add(length)
        .ok_or_else(|| TheseusError::Shape("NURBS geometry range overflow".into()))?;
    if end > data.len() {
        return Err(TheseusError::Shape(format!(
            "NURBS support {index} geometry range is out of bounds"
        )));
    }
    Ok(&data[offset..end])
}

fn read_usize(value: f64, label: &str) -> Result<usize, TheseusError> {
    if !value.is_finite() || value < 0.0 || value.fract().abs() > 1e-9 {
        return Err(TheseusError::Shape(format!(
            "{label} must be encoded as a nonnegative integer"
        )));
    }
    Ok(value as usize)
}

fn parse_nurbs_curve_support(data: &[f64]) -> Result<VariableSupportKind, TheseusError> {
    if data.len() < 6 {
        return Err(TheseusError::Shape(
            "NURBS curve support data is too short".into(),
        ));
    }
    let degree = read_usize(data[0], "NURBS curve degree")?;
    let knot_count = read_usize(data[1], "NURBS curve knot count")?;
    let cp_count = read_usize(data[2], "NURBS curve control point count")?;
    let domain = [data[3], data[4]];
    let initial_t = data[5];
    let expected = 6 + knot_count + cp_count * 4;
    if data.len() != expected {
        return Err(TheseusError::Shape(format!(
            "NURBS curve support data length got {}, expected {expected}",
            data.len()
        )));
    }
    let knots = data[6..6 + knot_count].to_vec();
    let cp_start = 6 + knot_count;
    let mut control_points = Vec::with_capacity(cp_count);
    for i in 0..cp_count {
        let j = cp_start + i * 4;
        control_points.push([data[j], data[j + 1], data[j + 2], data[j + 3]]);
    }
    let curve = nurbsbook::NurbsCurve::new(degree, knots, control_points)
        .map_err(|e| TheseusError::Shape(format!("invalid NURBS curve support: {e}")))?;
    Ok(VariableSupportKind::NurbsCurve {
        curve,
        domain,
        initial_t,
    })
}

fn parse_nurbs_surface_support(data: &[f64]) -> Result<VariableSupportKind, TheseusError> {
    if data.len() < 12 {
        return Err(TheseusError::Shape(
            "NURBS surface support data is too short".into(),
        ));
    }
    let degree_u = read_usize(data[0], "NURBS surface U degree")?;
    let degree_v = read_usize(data[1], "NURBS surface V degree")?;
    let count_u = read_usize(data[2], "NURBS surface U control point count")?;
    let count_v = read_usize(data[3], "NURBS surface V control point count")?;
    let knot_u_count = read_usize(data[4], "NURBS surface U knot count")?;
    let knot_v_count = read_usize(data[5], "NURBS surface V knot count")?;
    let domain_u = [data[6], data[7]];
    let domain_v = [data[8], data[9]];
    let initial_uv = [data[10], data[11]];
    let cp_count = count_u
        .checked_mul(count_v)
        .ok_or_else(|| TheseusError::Shape("NURBS surface control point count overflow".into()))?;
    let expected = 12 + knot_u_count + knot_v_count + cp_count * 4;
    if data.len() != expected {
        return Err(TheseusError::Shape(format!(
            "NURBS surface support data length got {}, expected {expected}",
            data.len()
        )));
    }
    let ku_start = 12;
    let kv_start = ku_start + knot_u_count;
    let cp_start = kv_start + knot_v_count;
    let knots_u = data[ku_start..ku_start + knot_u_count].to_vec();
    let knots_v = data[kv_start..kv_start + knot_v_count].to_vec();
    let mut control_points = Vec::with_capacity(cp_count);
    for i in 0..cp_count {
        let j = cp_start + i * 4;
        control_points.push([data[j], data[j + 1], data[j + 2], data[j + 3]]);
    }
    let surface = nurbsbook::NurbsSurface::new(
        degree_u,
        degree_v,
        knots_u,
        knots_v,
        count_u,
        count_v,
        control_points,
    )
    .map_err(|e| TheseusError::Shape(format!("invalid NURBS surface support: {e}")))?;
    Ok(VariableSupportKind::NurbsSurface {
        surface,
        domain_u,
        domain_v,
        initial_uv,
    })
}

unsafe fn create_inner_with_variable_supports(
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
    loads: *const f64,
    fixed_positions: *const f64,
    q_init: *const f64,
    lower_bounds: *const f64,
    upper_bounds: *const f64,
    num_variable_supports: usize,
    variable_node_indices: *const usize,
    support_kinds: *const i32,
    sphere_radii: *const f64,
    roller_enabled: *const u8,
    roller_lower: *const f64,
    roller_upper: *const f64,
    rail_start: *const f64,
    rail_end: *const f64,
    nurbs_offsets: *const usize,
    nurbs_lengths: *const usize,
    nurbs_data: *const f64,
    nurbs_data_len: usize,
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
    let incidence = SparseColMatOwned::from_coo(num_edges, num_nodes, rows, cols, vals)
        .map_err(|e| TheseusError::Shape(e))?;

    let free_inc = incidence.extract_columns(&free_idx);
    let fixed_inc = incidence.extract_columns(&fixed_idx);

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

    let (anchors, latent_init) = if num_variable_supports == 0 {
        (
            AnchorInfo::all_fixed(fixed_node_positions.clone()),
            Vec::new(),
        )
    } else {
        let var_nodes =
            slice::from_raw_parts(variable_node_indices, num_variable_supports).to_vec();
        let kinds = slice::from_raw_parts(support_kinds, num_variable_supports);
        let radii = slice::from_raw_parts(sphere_radii, num_variable_supports);
        let roller_en = slice::from_raw_parts(roller_enabled, num_variable_supports * 3);
        let roller_lb = slice::from_raw_parts(roller_lower, num_variable_supports * 3);
        let roller_ub = slice::from_raw_parts(roller_upper, num_variable_supports * 3);
        let rail_s = slice::from_raw_parts(rail_start, num_variable_supports * 3);
        let rail_e = slice::from_raw_parts(rail_end, num_variable_supports * 3);
        let nurbs_offsets_slice = if nurbs_offsets.is_null() {
            &[][..]
        } else {
            slice::from_raw_parts(nurbs_offsets, num_variable_supports)
        };
        let nurbs_lengths_slice = if nurbs_lengths.is_null() {
            &[][..]
        } else {
            slice::from_raw_parts(nurbs_lengths, num_variable_supports)
        };
        let nurbs_data_slice = if nurbs_data.is_null() {
            &[][..]
        } else {
            slice::from_raw_parts(nurbs_data, nurbs_data_len)
        };

        // Validate variable nodes are fixed support nodes.
        for &node in &var_nodes {
            if !fixed_idx.contains(&node) {
                return Err(TheseusError::Shape(format!(
                    "variable support node index {node} is not in fixed_node_indices"
                )));
            }
        }
        {
            let mut dedup = var_nodes.clone();
            dedup.sort_unstable();
            dedup.dedup();
            if dedup.len() != var_nodes.len() {
                return Err(TheseusError::Shape(
                    "variable support node indices must be unique".into(),
                ));
            }
        }

        let mut supports = Vec::with_capacity(num_variable_supports);
        let mut init_positions = Array2::zeros((num_variable_supports, 3));
        for i in 0..num_variable_supports {
            let node = var_nodes[i];
            let fixed_row = fixed_idx.iter().position(|&n| n == node).ok_or_else(|| {
                TheseusError::Shape(format!(
                    "variable support node index {node} not found in fixed index lookup"
                ))
            })?;
            let x0 = [
                fixed_node_positions[[fixed_row, 0]],
                fixed_node_positions[[fixed_row, 1]],
                fixed_node_positions[[fixed_row, 2]],
            ];
            init_positions[[i, 0]] = x0[0];
            init_positions[[i, 1]] = x0[1];
            init_positions[[i, 2]] = x0[2];

            let kind = match kinds[i] {
                0 => {
                    if !radii[i].is_finite() || radii[i] <= 0.0 {
                        return Err(TheseusError::Shape(format!(
                            "sphere support at index {i} must have positive finite radius"
                        )));
                    }
                    VariableSupportKind::Sphere { radius: radii[i] }
                }
                1 => {
                    let e0 = roller_en[i * 3] != 0;
                    let e1 = roller_en[i * 3 + 1] != 0;
                    let e2 = roller_en[i * 3 + 2] != 0;
                    if !(e0 || e1 || e2) {
                        return Err(TheseusError::Shape(format!(
                            "roller support at index {i} must enable at least one axis"
                        )));
                    }
                    VariableSupportKind::Roller {
                        enabled: [e0, e1, e2],
                        lower: [roller_lb[i * 3], roller_lb[i * 3 + 1], roller_lb[i * 3 + 2]],
                        upper: [roller_ub[i * 3], roller_ub[i * 3 + 1], roller_ub[i * 3 + 2]],
                    }
                }
                2 => VariableSupportKind::Rail {
                    start: [rail_s[i * 3], rail_s[i * 3 + 1], rail_s[i * 3 + 2]],
                    end: [rail_e[i * 3], rail_e[i * 3 + 1], rail_e[i * 3 + 2]],
                },
                3 => {
                    let data = nurbs_support_data(
                        i,
                        nurbs_offsets_slice,
                        nurbs_lengths_slice,
                        nurbs_data_slice,
                    )?;
                    parse_nurbs_curve_support(data)?
                }
                4 => {
                    let data = nurbs_support_data(
                        i,
                        nurbs_offsets_slice,
                        nurbs_lengths_slice,
                        nurbs_data_slice,
                    )?;
                    parse_nurbs_surface_support(data)?
                }
                k => {
                    return Err(TheseusError::Shape(format!(
                        "unknown variable support kind code {k} at index {i}"
                    )));
                }
            };
            supports.push(VariableSupport {
                node_index: node,
                reference_position: x0,
                kind,
            });
        }

        let anchors = AnchorInfo {
            variable_indices: var_nodes,
            fixed_indices: fixed_idx.clone(),
            reference_positions: fixed_node_positions.clone(),
            initial_variable_positions: init_positions,
            variable_supports: supports,
        };
        let latent_init = crate::variable_supports::initial_latents(
            &anchors,
            SolverOptions::default().anchor_saturation_lambda,
        )?;
        (anchors, latent_init)
    };

    let bounds = Bounds {
        lower: lb_slice.to_vec(),
        upper: ub_slice.to_vec(),
    };

    let initial_var_positions = anchors.initial_variable_positions.clone();
    let problem = Problem {
        topology,
        free_node_loads,
        fixed_node_positions,
        anchors,
        objectives: Vec::new(),
        bounds,
        solver: SolverOptions::default(),
        self_weight: None,
        pressure: None,
    };

    let mut state = OptimizationState::new(q_slice.to_vec(), initial_var_positions);
    if !latent_init.is_empty() {
        state.variable_anchor_latents = latent_init.clone();
        state.variable_anchor_positions =
            crate::variable_supports::map_latents_to_positions(&problem, &latent_init)?;
    }

    Ok(Box::into_raw(Box::new(TheseusHandle {
        problem,
        state,
        progress_callback: None,
        report_frequency: 1,
        last_termination_reason: "not run".to_string(),
        cancel_flag: AtomicBool::new(false),
    })))
}

/// Free a handle.
///
/// # Safety
/// `handle` must be a pointer returned by `theseus_create`, or null.
#[no_mangle]
pub unsafe extern "C" fn theseus_free(handle: *mut TheseusHandle) {
    if handle.is_null() {
        return;
    }
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
        let h = require_handle(handle)?;
        let idx = slice::from_raw_parts(node_indices, num_nodes).to_vec();
        let target = Array2::from_shape_vec(
            (num_nodes, 3),
            slice::from_raw_parts(target_xyz, num_nodes * 3).to_vec(),
        )
        .map_err(|e| TheseusError::Shape(format!("target_xyz: {e}")))?;
        h.problem.objectives.push(Box::new(TargetXYZ {
            weight,
            node_indices: idx,
            target,
        }));
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
        let h = require_handle(handle)?;
        let idx = slice::from_raw_parts(edge_indices, num_edges).to_vec();
        let tgt = slice::from_raw_parts(targets, num_edges).to_vec();
        h.problem.objectives.push(Box::new(TargetLength {
            weight,
            edge_indices: idx,
            target: tgt,
        }));
        Ok(())
    }))
}

/// Add a TargetForce objective.  Returns 0 on success.
///
/// # Safety
/// Valid handle and arrays.
#[no_mangle]
pub unsafe extern "C" fn theseus_add_target_force(
    handle: *mut TheseusHandle,
    weight: f64,
    edge_indices: *const usize,
    num_edges: usize,
    targets: *const f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        let idx = slice::from_raw_parts(edge_indices, num_edges).to_vec();
        let tgt = slice::from_raw_parts(targets, num_edges).to_vec();
        h.problem.objectives.push(Box::new(TargetForce {
            weight,
            edge_indices: idx,
            target: tgt,
        }));
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
        let h = require_handle(handle)?;
        let idx = slice::from_raw_parts(edge_indices, num_edges).to_vec();
        let thr = slice::from_raw_parts(thresholds, num_edges).to_vec();
        h.problem.objectives.push(Box::new(MinLength {
            weight,
            edge_indices: idx,
            threshold: thr,
            sharpness,
        }));
        Ok(())
    }))
}

/// Add a TargetXY objective (XY plane only).  Returns 0 on success.
///
/// `target_xy` must be exactly `num_nodes * 3` doubles, row-major (X,Y,Z per node).
/// The Z component is ignored by the loss but must be present so the layout matches
/// TargetXYZ and the FFI does not read past the buffer.
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
        let h = require_handle(handle)?;
        let idx = slice::from_raw_parts(node_indices, num_nodes).to_vec();
        let target = Array2::from_shape_vec(
            (num_nodes, 3),
            slice::from_raw_parts(target_xy, num_nodes * 3).to_vec(),
        )
        .map_err(|e| TheseusError::Shape(format!("target_xy: {e}")))?;
        h.problem.objectives.push(Box::new(TargetXY {
            weight,
            node_indices: idx,
            target,
        }));
        Ok(())
    }))
}

/// Add a TargetPlane objective (projection onto an arbitrary plane).
///
/// `target_xyz` is row-major `num_nodes × 3` world positions.
/// `origin`, `x_axis`, `y_axis` are 3-element arrays in world coordinates;
/// axes should be unit and orthogonal (e.g. Rhino plane Origin, XAxis, YAxis).
///
/// # Safety
/// Valid handle and arrays; origin/x_axis/y_axis must each point to 3 doubles.
#[no_mangle]
pub unsafe extern "C" fn theseus_add_target_plane(
    handle: *mut TheseusHandle,
    weight: f64,
    node_indices: *const usize,
    num_nodes: usize,
    target_xyz: *const f64,
    origin: *const f64,
    x_axis: *const f64,
    y_axis: *const f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        let idx = slice::from_raw_parts(node_indices, num_nodes).to_vec();
        let target = Array2::from_shape_vec(
            (num_nodes, 3),
            slice::from_raw_parts(target_xyz, num_nodes * 3).to_vec(),
        )
        .map_err(|e| TheseusError::Shape(format!("target_plane target_xyz: {e}")))?;
        let origin_arr: [f64; 3] = [*origin.add(0), *origin.add(1), *origin.add(2)];
        let x_axis_arr: [f64; 3] = [*x_axis.add(0), *x_axis.add(1), *x_axis.add(2)];
        let y_axis_arr: [f64; 3] = [*y_axis.add(0), *y_axis.add(1), *y_axis.add(2)];
        h.problem.objectives.push(Box::new(TargetPlane {
            weight,
            node_indices: idx,
            target,
            origin: origin_arr,
            x_axis: x_axis_arr,
            y_axis: y_axis_arr,
        }));
        Ok(())
    }))
}

/// Add a PlanarConstraintAlongDirection objective (pull nodes onto a plane along a direction).
///
/// No target positions — minimizes squared distance along `direction` to the plane.
/// `origin`, `x_axis`, `y_axis`, `direction` are 3-element arrays in world coordinates.
/// Direction must not be parallel to the plane (n·d ≠ 0).
///
/// # Safety
/// Valid handle and arrays; origin/x_axis/y_axis/direction must each point to 3 doubles.
#[no_mangle]
pub unsafe extern "C" fn theseus_add_planar_constraint_along_direction(
    handle: *mut TheseusHandle,
    weight: f64,
    node_indices: *const usize,
    num_nodes: usize,
    origin: *const f64,
    x_axis: *const f64,
    y_axis: *const f64,
    direction: *const f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        let idx = slice::from_raw_parts(node_indices, num_nodes).to_vec();
        let origin_arr: [f64; 3] = [*origin.add(0), *origin.add(1), *origin.add(2)];
        let x_axis_arr: [f64; 3] = [*x_axis.add(0), *x_axis.add(1), *x_axis.add(2)];
        let y_axis_arr: [f64; 3] = [*y_axis.add(0), *y_axis.add(1), *y_axis.add(2)];
        let direction_arr: [f64; 3] = [*direction.add(0), *direction.add(1), *direction.add(2)];
        crate::objectives::planar_constraint_n_dot_d(&x_axis_arr, &y_axis_arr, &direction_arr)?;
        h.problem
            .objectives
            .push(Box::new(PlanarConstraintAlongDirection {
                weight,
                node_indices: idx,
                origin: origin_arr,
                x_axis: x_axis_arr,
                y_axis: y_axis_arr,
                direction: direction_arr,
            }));
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
    use_normalized_variance: u8,
    normalization_strategy: i32,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        let idx = slice::from_raw_parts(edge_indices, num_edges).to_vec();
        let normalization_strategy =
            LengthVarianceNormalizationStrategy::try_from(normalization_strategy)?;
        h.problem.objectives.push(Box::new(LengthVariation {
            weight,
            edge_indices: idx,
            sharpness,
            use_normalized_variance: use_normalized_variance != 0,
            normalization_strategy,
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
    use_normalized_variance: u8,
    normalization_strategy: i32,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        let idx = slice::from_raw_parts(edge_indices, num_edges).to_vec();
        let normalization_strategy =
            ForceVarianceNormalizationStrategy::try_from(normalization_strategy)?;
        h.problem.objectives.push(Box::new(ForceVariation {
            weight,
            edge_indices: idx,
            sharpness,
            use_normalized_variance: use_normalized_variance != 0,
            normalization_strategy,
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
        let h = require_handle(handle)?;
        let idx = slice::from_raw_parts(edge_indices, num_edges).to_vec();
        h.problem.objectives.push(Box::new(SumForceLength {
            weight,
            edge_indices: idx,
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
        let h = require_handle(handle)?;
        let idx = slice::from_raw_parts(edge_indices, num_edges).to_vec();
        let thr = slice::from_raw_parts(thresholds, num_edges).to_vec();
        h.problem.objectives.push(Box::new(MaxLength {
            weight,
            edge_indices: idx,
            threshold: thr,
            sharpness,
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
        let h = require_handle(handle)?;
        let idx = slice::from_raw_parts(edge_indices, num_edges).to_vec();
        let thr = slice::from_raw_parts(thresholds, num_edges).to_vec();
        h.problem.objectives.push(Box::new(MinForce {
            weight,
            edge_indices: idx,
            threshold: thr,
            sharpness,
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
        let h = require_handle(handle)?;
        let idx = slice::from_raw_parts(edge_indices, num_edges).to_vec();
        let thr = slice::from_raw_parts(thresholds, num_edges).to_vec();
        h.problem.objectives.push(Box::new(MaxForce {
            weight,
            edge_indices: idx,
            threshold: thr,
            sharpness,
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
        let h = require_handle(handle)?;
        let idx = slice::from_raw_parts(node_indices, num_nodes).to_vec();
        let target = Array2::from_shape_vec(
            (num_nodes, 3),
            slice::from_raw_parts(target_xyz, num_nodes * 3).to_vec(),
        )
        .map_err(|e| TheseusError::Shape(format!("rigid_set target: {e}")))?;
        h.problem.objectives.push(Box::new(RigidSetCompare {
            weight,
            node_indices: idx,
            target,
        }));
        Ok(())
    }))
}

fn validate_reaction_sign_dirs(
    dirs: &Array2<f64>,
    sign: ReactionMagnitudeSign,
) -> Result<(), TheseusError> {
    if sign != ReactionMagnitudeSign::SignedProjected {
        return Ok(());
    }

    for row in 0..dirs.nrows() {
        let norm =
            (dirs[[row, 0]].powi(2) + dirs[[row, 1]].powi(2) + dirs[[row, 2]].powi(2)).sqrt();
        if norm < f64::EPSILON {
            return Err(TheseusError::Shape(format!(
                "signed projected reaction magnitude requires non-zero target direction at row {row}"
            )));
        }
    }
    Ok(())
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
        let h = require_handle(handle)?;
        let idx = slice::from_raw_parts(anchor_indices, num_anchors).to_vec();
        let dirs = Array2::from_shape_vec(
            (num_anchors, 3),
            slice::from_raw_parts(target_dirs, num_anchors * 3).to_vec(),
        )
        .map_err(|e| TheseusError::Shape(format!("reaction_direction targets: {e}")))?;
        h.problem.objectives.push(Box::new(ReactionDirection {
            weight,
            anchor_indices: idx,
            target_directions: dirs,
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
        let h = require_handle(handle)?;
        let idx = slice::from_raw_parts(anchor_indices, num_anchors).to_vec();
        let dirs = Array2::from_shape_vec(
            (num_anchors, 3),
            slice::from_raw_parts(target_dirs, num_anchors * 3).to_vec(),
        )
        .map_err(|e| TheseusError::Shape(format!("reaction_dir_mag targets: {e}")))?;
        let mags = slice::from_raw_parts(target_mags, num_anchors).to_vec();
        h.problem
            .objectives
            .push(Box::new(ReactionDirectionMagnitude {
                weight,
                anchor_indices: idx,
                target_directions: dirs,
                target_magnitudes: mags,
                behavior: ReactionMagnitudeBehavior::Max,
                sign: ReactionMagnitudeSign::Unsigned,
            }));
        Ok(())
    }))
}

/// Add a ReactionMagnitude objective (constrain anchor reaction magnitudes only).
///
/// `target_dirs` is used as the sign axis when `sign_semantics` is
/// SignedProjected; it is ignored for Unsigned magnitude mode.
///
/// # Safety
/// Valid handle and arrays.
#[no_mangle]
pub unsafe extern "C" fn theseus_add_reaction_magnitude(
    handle: *mut TheseusHandle,
    weight: f64,
    anchor_indices: *const usize,
    num_anchors: usize,
    target_dirs: *const f64,
    target_mags: *const f64,
    behavior: i32,
    sign_semantics: i32,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        let behavior = ReactionMagnitudeBehavior::try_from(behavior)?;
        let sign = ReactionMagnitudeSign::try_from(sign_semantics)?;
        let idx = slice::from_raw_parts(anchor_indices, num_anchors).to_vec();
        let dirs = Array2::from_shape_vec(
            (num_anchors, 3),
            slice::from_raw_parts(target_dirs, num_anchors * 3).to_vec(),
        )
        .map_err(|e| TheseusError::Shape(format!("reaction_magnitude dirs: {e}")))?;
        validate_reaction_sign_dirs(&dirs, sign)?;
        let mags = slice::from_raw_parts(target_mags, num_anchors).to_vec();
        h.problem.objectives.push(Box::new(ReactionMagnitude {
            weight,
            anchor_indices: idx,
            target_directions: dirs,
            target_magnitudes: mags,
            behavior,
            sign,
        }));
        Ok(())
    }))
}

/// Add a ReactionDirectionMagnitude objective with configurable magnitude behavior.
///
/// # Safety
/// Valid handle and arrays.
#[no_mangle]
pub unsafe extern "C" fn theseus_add_reaction_direction_magnitude_with_options(
    handle: *mut TheseusHandle,
    weight: f64,
    anchor_indices: *const usize,
    num_anchors: usize,
    target_dirs: *const f64,
    target_mags: *const f64,
    behavior: i32,
    sign_semantics: i32,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        let behavior = ReactionMagnitudeBehavior::try_from(behavior)?;
        let sign = ReactionMagnitudeSign::try_from(sign_semantics)?;
        let idx = slice::from_raw_parts(anchor_indices, num_anchors).to_vec();
        let dirs = Array2::from_shape_vec(
            (num_anchors, 3),
            slice::from_raw_parts(target_dirs, num_anchors * 3).to_vec(),
        )
        .map_err(|e| TheseusError::Shape(format!("reaction_dir_mag targets: {e}")))?;
        validate_reaction_sign_dirs(&dirs, sign)?;
        let mags = slice::from_raw_parts(target_mags, num_anchors).to_vec();
        h.problem
            .objectives
            .push(Box::new(ReactionDirectionMagnitude {
                weight,
                anchor_indices: idx,
                target_directions: dirs,
                target_magnitudes: mags,
                behavior,
                sign,
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
    anchor_saturation_lambda: f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        h.problem.solver.max_iterations = max_iterations;
        h.problem.solver.absolute_tolerance = abs_tol;
        h.problem.solver.relative_tolerance = rel_tol;
        h.problem.solver.report_frequency = 1;
        h.problem.solver.barrier_weight = barrier_weight;
        h.problem.solver.barrier_sharpness = barrier_sharpness;
        h.problem.solver.anchor_saturation_lambda = anchor_saturation_lambda.max(1e-12);
        let latent_dim = crate::variable_supports::latent_dim(&h.problem);
        if latent_dim > 0 && !h.problem.anchors.variable_supports.is_empty() {
            let latents = crate::variable_supports::initial_latents(
                &h.problem.anchors,
                h.problem.solver.anchor_saturation_lambda,
            )?;
            h.state.variable_anchor_latents = latents;
            h.state.variable_anchor_positions = crate::variable_supports::map_latents_to_positions(
                &h.problem,
                &h.state.variable_anchor_latents,
            )?;
        }
        Ok(())
    }))
}

/// Configure q parameterization mode.
///
/// mode:
///   0 = DirectSoftBounds (default; optimize q directly + soft bounds)
///   2 = DirectBoxBounds  (optimize physical q with hard L-BFGS-B boxes)
///
/// # Safety
/// Valid handle.
#[no_mangle]
pub unsafe extern "C" fn theseus_set_q_parameterization_mode(
    handle: *mut TheseusHandle,
    mode: i32,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        h.problem.solver.q_parameterization_mode = QParameterizationMode::try_from(mode)?;
        Ok(())
    }))
}

// ─────────────────────────────────────────────────────────────
//  Self-weight configuration
// ─────────────────────────────────────────────────────────────

/// Configure prescribed-density self-weight.
///
/// `linear_densities` must point to `num_edges` doubles (mass/length per edge).
/// `gravity` must point to 3 doubles (gravity vector, e.g. [0, 0, -9.81]).
///
/// # Safety
/// Valid handle and buffers.
#[no_mangle]
pub unsafe extern "C" fn theseus_set_self_weight_prescribed(
    handle: *mut TheseusHandle,
    linear_densities: *const f64,
    gravity: *const f64,
    max_iters: usize,
    tolerance: f64,
    relaxation: f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        let ne = h.problem.topology.num_edges;
        let mu = slice::from_raw_parts(linear_densities, ne).to_vec();
        let g = slice::from_raw_parts(gravity, 3);
        h.problem.self_weight = Some(SelfWeightParams::Prescribed {
            linear_densities: mu,
            gravity: [g[0], g[1], g[2]],
            max_iters,
            tolerance,
            relaxation,
        });
        Ok(())
    }))
}

/// Configure force-based sizing self-weight.
///
/// Cross-section area is derived as `A_k = |F_k| / sigma`, giving
/// linear density `mu_k = rho * A_k`.
///
/// # Safety
/// Valid handle.
#[no_mangle]
pub unsafe extern "C" fn theseus_set_self_weight_sizing(
    handle: *mut TheseusHandle,
    rho: f64,
    sigma: f64,
    gravity: *const f64,
    max_iters: usize,
    tolerance: f64,
    relaxation: f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        let g = slice::from_raw_parts(gravity, 3);
        h.problem.self_weight = Some(SelfWeightParams::Sizing {
            rho,
            sigma,
            gravity: [g[0], g[1], g[2]],
            max_iters,
            tolerance,
            relaxation,
        });
        Ok(())
    }))
}

/// Clear self-weight configuration (revert to no self-weight).
///
/// # Safety
/// Valid handle.
#[no_mangle]
pub unsafe extern "C" fn theseus_clear_self_weight(handle: *mut TheseusHandle) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        h.problem.self_weight = None;
        Ok(())
    }))
}

// ─────────────────────────────────────────────────────────────
//  Pressure load configuration
// ─────────────────────────────────────────────────────────────

/// Configure pressure loads on faces.
///
/// `face_offsets` is a prefix-sum array of length `num_faces + 1`:
///   face `f` has vertices `face_vertices[face_offsets[f] .. face_offsets[f+1]]`.
/// `face_vertices` contains the ordered vertex indices (global node indices).
/// `pressures` has length `num_faces`.
///
/// # Safety
/// Valid handle and buffers.
#[no_mangle]
pub unsafe extern "C" fn theseus_set_pressure(
    handle: *mut TheseusHandle,
    num_faces: usize,
    face_offsets: *const usize,
    face_vertices: *const usize,
    pressures: *const f64,
    max_iters: usize,
    tolerance: f64,
    relaxation: f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        let offsets = slice::from_raw_parts(face_offsets, num_faces + 1);
        let total_verts = offsets[num_faces];
        let verts = slice::from_raw_parts(face_vertices, total_verts);
        let press = slice::from_raw_parts(pressures, num_faces);

        let mut faces = Vec::with_capacity(num_faces);
        for f in 0..num_faces {
            let start = offsets[f];
            let end = offsets[f + 1];
            faces.push(verts[start..end].to_vec());
        }

        h.problem.pressure = Some(PressureParams::Normal {
            face_topology: FaceTopology { faces },
            pressures: press.to_vec(),
            max_iters,
            tolerance,
            relaxation,
        });
        Ok(())
    }))
}

/// Configure hydrostatic pressure loads on faces.
///
/// Pressure varies linearly with depth below `z_datum`:
///   `p_f = rho_fluid * g_magnitude * max(0, z_datum - centroid_depth)`.
///
/// `up_direction` must point to 3 doubles (unit "up" vector, e.g. `[0,0,1]`).
///
/// # Safety
/// Valid handle and buffers.
#[no_mangle]
pub unsafe extern "C" fn theseus_set_pressure_hydrostatic(
    handle: *mut TheseusHandle,
    num_faces: usize,
    face_offsets: *const usize,
    face_vertices: *const usize,
    rho_fluid: f64,
    g_magnitude: f64,
    z_datum: f64,
    up_direction: *const f64,
    max_iters: usize,
    tolerance: f64,
    relaxation: f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        let offsets = slice::from_raw_parts(face_offsets, num_faces + 1);
        let total_verts = offsets[num_faces];
        let verts = slice::from_raw_parts(face_vertices, total_verts);
        let up = slice::from_raw_parts(up_direction, 3);

        let mut faces = Vec::with_capacity(num_faces);
        for f in 0..num_faces {
            faces.push(verts[offsets[f]..offsets[f + 1]].to_vec());
        }

        h.problem.pressure = Some(PressureParams::Hydrostatic {
            face_topology: FaceTopology { faces },
            rho_fluid,
            g_magnitude,
            z_datum,
            up_direction: [up[0], up[1], up[2]],
            max_iters,
            tolerance,
            relaxation,
        });
        Ok(())
    }))
}

/// Configure directional pressure loads on faces.
///
/// Load on each face is `p * max(0, n_f · d_hat) * d_hat`, proportional
/// to the face's projected area in the given direction.
///
/// `direction` must point to 3 doubles (unit load direction).
/// `pressures` has length `num_faces`.
///
/// # Safety
/// Valid handle and buffers.
#[no_mangle]
pub unsafe extern "C" fn theseus_set_pressure_directional(
    handle: *mut TheseusHandle,
    num_faces: usize,
    face_offsets: *const usize,
    face_vertices: *const usize,
    pressures: *const f64,
    direction: *const f64,
    max_iters: usize,
    tolerance: f64,
    relaxation: f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        let offsets = slice::from_raw_parts(face_offsets, num_faces + 1);
        let total_verts = offsets[num_faces];
        let verts = slice::from_raw_parts(face_vertices, total_verts);
        let press = slice::from_raw_parts(pressures, num_faces);
        let dir = slice::from_raw_parts(direction, 3);

        let mut faces = Vec::with_capacity(num_faces);
        for f in 0..num_faces {
            faces.push(verts[offsets[f]..offsets[f + 1]].to_vec());
        }

        h.problem.pressure = Some(PressureParams::Directional {
            face_topology: FaceTopology { faces },
            pressures: press.to_vec(),
            direction: [dir[0], dir[1], dir[2]],
            max_iters,
            tolerance,
            relaxation,
        });
        Ok(())
    }))
}

/// Clear pressure load configuration.
///
/// # Safety
/// Valid handle.
#[no_mangle]
pub unsafe extern "C" fn theseus_clear_pressure(handle: *mut TheseusHandle) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        h.problem.pressure = None;
        Ok(())
    }))
}

// ─────────────────────────────────────────────────────────────
//  Progress callback
// ─────────────────────────────────────────────────────────────

/// Register a progress callback invoked every `frequency` accepted L-BFGS iterations.
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
        let h = require_handle(handle)?;
        h.progress_callback = callback;
        h.report_frequency = if frequency == 0 { 1 } else { frequency };
        Ok(())
    }))
}

/// Retrieve the termination reason from the most recent optimization run.
///
/// Returns the number of bytes written (excluding null terminator), or -1 if
/// the handle is null or the buffer is too small.
///
/// # Safety
/// `buf` must point to at least `buf_len` writable bytes.
#[no_mangle]
pub unsafe extern "C" fn theseus_get_termination_reason(
    handle: *const TheseusHandle,
    buf: *mut u8,
    buf_len: usize,
) -> i32 {
    if handle.is_null() {
        set_last_error("null handle");
        return -1;
    }
    let msg = &(*handle).last_termination_reason;
    let bytes = msg.as_bytes();
    if buf_len < bytes.len() + 1 {
        return -1;
    }
    let out = slice::from_raw_parts_mut(buf, buf_len);
    out[..bytes.len()].copy_from_slice(bytes);
    out[bytes.len()] = 0;
    bytes.len() as i32
}

// ─────────────────────────────────────────────────────────────
//  Cancellation
// ─────────────────────────────────────────────────────────────

/// Request cancellation of an in-progress optimization.
///
/// The flag is checked after each objective evaluation. Safe to call
/// from any thread while `theseus_optimize` is running.
///
/// # Safety
/// `handle` must be a valid pointer returned by `theseus_create`.
#[no_mangle]
pub unsafe extern "C" fn theseus_cancel(handle: *mut TheseusHandle) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        h.cancel_flag.store(true, Ordering::Relaxed);
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
        let h = require_handle(handle)?;
        let cb = h.progress_callback;
        let freq = h.report_frequency;
        let result = optimizer::optimize(&h.problem, &mut h.state, cb, freq, &h.cancel_flag)?;
        h.last_termination_reason = result.termination_reason.clone();

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

/// Return the number of loss values recorded by the last optimization.
///
/// # Safety
/// `handle` must be a valid pointer returned by `theseus_create`.
#[no_mangle]
pub unsafe extern "C" fn theseus_get_loss_trace_len(handle: *mut TheseusHandle) -> usize {
    let Ok(h) = require_handle(handle) else {
        return 0;
    };
    h.state.loss_trace.len()
}

/// Copy the loss trace from the last optimization into `out_loss_trace`.
///
/// Returns the number of values copied, capped by `out_len`.
///
/// # Safety
/// `handle` must be valid and `out_loss_trace` must point to a buffer of at
/// least `out_len` doubles when `out_len > 0`.
#[no_mangle]
pub unsafe extern "C" fn theseus_get_loss_trace(
    handle: *mut TheseusHandle,
    out_loss_trace: *mut f64,
    out_len: usize,
) -> usize {
    let Ok(h) = require_handle(handle) else {
        return 0;
    };
    let n = h.state.loss_trace.len().min(out_len);
    if n == 0 {
        return 0;
    }
    slice::from_raw_parts_mut(out_loss_trace, n).copy_from_slice(&h.state.loss_trace[..n]);
    n
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
        let h = require_handle(handle)?;
        let mut cache = FdmCache::new(&h.problem)?;
        let anchors = crate::variable_supports::map_latents_to_positions(
            &h.problem,
            &h.state.variable_anchor_latents,
        )?;

        crate::fdm::solve_fdm_with_loads(
            &mut cache,
            &h.state.force_densities,
            &h.problem,
            &anchors,
            1e-12,
        )?;

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
//  Inverse solvers  (experimental)
// ─────────────────────────────────────────────────────────────

/// Solve for force densities via pseudoinverse of the equilibrium system.
///
/// Given target free-node positions, finds q that best satisfies M q = p.
/// When `use_l2` is non-zero, uses L2 (least-squares) via Tikhonov-regularised
/// normal equations (single solve).  When zero, uses L1 (sum of absolute
/// residuals) via IRLS for up to `max_l1_iter` iterations.
///
/// When `use_augmented` is non-zero, uses the augmented saddle-point system
/// instead of forming M^T M.  This avoids fill-in explosion and is faster
/// for large meshes (>50k edges).  Requires `regularization > 0`.
///
/// Then performs a forward FDM solve with the resulting q to produce final geometry.
///
/// `target_free_xyz` must point to `num_free * 3` doubles (row-major).
///
/// Returns 0 on success, -1 on error, -2 on internal panic.
///
/// # Safety
/// Valid handle and output buffers.
#[no_mangle]
pub unsafe extern "C" fn theseus_solve_pseudoinverse(
    handle: *mut TheseusHandle,
    target_free_xyz: *const f64,
    regularization: f64,
    use_l2: i32,
    max_l1_iter: usize,
    use_augmented: i32,
    enforce_zero_rx: i32,
    enforce_zero_ry: i32,
    enforce_zero_rz: i32,
    solve_for_q: i32,
    out_q: *mut f64,
    out_xyz: *mut f64,
    out_lengths: *mut f64,
    out_forces: *mut f64,
    out_reactions: *mut f64,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        let nn_free = h.problem.topology.free_node_indices.len();
        let nn = h.problem.topology.num_nodes;
        let ne = h.problem.topology.num_edges;

        let target = Array2::from_shape_vec(
            (nn_free, 3),
            slice::from_raw_parts(target_free_xyz, nn_free * 3).to_vec(),
        )
        .map_err(|e| TheseusError::Shape(format!("target_free_xyz: {e}")))?;

        let q = crate::inverse::solve_pseudoinverse_dispatch(
            &h.problem,
            &target,
            regularization,
            use_l2 != 0,
            max_l1_iter,
            use_augmented != 0,
            enforce_zero_rx != 0,
            enforce_zero_ry != 0,
            enforce_zero_rz != 0,
            solve_for_q != 0,
        )?;

        // Forward solve with the computed q (with self-weight/pressure if active)
        let mut cache = FdmCache::new(&h.problem)?;
        let anchors = crate::variable_supports::map_latents_to_positions(
            &h.problem,
            &h.state.variable_anchor_latents,
        )?;
        crate::fdm::solve_fdm_with_loads(&mut cache, &q, &h.problem, &anchors, 1e-12)?;

        // Copy outputs
        slice::from_raw_parts_mut(out_q, ne).copy_from_slice(&q);

        let xyz_out = slice::from_raw_parts_mut(out_xyz, nn * 3);
        for i in 0..nn {
            for d in 0..3 {
                xyz_out[i * 3 + d] = cache.nf[[i, d]];
            }
        }
        slice::from_raw_parts_mut(out_lengths, ne).copy_from_slice(&cache.member_lengths);
        slice::from_raw_parts_mut(out_forces, ne).copy_from_slice(&cache.member_forces);

        let r_out = slice::from_raw_parts_mut(out_reactions, nn * 3);
        for i in 0..nn {
            for d in 0..3 {
                r_out[i * 3 + d] = cache.reactions[[i, d]];
            }
        }

        Ok(())
    }))
}

/// Solve for non-negative force densities via NNLS (spectral projected gradient).
///
/// Given target free-node positions, finds q ≥ 0 that minimises ‖Mq − p‖².
/// Then performs a forward FDM solve with the resulting q.
///
/// `target_free_xyz` must point to `num_free * 3` doubles (row-major).
///
/// Returns 0 on success, -1 on error, -2 on internal panic.
///
/// # Safety
/// Valid handle and output buffers.
#[no_mangle]
pub unsafe extern "C" fn theseus_solve_nnls(
    handle: *mut TheseusHandle,
    target_free_xyz: *const f64,
    max_iter: usize,
    tol: f64,
    out_q: *mut f64,
    out_xyz: *mut f64,
    out_lengths: *mut f64,
    out_forces: *mut f64,
    out_reactions: *mut f64,
    out_iterations: *mut usize,
    out_converged: *mut bool,
) -> i32 {
    ffi_guard(AssertUnwindSafe(|| {
        let h = require_handle(handle)?;
        let nn_free = h.problem.topology.free_node_indices.len();
        let nn = h.problem.topology.num_nodes;
        let ne = h.problem.topology.num_edges;

        let target = Array2::from_shape_vec(
            (nn_free, 3),
            slice::from_raw_parts(target_free_xyz, nn_free * 3).to_vec(),
        )
        .map_err(|e| TheseusError::Shape(format!("target_free_xyz: {e}")))?;

        let result = crate::inverse::solve_nnls(&h.problem, &target, max_iter, tol)?;
        let q = result.q;

        // Forward solve with the computed q (with self-weight/pressure if active)
        let mut cache = FdmCache::new(&h.problem)?;
        let anchors = crate::variable_supports::map_latents_to_positions(
            &h.problem,
            &h.state.variable_anchor_latents,
        )?;
        crate::fdm::solve_fdm_with_loads(&mut cache, &q, &h.problem, &anchors, 1e-12)?;

        // Copy outputs
        slice::from_raw_parts_mut(out_q, ne).copy_from_slice(&q);

        let xyz_out = slice::from_raw_parts_mut(out_xyz, nn * 3);
        for i in 0..nn {
            for d in 0..3 {
                xyz_out[i * 3 + d] = cache.nf[[i, d]];
            }
        }
        slice::from_raw_parts_mut(out_lengths, ne).copy_from_slice(&cache.member_lengths);
        slice::from_raw_parts_mut(out_forces, ne).copy_from_slice(&cache.member_forces);

        let r_out = slice::from_raw_parts_mut(out_reactions, nn * 3);
        for i in 0..nn {
            for d in 0..3 {
                r_out[i * 3 + d] = cache.reactions[[i, d]];
            }
        }

        *out_iterations = result.iterations;
        *out_converged = result.converged;

        Ok(())
    }))
}
