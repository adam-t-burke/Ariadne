//! FFI round-trip tests — call the `extern "C"` functions directly from Rust
//! to catch marshalling bugs before C# enters the picture.
//!
//! These tests mirror the safe-Rust integration tests in `integration.rs`
//! but go through the raw pointer / handle-based FFI boundary.

use std::ptr;

// Re-export the FFI functions from the crate (cdylib symbols).
// Because the crate also builds as `rlib`, we can link them directly.
use theseus::ffi::*;

// ─────────────────────────────────────────────────────────────
//  Shared arch-network data (same 7-node / 8-edge arch)
// ─────────────────────────────────────────────────────────────

struct ArchData {
    num_edges: usize,
    num_nodes: usize,
    num_free: usize,
    num_fixed: usize,
    coo_rows: Vec<usize>,
    coo_cols: Vec<usize>,
    coo_vals: Vec<f64>,
    free_idx: Vec<usize>,
    fixed_idx: Vec<usize>,
    loads: Vec<f64>,
    fixed_pos: Vec<f64>,
    q_init: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
}

fn arch_data() -> ArchData {
    let edges: Vec<(usize, usize)> = vec![
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6),
        (1, 5), (2, 4),
    ];
    let num_edges = edges.len();
    let num_nodes = 7;
    let num_free = 5;
    let num_fixed = 2;

    let mut coo_rows = Vec::new();
    let mut coo_cols = Vec::new();
    let mut coo_vals = Vec::new();
    for (e, &(s, t)) in edges.iter().enumerate() {
        coo_rows.push(e); coo_cols.push(s); coo_vals.push(-1.0);
        coo_rows.push(e); coo_cols.push(t); coo_vals.push(1.0);
    }

    let free_idx = vec![1, 2, 3, 4, 5];
    let fixed_idx = vec![0, 6];

    // Loads: 5 free nodes × 3 (unit downward)
    let loads = vec![
        0.0, 0.0, -1.0,
        0.0, 0.0, -1.0,
        0.0, 0.0, -2.0,
        0.0, 0.0, -1.0,
        0.0, 0.0, -1.0,
    ];

    // Fixed positions: node 0 at origin, node 6 at (6,0,0)
    let fixed_pos = vec![0.0, 0.0, 0.0, 6.0, 0.0, 0.0];

    let q_init = vec![1.0; num_edges];
    let lower = vec![0.1; num_edges];
    let upper = vec![100.0; num_edges];

    ArchData {
        num_edges, num_nodes, num_free, num_fixed,
        coo_rows, coo_cols, coo_vals,
        free_idx, fixed_idx,
        loads, fixed_pos, q_init, lower, upper,
    }
}

/// Create a handle via FFI; panics if null.
unsafe fn create_handle(d: &ArchData) -> *mut TheseusHandle {
    let h = theseus_create(
        d.num_edges, d.num_nodes, d.num_free,
        d.coo_rows.as_ptr(), d.coo_cols.as_ptr(), d.coo_vals.as_ptr(),
        d.coo_rows.len(),
        d.free_idx.as_ptr(), d.fixed_idx.as_ptr(), d.num_fixed,
        d.loads.as_ptr(), d.fixed_pos.as_ptr(),
        d.q_init.as_ptr(), d.lower.as_ptr(), d.upper.as_ptr(),
    );
    assert!(!h.is_null(), "theseus_create returned null: {}", get_last_error());
    h
}

fn get_last_error() -> String {
    let mut buf = vec![0u8; 1024];
    let n = unsafe { theseus_last_error(buf.as_mut_ptr(), buf.len()) };
    if n <= 0 { return String::from("(no error)"); }
    String::from_utf8_lossy(&buf[..n as usize]).to_string()
}

// ─────────────────────────────────────────────────────────────
//  Test: create / free round-trip
// ─────────────────────────────────────────────────────────────

#[test]
fn ffi_create_and_free() {
    let d = arch_data();
    unsafe {
        let h = create_handle(&d);
        theseus_free(h);
    }
    // Also verify freeing null is safe
    unsafe { theseus_free(ptr::null_mut()); }
}

// ─────────────────────────────────────────────────────────────
//  Test: forward solve through FFI
// ─────────────────────────────────────────────────────────────

#[test]
fn ffi_forward_solve() {
    let d = arch_data();
    unsafe {
        let h = create_handle(&d);

        let mut xyz = vec![0.0; d.num_nodes * 3];
        let mut lengths = vec![0.0; d.num_edges];
        let mut forces = vec![0.0; d.num_edges];
        let mut q_out = vec![0.0; d.num_edges];
        let mut reactions = vec![0.0; d.num_nodes * 3];

        let rc = theseus_solve_forward(
            h,
            xyz.as_mut_ptr(),
            lengths.as_mut_ptr(),
            forces.as_mut_ptr(),
            q_out.as_mut_ptr(),
            reactions.as_mut_ptr(),
        );
        assert_eq!(rc, 0, "forward solve failed: {}", get_last_error());

        // Anchors preserved
        assert!((xyz[0 * 3] - 0.0).abs() < 1e-12, "anchor 0 x");
        assert!((xyz[6 * 3] - 6.0).abs() < 1e-12, "anchor 6 x");

        // All positions finite
        for v in &xyz { assert!(v.is_finite(), "non-finite xyz: {v}"); }
        // All lengths positive
        for &l in &lengths { assert!(l > 0.0 && l.is_finite(), "bad length: {l}"); }
        // q values should match input
        for &qi in &q_out { assert!((qi - 1.0).abs() < 1e-12, "q should match init: {qi}"); }
        // Reactions should be finite
        for &r in &reactions { assert!(r.is_finite(), "non-finite reaction: {r}"); }

        theseus_free(h);
    }
}

// ─────────────────────────────────────────────────────────────
//  Test: optimise with TargetXYZ through FFI
// ─────────────────────────────────────────────────────────────

#[test]
fn ffi_optimize_target_xyz() {
    let d = arch_data();
    unsafe {
        let h = create_handle(&d);

        // Target positions for 5 free nodes
        let target_indices: Vec<usize> = vec![1, 2, 3, 4, 5];
        let target_xyz: Vec<f64> = vec![
            1.0, 0.0, 1.0,
            2.0, 0.0, 2.0,
            3.0, 0.0, 2.5,
            4.0, 0.0, 2.0,
            5.0, 0.0, 1.0,
        ];

        let rc = theseus_add_target_xyz(
            h, 1.0,
            target_indices.as_ptr(), target_indices.len(),
            target_xyz.as_ptr(),
        );
        assert_eq!(rc, 0, "add_target_xyz failed: {}", get_last_error());

        let rc = theseus_set_solver_options(h, 200, 1e-6, 1e-6, 1000.0, 10.0);
        assert_eq!(rc, 0, "set_solver_options failed: {}", get_last_error());

        let mut xyz = vec![0.0; d.num_nodes * 3];
        let mut lengths = vec![0.0; d.num_edges];
        let mut forces = vec![0.0; d.num_edges];
        let mut q = vec![0.0; d.num_edges];
        let mut reactions = vec![0.0; d.num_nodes * 3];
        let mut iterations: usize = 0;
        let mut converged: bool = false;

        let rc = theseus_optimize(
            h,
            xyz.as_mut_ptr(),
            lengths.as_mut_ptr(),
            forces.as_mut_ptr(),
            q.as_mut_ptr(),
            reactions.as_mut_ptr(),
            &mut iterations as *mut usize,
            &mut converged as *mut bool,
        );
        assert_eq!(rc, 0, "optimize failed: {}", get_last_error());
        assert!(iterations > 0, "should run at least 1 iteration");

        // All geometry finite and positive
        for &l in &lengths {
            assert!(l.is_finite() && l > 0.0, "bad length: {l}");
        }
        for &qi in &q {
            assert!(qi.is_finite() && qi > 0.0, "bad q: {qi}");
        }

        // Positions should move toward target
        let free_idx = &[1usize, 2, 3, 4, 5];
        let mut total_error = 0.0;
        for (i, &node) in free_idx.iter().enumerate() {
            for d in 0..3 {
                let diff = xyz[node * 3 + d] - target_xyz[i * 3 + d];
                total_error += diff * diff;
            }
        }
        assert!(total_error < 20.0, "total squared error = {total_error:.4}");

        theseus_free(h);
    }
}

// ─────────────────────────────────────────────────────────────
//  Test: combined objectives through FFI
// ─────────────────────────────────────────────────────────────

#[test]
fn ffi_optimize_combined() {
    let d = arch_data();
    unsafe {
        let h = create_handle(&d);

        // TargetXYZ
        let indices: Vec<usize> = vec![1, 2, 3, 4, 5];
        let target: Vec<f64> = vec![
            1.0, 0.0, 0.8,
            2.0, 0.0, 1.5,
            3.0, 0.0, 2.0,
            4.0, 0.0, 1.5,
            5.0, 0.0, 0.8,
        ];
        assert_eq!(0, theseus_add_target_xyz(h, 1.0, indices.as_ptr(), indices.len(), target.as_ptr()));

        // LengthVariation
        let all_edges: Vec<usize> = (0..d.num_edges).collect();
        assert_eq!(0, theseus_add_length_variation(h, 0.5, all_edges.as_ptr(), all_edges.len(), 20.0));

        // SumForceLength
        assert_eq!(0, theseus_add_sum_force_length(h, 0.01, all_edges.as_ptr(), all_edges.len()));

        assert_eq!(0, theseus_set_solver_options(h, 200, 1e-6, 1e-6, 1000.0, 10.0));

        let mut xyz = vec![0.0; d.num_nodes * 3];
        let mut lengths = vec![0.0; d.num_edges];
        let mut forces = vec![0.0; d.num_edges];
        let mut q_out = vec![0.0; d.num_edges];
        let mut reactions = vec![0.0; d.num_nodes * 3];
        let mut iterations: usize = 0;
        let mut converged: bool = false;

        let rc = theseus_optimize(
            h,
            xyz.as_mut_ptr(), lengths.as_mut_ptr(), forces.as_mut_ptr(),
            q_out.as_mut_ptr(), reactions.as_mut_ptr(),
            &mut iterations, &mut converged,
        );
        assert_eq!(rc, 0, "optimize failed: {}", get_last_error());
        assert!(iterations > 0);

        for &l in &lengths {
            assert!(l.is_finite() && l > 0.0);
        }

        theseus_free(h);
    }
}

// ─────────────────────────────────────────────────────────────
//  Test: all objective registration functions accept valid input
// ─────────────────────────────────────────────────────────────

#[test]
fn ffi_register_all_objectives() {
    let d = arch_data();
    unsafe {
        let h = create_handle(&d);

        let node_idx: Vec<usize> = vec![1, 2, 3];
        let edge_idx: Vec<usize> = vec![0, 1, 2];
        let anchor_idx: Vec<usize> = vec![0, 1]; // indices into fixed_node_indices

        let target_3x3: Vec<f64> = vec![
            1.0, 0.0, 1.0,
            2.0, 0.0, 2.0,
            3.0, 0.0, 3.0,
        ];
        let thresholds_3 = vec![0.5; 3];
        let targets_3 = vec![1.0; 3];

        let dirs_2x3: Vec<f64> = vec![
            0.0, 0.0, 1.0,
            0.0, 0.0, 1.0,
        ];
        let mags_2 = vec![5.0; 2];

        // All objective types
        assert_eq!(0, theseus_add_target_xyz(h, 1.0, node_idx.as_ptr(), node_idx.len(), target_3x3.as_ptr()));
        assert_eq!(0, theseus_add_target_xy(h, 1.0, node_idx.as_ptr(), node_idx.len(), target_3x3.as_ptr()));
        assert_eq!(0, theseus_add_target_length(h, 1.0, edge_idx.as_ptr(), edge_idx.len(), targets_3.as_ptr()));
        assert_eq!(0, theseus_add_length_variation(h, 0.5, edge_idx.as_ptr(), edge_idx.len(), 20.0));
        assert_eq!(0, theseus_add_force_variation(h, 0.5, edge_idx.as_ptr(), edge_idx.len(), 20.0));
        assert_eq!(0, theseus_add_sum_force_length(h, 0.01, edge_idx.as_ptr(), edge_idx.len()));
        assert_eq!(0, theseus_add_min_length(h, 1.0, edge_idx.as_ptr(), edge_idx.len(), thresholds_3.as_ptr(), 10.0));
        assert_eq!(0, theseus_add_max_length(h, 1.0, edge_idx.as_ptr(), edge_idx.len(), thresholds_3.as_ptr(), 10.0));
        assert_eq!(0, theseus_add_min_force(h, 1.0, edge_idx.as_ptr(), edge_idx.len(), thresholds_3.as_ptr(), 10.0));
        assert_eq!(0, theseus_add_max_force(h, 1.0, edge_idx.as_ptr(), edge_idx.len(), thresholds_3.as_ptr(), 10.0));
        assert_eq!(0, theseus_add_rigid_set_compare(h, 1.0, node_idx.as_ptr(), node_idx.len(), target_3x3.as_ptr()));
        assert_eq!(0, theseus_add_reaction_direction(h, 1.0, anchor_idx.as_ptr(), anchor_idx.len(), dirs_2x3.as_ptr()));
        assert_eq!(0, theseus_add_reaction_direction_magnitude(h, 1.0, anchor_idx.as_ptr(), anchor_idx.len(), dirs_2x3.as_ptr(), mags_2.as_ptr()));

        theseus_free(h);
    }
}

// ─────────────────────────────────────────────────────────────
//  Test: error reporting works
// ─────────────────────────────────────────────────────────────

#[test]
fn ffi_error_reporting() {
    // theseus_last_error should return 0 for "no error" initially
    let mut buf = vec![0u8; 256];
    let n = unsafe { theseus_last_error(buf.as_mut_ptr(), buf.len()) };
    // n == 0 means no error recorded (or the error was cleared)
    assert!(n >= 0, "unexpected negative from theseus_last_error");
}
