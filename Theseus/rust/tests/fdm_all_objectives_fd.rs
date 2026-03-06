//! Finite-difference gradient tests for ALL 13 FDM objective types.
//!
//! Each objective is tested in isolation on the 7-node arch network,
//! verifying every component of the analytic gradient against a
//! central-difference estimate:
//!
//!     dJ/dθ_i  ≈  [ J(θ + h eᵢ) − J(θ − h eᵢ) ] / 2h
//!
//! Objectives tested:
//!   1.  TargetXYZ
//!   2.  TargetXY
//!   3.  TargetPlane
//!   4.  PlanarConstraintAlongDirection
//!   5.  TargetLength
//!   6.  LengthVariation
//!   7.  ForceVariation
//!   8.  SumForceLength
//!   9.  MinLength
//!  10.  MaxLength
//!  11.  MinForce
//!  12.  MaxForce
//!  13.  RigidSetCompare
//!  14.  ReactionDirection
//!  15.  ReactionDirectionMagnitude
//!  16.  Kitchen-sink combined (all objectives at once)

use ndarray::Array2;
use theseus::sparse::SparseColMatOwned;
use theseus::types::*;

// ─────────────────────────────────────────────────────────────
//  Helpers: reuse the 7-node arch from fd_gradient.rs
// ─────────────────────────────────────────────────────────────

fn build_incidence(edges: &[(usize, usize)], num_nodes: usize) -> SparseColMatOwned {
    let ne = edges.len();
    let mut rows = Vec::with_capacity(ne * 2);
    let mut cols = Vec::with_capacity(ne * 2);
    let mut vals = Vec::with_capacity(ne * 2);
    for (e, &(s, t)) in edges.iter().enumerate() {
        rows.push(e);
        cols.push(s);
        vals.push(-1.0);
        rows.push(e);
        cols.push(t);
        vals.push(1.0);
    }
    SparseColMatOwned::from_coo(ne, num_nodes, &rows, &cols, &vals).unwrap()
}

fn make_arch_problem(bounds: Bounds, objectives: Vec<Box<dyn ObjectiveTrait>>) -> Problem {
    let num_nodes = 7;
    let num_edges = 8;

    let edges = vec![
        (0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6),
        (1, 5), (2, 4),
    ];

    let free_idx: Vec<usize> = vec![1, 2, 3, 4, 5];
    let fixed_idx: Vec<usize> = vec![0, 6];

    let incidence = build_incidence(&edges, num_nodes);
    let free_inc = incidence.extract_columns(&free_idx);
    let fixed_inc = incidence.extract_columns(&fixed_idx);

    let topology = NetworkTopology {
        incidence,
        free_incidence: free_inc,
        fixed_incidence: fixed_inc,
        num_edges,
        num_nodes,
        free_node_indices: free_idx,
        fixed_node_indices: fixed_idx,
    };

    let nn_free = 5;
    let free_node_loads = Array2::from_shape_vec(
        (nn_free, 3),
        vec![
            0.0, 0.0, -1.0,
            0.0, 0.0, -1.0,
            0.0, 0.0, -2.0,
            0.0, 0.0, -1.0,
            0.0, 0.0, -1.0,
        ],
    ).unwrap();

    let fixed_node_positions = Array2::from_shape_vec(
        (2, 3),
        vec![0.0, 0.0, 0.0, 6.0, 0.0, 0.0],
    ).unwrap();

    let anchors = AnchorInfo::all_fixed(fixed_node_positions.clone());

    Problem {
        topology,
        free_node_loads,
        fixed_node_positions,
        anchors,
        objectives,
        bounds,
        solver: SolverOptions::default(),
    }
}

fn default_bounds() -> Bounds {
    Bounds { lower: vec![0.1; 8], upper: vec![100.0; 8] }
}

fn default_theta() -> Vec<f64> {
    vec![2.0, 3.0, 1.5, 2.5, 1.0, 3.5, 2.0, 1.8]
}

// ─────────────────────────────────────────────────────────────
//  FD test driver
// ─────────────────────────────────────────────────────────────

fn eval_loss(problem: &Problem, theta: &[f64], lb: &[f64], ub: &[f64], lb_idx: &[usize], ub_idx: &[usize]) -> f64 {
    let mut cache = FdmCache::new(problem).unwrap();
    let mut grad = vec![0.0; theta.len()];
    theseus::gradients::value_and_gradient(
        &mut cache, problem, theta, &mut grad, lb, ub, lb_idx, ub_idx,
    ).unwrap()
}

fn fd_check(
    problem: &Problem,
    theta: &[f64],
    h: f64,
    tol_abs: f64,
    tol_rel: f64,
    label: &str,
) {
    let ne = problem.topology.num_edges;
    let n = theta.len();

    let lb: Vec<f64> = problem.bounds.lower.iter()
        .chain(std::iter::repeat(&f64::NEG_INFINITY).take(n.saturating_sub(ne)))
        .take(n).copied().collect();
    let ub: Vec<f64> = problem.bounds.upper.iter()
        .chain(std::iter::repeat(&f64::INFINITY).take(n.saturating_sub(ne)))
        .take(n).copied().collect();
    let lb_idx: Vec<usize> = (0..n).filter(|&i| lb[i].is_finite()).collect();
    let ub_idx: Vec<usize> = (0..n).filter(|&i| ub[i].is_finite()).collect();

    let mut cache = FdmCache::new(problem).unwrap();
    let mut grad_a = vec![0.0; n];
    let loss = theseus::gradients::value_and_gradient(
        &mut cache, problem, theta, &mut grad_a, &lb, &ub, &lb_idx, &ub_idx,
    ).unwrap();

    let mut grad_fd = vec![0.0; n];
    let mut tp = theta.to_vec();
    let mut tm = theta.to_vec();
    for i in 0..n {
        tp[i] = theta[i] + h;
        tm[i] = theta[i] - h;
        let fp = eval_loss(problem, &tp, &lb, &ub, &lb_idx, &ub_idx);
        let fm = eval_loss(problem, &tm, &lb, &ub, &lb_idx, &ub_idx);
        grad_fd[i] = (fp - fm) / (2.0 * h);
        tp[i] = theta[i];
        tm[i] = theta[i];
    }

    eprintln!("══════════════════════════════════════════════");
    eprintln!("FDM FD check: {label}  (h={h:.1e}, loss={loss:.6e})");
    let mut max_abs = 0.0_f64;
    let mut max_rel = 0.0_f64;
    let mut worst = 0;
    for i in 0..n {
        let ae = (grad_a[i] - grad_fd[i]).abs();
        let den = grad_fd[i].abs().max(grad_a[i].abs()).max(1e-14);
        let re = ae / den;
        if ae > max_abs { max_abs = ae; worst = i; }
        max_rel = max_rel.max(re);
        let flag = if ae > tol_abs && re > tol_rel { " <<<" } else { "" };
        eprintln!("  θ[{i:>2}]  a={:+13.6e}  fd={:+13.6e}  abs={:.2e}  rel={:.2e}{flag}",
            grad_a[i], grad_fd[i], ae, re);
    }
    eprintln!("  max_abs={max_abs:.3e} at θ[{worst}], max_rel={max_rel:.3e}");
    eprintln!("══════════════════════════════════════════════");

    for i in 0..n {
        let ae = (grad_a[i] - grad_fd[i]).abs();
        let den = grad_fd[i].abs().max(grad_a[i].abs()).max(1e-14);
        let re = ae / den;
        assert!(ae < tol_abs || re < tol_rel,
            "[{label}] θ[{i}]: a={:.8e}, fd={:.8e}, abs={:.3e}, rel={:.3e}",
            grad_a[i], grad_fd[i], ae, re);
    }
}

// ─────────────────────────────────────────────────────────────
//  1. TargetXYZ
// ─────────────────────────────────────────────────────────────

#[test]
fn fdm_fd_target_xyz() {
    let target = Array2::from_shape_vec((5, 3), vec![
        1.0,0.0,1.0, 2.0,0.0,2.0, 3.0,0.0,2.5, 4.0,0.0,2.0, 5.0,0.0,1.0,
    ]).unwrap();
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![Box::new(TargetXYZ {
        weight: 1.0, node_indices: vec![1,2,3,4,5], target,
    })];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "TargetXYZ");
}

// ─────────────────────────────────────────────────────────────
//  2. TargetXY
// ─────────────────────────────────────────────────────────────

#[test]
fn fdm_fd_target_xy() {
    let target = Array2::from_shape_vec((5, 3), vec![
        1.0,0.0,0.0, 2.0,0.0,0.0, 3.0,0.0,0.0, 4.0,0.0,0.0, 5.0,0.0,0.0,
    ]).unwrap();
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![Box::new(TargetXY {
        weight: 1.0, node_indices: vec![1,2,3,4,5], target,
    })];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "TargetXY");
}

// ─────────────────────────────────────────────────────────────
//  3. TargetPlane
// ─────────────────────────────────────────────────────────────

#[test]
fn fdm_fd_target_plane() {
    let target = Array2::from_shape_vec((5, 3), vec![
        1.0,0.0,1.0, 2.0,0.0,2.0, 3.0,0.0,2.5, 4.0,0.0,2.0, 5.0,0.0,1.0,
    ]).unwrap();
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![Box::new(TargetPlane {
        weight: 1.0,
        node_indices: vec![1,2,3,4,5],
        target,
        origin: [0.0, 0.0, 0.0],
        x_axis: [1.0, 0.0, 0.0],
        y_axis: [0.0, 0.0, 1.0],
    })];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "TargetPlane");
}

// ─────────────────────────────────────────────────────────────
//  4. PlanarConstraintAlongDirection
// ─────────────────────────────────────────────────────────────

#[test]
fn fdm_fd_planar_constraint() {
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![Box::new(PlanarConstraintAlongDirection {
        weight: 1.0,
        node_indices: vec![1,2,3,4,5],
        origin: [3.0, 0.0, 0.0],
        x_axis: [1.0, 0.0, 0.0],
        y_axis: [0.0, 1.0, 0.0],
        direction: [0.0, 0.0, 1.0],
    })];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "PlanarConstraintAlongDirection");
}

// ─────────────────────────────────────────────────────────────
//  5. TargetLength
// ─────────────────────────────────────────────────────────────

#[test]
fn fdm_fd_target_length() {
    let ne = 8;
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![Box::new(TargetLength {
        weight: 1.0, edge_indices: (0..ne).collect(), target: vec![1.5; ne],
    })];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "TargetLength");
}

// ─────────────────────────────────────────────────────────────
//  6. LengthVariation
// ─────────────────────────────────────────────────────────────

#[test]
fn fdm_fd_length_variation() {
    let ne = 8;
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![Box::new(LengthVariation {
        weight: 1.0, edge_indices: (0..ne).collect(), sharpness: 20.0,
    })];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "LengthVariation");
}

// ─────────────────────────────────────────────────────────────
//  7. ForceVariation
// ─────────────────────────────────────────────────────────────

#[test]
fn fdm_fd_force_variation() {
    let ne = 8;
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![Box::new(ForceVariation {
        weight: 1.0, edge_indices: (0..ne).collect(), sharpness: 20.0,
    })];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "ForceVariation");
}

// ─────────────────────────────────────────────────────────────
//  8. SumForceLength
// ─────────────────────────────────────────────────────────────

#[test]
fn fdm_fd_sum_force_length() {
    let ne = 8;
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![Box::new(SumForceLength {
        weight: 1.0, edge_indices: (0..ne).collect(),
    })];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "SumForceLength");
}

// ─────────────────────────────────────────────────────────────
//  9. MinLength
// ─────────────────────────────────────────────────────────────

#[test]
fn fdm_fd_min_length() {
    let ne = 8;
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![Box::new(MinLength {
        weight: 1.0,
        edge_indices: (0..ne).collect(),
        threshold: vec![0.5; ne],
        sharpness: 10.0,
    })];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "MinLength");
}

// ─────────────────────────────────────────────────────────────
//  10. MaxLength
// ─────────────────────────────────────────────────────────────

#[test]
fn fdm_fd_max_length() {
    let ne = 8;
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![Box::new(MaxLength {
        weight: 1.0,
        edge_indices: (0..ne).collect(),
        threshold: vec![3.0; ne],
        sharpness: 10.0,
    })];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "MaxLength");
}

// ─────────────────────────────────────────────────────────────
//  11. MinForce
// ─────────────────────────────────────────────────────────────

#[test]
fn fdm_fd_min_force() {
    let ne = 8;
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![Box::new(MinForce {
        weight: 1.0,
        edge_indices: (0..ne).collect(),
        threshold: vec![0.5; ne],
        sharpness: 10.0,
    })];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "MinForce");
}

// ─────────────────────────────────────────────────────────────
//  12. MaxForce
// ─────────────────────────────────────────────────────────────

#[test]
fn fdm_fd_max_force() {
    let ne = 8;
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![Box::new(MaxForce {
        weight: 1.0,
        edge_indices: (0..ne).collect(),
        threshold: vec![5.0; ne],
        sharpness: 10.0,
    })];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "MaxForce");
}

// ─────────────────────────────────────────────────────────────
//  13. RigidSetCompare
// ─────────────────────────────────────────────────────────────

#[test]
fn fdm_fd_rigid_set_compare() {
    let target = Array2::from_shape_vec((3, 3), vec![
        1.0, 0.0, 1.0,
        3.0, 0.0, 2.5,
        5.0, 0.0, 1.0,
    ]).unwrap();
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![Box::new(RigidSetCompare {
        weight: 1.0, node_indices: vec![1, 3, 5], target,
    })];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "RigidSetCompare");
}

// ─────────────────────────────────────────────────────────────
//  14. ReactionDirection
// ─────────────────────────────────────────────────────────────

#[test]
fn fdm_fd_reaction_direction() {
    let dirs = Array2::from_shape_vec((2, 3), vec![
        0.0, 0.0, 1.0,
        0.0, 0.0, 1.0,
    ]).unwrap();
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![Box::new(ReactionDirection {
        weight: 1.0, anchor_indices: vec![0, 6], target_directions: dirs,
    })];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "ReactionDirection");
}

// ─────────────────────────────────────────────────────────────
//  15. ReactionDirectionMagnitude
// ─────────────────────────────────────────────────────────────

#[test]
fn fdm_fd_reaction_direction_magnitude() {
    let dirs = Array2::from_shape_vec((2, 3), vec![
        0.0, 0.0, 1.0,
        0.0, 0.0, 1.0,
    ]).unwrap();
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![Box::new(ReactionDirectionMagnitude {
        weight: 1.0,
        anchor_indices: vec![0, 6],
        target_directions: dirs,
        target_magnitudes: vec![3.0, 3.0],
    })];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "ReactionDirectionMagnitude");
}

// ─────────────────────────────────────────────────────────────
//  16. Kitchen sink: ALL objectives combined
// ─────────────────────────────────────────────────────────────

#[test]
fn fdm_fd_all_objectives_combined() {
    let ne = 8;
    let target_xyz = Array2::from_shape_vec((5, 3), vec![
        1.0,0.0,1.0, 2.0,0.0,2.0, 3.0,0.0,2.5, 4.0,0.0,2.0, 5.0,0.0,1.0,
    ]).unwrap();
    let target_xy = Array2::from_shape_vec((5, 3), vec![
        1.0,0.0,0.0, 2.0,0.0,0.0, 3.0,0.0,0.0, 4.0,0.0,0.0, 5.0,0.0,0.0,
    ]).unwrap();
    let target_plane = Array2::from_shape_vec((3, 3), vec![
        1.0,0.0,1.0, 3.0,0.0,2.5, 5.0,0.0,1.0,
    ]).unwrap();
    let rigid_target = Array2::from_shape_vec((3, 3), vec![
        1.0,0.0,1.0, 3.0,0.0,2.5, 5.0,0.0,1.0,
    ]).unwrap();
    let reaction_dirs = Array2::from_shape_vec((2, 3), vec![
        0.0,0.0,1.0, 0.0,0.0,1.0,
    ]).unwrap();
    let reaction_dirs2 = Array2::from_shape_vec((2, 3), vec![
        0.0,0.0,1.0, 0.0,0.0,1.0,
    ]).unwrap();

    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![
        Box::new(TargetXYZ { weight: 1.0, node_indices: vec![1,2,3,4,5], target: target_xyz }),
        Box::new(TargetXY { weight: 0.5, node_indices: vec![1,2,3,4,5], target: target_xy }),
        Box::new(TargetPlane {
            weight: 0.3, node_indices: vec![1,3,5], target: target_plane,
            origin: [0.0,0.0,0.0], x_axis: [1.0,0.0,0.0], y_axis: [0.0,0.0,1.0],
        }),
        Box::new(PlanarConstraintAlongDirection {
            weight: 0.2, node_indices: vec![1,2,3,4,5],
            origin: [3.0,0.0,0.0], x_axis: [1.0,0.0,0.0], y_axis: [0.0,1.0,0.0],
            direction: [0.0,0.0,1.0],
        }),
        Box::new(TargetLength { weight: 0.4, edge_indices: (0..ne).collect(), target: vec![1.5; ne] }),
        Box::new(LengthVariation { weight: 0.3, edge_indices: (0..ne).collect(), sharpness: 20.0 }),
        Box::new(ForceVariation { weight: 0.2, edge_indices: (0..ne).collect(), sharpness: 15.0 }),
        Box::new(SumForceLength { weight: 0.1, edge_indices: (0..ne).collect() }),
        Box::new(MinLength { weight: 0.3, edge_indices: (0..ne).collect(), threshold: vec![0.5; ne], sharpness: 10.0 }),
        Box::new(MaxLength { weight: 0.3, edge_indices: (0..ne).collect(), threshold: vec![3.0; ne], sharpness: 10.0 }),
        Box::new(MinForce { weight: 0.2, edge_indices: (0..ne).collect(), threshold: vec![0.5; ne], sharpness: 10.0 }),
        Box::new(MaxForce { weight: 0.2, edge_indices: (0..ne).collect(), threshold: vec![5.0; ne], sharpness: 10.0 }),
        Box::new(RigidSetCompare { weight: 0.5, node_indices: vec![1,3,5], target: rigid_target }),
        Box::new(ReactionDirection { weight: 0.3, anchor_indices: vec![0,6], target_directions: reaction_dirs }),
        Box::new(ReactionDirectionMagnitude {
            weight: 0.2, anchor_indices: vec![0,6], target_directions: reaction_dirs2,
            target_magnitudes: vec![3.0, 3.0],
        }),
    ];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "ALL 15 objectives combined");
}

// ─────────────────────────────────────────────────────────────
//  Subset-edge tests (verify partial edge index lists)
// ─────────────────────────────────────────────────────────────

#[test]
fn fdm_fd_target_length_subset() {
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![Box::new(TargetLength {
        weight: 1.0, edge_indices: vec![0, 2, 5, 7], target: vec![1.0, 1.5, 2.0, 1.2],
    })];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "TargetLength (subset)");
}

#[test]
fn fdm_fd_min_length_subset() {
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![Box::new(MinLength {
        weight: 1.0, edge_indices: vec![1, 3, 6], threshold: vec![0.3, 0.5, 0.4], sharpness: 10.0,
    })];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "MinLength (subset)");
}

#[test]
fn fdm_fd_max_force_subset() {
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![Box::new(MaxForce {
        weight: 1.0, edge_indices: vec![0, 4, 7], threshold: vec![4.0, 3.0, 5.0], sharpness: 10.0,
    })];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "MaxForce (subset)");
}

// ─────────────────────────────────────────────────────────────
//  Non-unit weights
// ─────────────────────────────────────────────────────────────

#[test]
fn fdm_fd_weighted_objectives() {
    let ne = 8;
    let target = Array2::from_shape_vec((5, 3), vec![
        1.0,0.0,1.0, 2.0,0.0,2.0, 3.0,0.0,2.5, 4.0,0.0,2.0, 5.0,0.0,1.0,
    ]).unwrap();
    let obj: Vec<Box<dyn ObjectiveTrait>> = vec![
        Box::new(TargetXYZ { weight: 3.7, node_indices: vec![1,2,3,4,5], target }),
        Box::new(SumForceLength { weight: 0.05, edge_indices: (0..ne).collect() }),
    ];
    let p = make_arch_problem(default_bounds(), obj);
    fd_check(&p, &default_theta(), 1e-6, 1e-4, 1e-3, "Weighted TargetXYZ + SumForceLength");
}
