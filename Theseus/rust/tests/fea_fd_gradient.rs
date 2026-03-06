//! Finite-difference gradient tests for the FEA solver.
//!
//! Builds a small 3D truss and verifies every component of the analytic
//! gradient against a central-difference estimate:
//!
//!     dJ/dθ_i  ≈  [ J(θ + h eᵢ) − J(θ − h eᵢ) ] / 2h
//!
//! Tests cover:
//!   - Node-position design variables (shape optimisation)
//!   - Cross-section area design variables (sizing optimisation)
//!   - Combined node-position + area variables
//!   - All six objective types: Compliance, MaxDisplacement,
//!     TargetDisplacement, MinWeight, MaxStress, TargetGeometry
//!   - Self-weight loads

use ndarray::Array2;
use theseus::fea_types::*;
use theseus::fea_objectives::*;
use theseus::fea_gradients::fea_value_and_gradient;

// ─────────────────────────────────────────────────────────────
//  Test truss geometry: 3D tower
// ─────────────────────────────────────────────────────────────
//
//  An 8-node 3D tower truss with 4 pinned supports (nodes 0-3)
//  at the base and 4 free nodes (4-7) at the top.  12 bar elements.
//
//  Base (z=0):   (0)---(1)      Top (z=2):   (4)---(5)
//                 |  X  |                      |  X  |
//                (2)---(3)                    (6)---(7)
//
//  Verticals: 0-4, 1-5, 2-6, 3-7
//  Base:      0-1, 0-2, 1-3, 2-3
//  Top:       4-5, 4-6, 5-7, 6-7
//
//  Supports: nodes 0,1,2,3 pinned (all DOFs fixed)

fn make_node_positions() -> Array2<f64> {
    Array2::from_shape_vec(
        (8, 3),
        vec![
            0.0, 0.0, 0.0,  // node 0 (support)
            2.0, 0.0, 0.0,  // node 1 (support)
            0.0, 2.0, 0.0,  // node 2 (support)
            2.0, 2.0, 0.0,  // node 3 (support)
            0.0, 0.0, 2.0,  // node 4 (free)
            2.0, 0.0, 2.0,  // node 5 (free)
            0.0, 2.0, 2.0,  // node 6 (free)
            2.0, 2.0, 2.0,  // node 7 (free)
        ],
    )
    .unwrap()
}

fn make_edge_nodes() -> Vec<(usize, usize)> {
    vec![
        // Verticals (0-3)
        (0, 4), (1, 5), (2, 6), (3, 7),
        // Top horizontals (4-7)
        (4, 5), (4, 6), (5, 7), (6, 7),
        // Diagonal braces (8-11): each face gets one diagonal
        (0, 5), (1, 7), (2, 4), (3, 6),
    ]
}

fn make_fea_problem(
    objectives: Vec<Box<dyn FeaObjectiveTrait>>,
    include_self_weight: bool,
) -> FeaProblem {
    let num_nodes = 8;
    let num_elements = 12;
    let node_positions = make_node_positions();
    let edge_nodes = make_edge_nodes();

    let materials = vec![FeaMaterial::default_steel()];
    let sections = vec![
        FeaSection { area: 0.005 },
        FeaSection { area: 0.008 },
    ];
    // Verticals (0..4) use section 0, rest use section 1
    let element_props: Vec<FeaElementProps> = (0..num_elements)
        .map(|e| FeaElementProps {
            material_idx: 0,
            section_idx: if e < 4 { 0 } else { 1 },
        })
        .collect();

    let supports = vec![
        FeaSupport::pinned(0),
        FeaSupport::pinned(1),
        FeaSupport::pinned(2),
        FeaSupport::pinned(3),
    ];

    let loads = vec![
        FeaLoad { node_idx: 4, force: [1000.0, 0.0, -10000.0] },
        FeaLoad { node_idx: 5, force: [1000.0, 0.0, -10000.0] },
        FeaLoad { node_idx: 6, force: [0.0, 1000.0, -10000.0] },
        FeaLoad { node_idx: 7, force: [0.0, 1000.0, -10000.0] },
    ];

    let dof_map = DofMap::from_supports(num_nodes, &supports);

    FeaProblem {
        num_nodes,
        num_elements,
        materials,
        sections,
        element_props,
        supports,
        loads,
        node_positions,
        edge_nodes,
        dof_map,
        include_self_weight,
        gravity: [0.0, 0.0, -9.81],
        objectives,
        solver: FeaSolverOptions {
            barrier_weight: 0.0,
            barrier_sharpness: 10.0,
            ..Default::default()
        },
    }
}

// ─────────────────────────────────────────────────────────────
//  FD test driver
// ─────────────────────────────────────────────────────────────

fn eval_fea_loss(
    problem: &FeaProblem,
    theta: &[f64],
    mask: &FeaVariableMask,
    lb: &[f64],
    ub: &[f64],
    lb_idx: &[usize],
    ub_idx: &[usize],
) -> f64 {
    let mut cache = FeaCache::new(problem).unwrap();
    let mut grad = vec![0.0; theta.len()];
    fea_value_and_gradient(
        &mut cache, problem, theta, &mut grad, mask, lb, ub, lb_idx, ub_idx,
    )
    .unwrap()
}

fn fea_fd_gradient_check(
    problem: &FeaProblem,
    theta: &[f64],
    mask: &FeaVariableMask,
    h: f64,
    tol_abs: f64,
    tol_rel: f64,
    label: &str,
) {
    let n = theta.len();
    let lb = vec![f64::NEG_INFINITY; n];
    let ub = vec![f64::INFINITY; n];
    let lb_idx: Vec<usize> = vec![];
    let ub_idx: Vec<usize> = vec![];

    // Analytic gradient
    let mut cache = FeaCache::new(problem).unwrap();
    let mut grad_analytic = vec![0.0; n];
    let loss = fea_value_and_gradient(
        &mut cache, problem, theta, &mut grad_analytic, mask, &lb, &ub, &lb_idx, &ub_idx,
    )
    .unwrap();

    // FD gradient
    let mut grad_fd = vec![0.0; n];
    let mut theta_p = theta.to_vec();
    let mut theta_m = theta.to_vec();

    for i in 0..n {
        theta_p[i] = theta[i] + h;
        theta_m[i] = theta[i] - h;

        let f_p = eval_fea_loss(problem, &theta_p, mask, &lb, &ub, &lb_idx, &ub_idx);
        let f_m = eval_fea_loss(problem, &theta_m, mask, &lb, &ub, &lb_idx, &ub_idx);
        grad_fd[i] = (f_p - f_m) / (2.0 * h);

        theta_p[i] = theta[i];
        theta_m[i] = theta[i];
    }

    // Diagnostics
    eprintln!("══════════════════════════════════════════════");
    eprintln!("FEA FD gradient check: {label}  (h = {h:.1e}, loss = {loss:.6e})");
    let mut max_abs = 0.0_f64;
    let mut max_rel = 0.0_f64;
    let mut worst_i = 0;
    for i in 0..n {
        let abs_err = (grad_analytic[i] - grad_fd[i]).abs();
        let denom = grad_fd[i].abs().max(grad_analytic[i].abs()).max(1e-14);
        let rel_err = abs_err / denom;
        if abs_err > max_abs {
            max_abs = abs_err;
            worst_i = i;
        }
        max_rel = max_rel.max(rel_err);
        let flag = if abs_err > tol_abs && rel_err > tol_rel { " <<<" } else { "" };
        eprintln!(
            "  θ[{i:>3}]  analytic={:+13.6e}  fd={:+13.6e}  abs={:.2e}  rel={:.2e}{flag}",
            grad_analytic[i], grad_fd[i], abs_err, rel_err,
        );
    }
    eprintln!("  max |g_a - g_fd| = {max_abs:.3e} at θ[{worst_i}]");
    eprintln!("  max relative err = {max_rel:.3e}");
    eprintln!("══════════════════════════════════════════════");

    // Assert
    for i in 0..n {
        let abs_err = (grad_analytic[i] - grad_fd[i]).abs();
        let denom = grad_fd[i].abs().max(grad_analytic[i].abs()).max(1e-14);
        let rel_err = abs_err / denom;
        assert!(
            abs_err < tol_abs || rel_err < tol_rel,
            "[{label}] θ[{i}]: analytic={:.8e}, fd={:.8e}, abs={:.3e}, rel={:.3e}",
            grad_analytic[i], grad_fd[i], abs_err, rel_err,
        );
    }
}

// ─────────────────────────────────────────────────────────────
//  Helper: pack initial theta from problem + mask
// ─────────────────────────────────────────────────────────────

fn pack_theta(problem: &FeaProblem, mask: &FeaVariableMask) -> Vec<f64> {
    let mut theta = Vec::new();
    if mask.node_positions {
        for &ni in &mask.node_indices {
            for d in 0..3 {
                theta.push(problem.node_positions[[ni, d]]);
            }
        }
    }
    if mask.cross_section_areas {
        for &si in &mask.area_element_indices {
            theta.push(problem.sections[si].area);
        }
    }
    if mask.support_positions {
        for &ni in &mask.support_node_indices {
            for d in 0..3 {
                theta.push(problem.node_positions[[ni, d]]);
            }
        }
    }
    theta
}

// ─────────────────────────────────────────────────────────────
//  Forward solve sanity check
// ─────────────────────────────────────────────────────────────

#[test]
fn fea_forward_solve_basic() {
    let problem = make_fea_problem(vec![], false);
    let mut cache = FeaCache::new(&problem).unwrap();
    let areas: Vec<f64> = problem.sections.iter().map(|s| s.area).collect();
    let result = theseus::fea_solve::solve_fea(
        &mut cache, &problem, &problem.node_positions, &areas,
    )
    .unwrap();

    // Supported nodes should have zero displacement
    for n in 0..4 {
        for d in 0..3 {
            assert!(result.displacements[n * 3 + d].abs() < 1e-12,
                "Support node {n} DOF {d} should have zero displacement, got {}",
                result.displacements[n * 3 + d]);
        }
    }

    // Loaded nodes should deflect downward (negative z)
    for n in 4..8 {
        assert!(result.displacements[n * 3 + 2] < 0.0,
            "Node {n} should deflect downward, got {}", result.displacements[n * 3 + 2]);
    }

    // All utilizations should be finite and non-negative
    for (e, &u) in result.utilization.iter().enumerate() {
        assert!(u.is_finite(), "Utilization at element {e} is not finite");
        assert!(u >= 0.0, "Utilization at element {e} is negative");
    }

    eprintln!("Forward solve: max displacement z = {:.6e}",
        (4..8).map(|n| result.displacements[n * 3 + 2]).fold(f64::INFINITY, f64::min));
}

// ─────────────────────────────────────────────────────────────
//  Compliance — node positions
// ─────────────────────────────────────────────────────────────

#[test]
fn fea_fd_compliance_node_pos() {
    let objectives: Vec<Box<dyn FeaObjectiveTrait>> = vec![
        Box::new(FeaCompliance { weight: 1.0 }),
    ];
    let problem = make_fea_problem(objectives, false);

    let mask = FeaVariableMask {
        node_positions: true,
        node_indices: vec![4, 5, 6, 7],
        ..Default::default()
    };

    let mut theta = pack_theta(&problem, &mask);
    theta[0] += 0.05;  // node 4 x
    theta[4] -= 0.03;  // node 5 y
    theta[8] += 0.02;  // node 6 z

    fea_fd_gradient_check(&problem, &theta, &mask, 1e-7, 1e-2, 1e-3, "Compliance / node pos");
}

// ─────────────────────────────────────────────────────────────
//  Compliance — cross-section areas
// ─────────────────────────────────────────────────────────────

#[test]
fn fea_fd_compliance_areas() {
    let objectives: Vec<Box<dyn FeaObjectiveTrait>> = vec![
        Box::new(FeaCompliance { weight: 1.0 }),
    ];
    let problem = make_fea_problem(objectives, false);

    let mask = FeaVariableMask {
        node_positions: false,
        cross_section_areas: true,
        area_element_indices: vec![0, 1],
        ..Default::default()
    };

    let theta = pack_theta(&problem, &mask);
    fea_fd_gradient_check(&problem, &theta, &mask, 1e-7, 1e-2, 1e-3, "Compliance / areas");
}

// ─────────────────────────────────────────────────────────────
//  MaxDisplacement — node positions
// ─────────────────────────────────────────────────────────────

#[test]
fn fea_fd_max_displacement_node_pos() {
    let objectives: Vec<Box<dyn FeaObjectiveTrait>> = vec![
        Box::new(FeaMaxDisplacement {
            weight: 1.0,
            node_indices: vec![4, 5, 6, 7],
        }),
    ];
    let problem = make_fea_problem(objectives, false);

    let mask = FeaVariableMask {
        node_positions: true,
        node_indices: vec![4, 5, 6, 7],
        ..Default::default()
    };

    let mut theta = pack_theta(&problem, &mask);
    theta[1] += 0.03;
    theta[7] -= 0.02;

    fea_fd_gradient_check(&problem, &theta, &mask, 1e-7, 1e-2, 1e-3, "MaxDisplacement / node pos");
}

// ─────────────────────────────────────────────────────────────
//  TargetDisplacement — node positions
// ─────────────────────────────────────────────────────────────

#[test]
fn fea_fd_target_displacement_node_pos() {
    let objectives: Vec<Box<dyn FeaObjectiveTrait>> = vec![
        Box::new(FeaTargetDisplacement {
            weight: 1.0,
            node_indices: vec![4, 5, 6, 7],
            targets: vec![
                [0.0, 0.0, -0.001],
                [0.0, 0.0, -0.001],
                [0.0, 0.0, -0.001],
                [0.0, 0.0, -0.001],
            ],
        }),
    ];
    let problem = make_fea_problem(objectives, false);

    let mask = FeaVariableMask {
        node_positions: true,
        node_indices: vec![4, 5, 6, 7],
        ..Default::default()
    };

    let mut theta = pack_theta(&problem, &mask);
    theta[2] += 0.04;
    theta[9] -= 0.03;

    fea_fd_gradient_check(&problem, &theta, &mask, 1e-7, 1e-2, 1e-3, "TargetDisplacement / node pos");
}

// ─────────────────────────────────────────────────────────────
//  MinWeight — node positions
// ─────────────────────────────────────────────────────────────

#[test]
fn fea_fd_min_weight_node_pos() {
    let objectives: Vec<Box<dyn FeaObjectiveTrait>> = vec![
        Box::new(FeaMinWeight { weight: 1.0 }),
    ];
    let problem = make_fea_problem(objectives, false);

    let mask = FeaVariableMask {
        node_positions: true,
        node_indices: vec![4, 5, 6, 7],
        ..Default::default()
    };

    let mut theta = pack_theta(&problem, &mask);
    theta[0] += 0.05;
    theta[5] -= 0.04;

    fea_fd_gradient_check(&problem, &theta, &mask, 1e-7, 1e-2, 1e-3, "MinWeight / node pos");
}

// ─────────────────────────────────────────────────────────────
//  MinWeight — cross-section areas
// ─────────────────────────────────────────────────────────────

#[test]
fn fea_fd_min_weight_areas() {
    let objectives: Vec<Box<dyn FeaObjectiveTrait>> = vec![
        Box::new(FeaMinWeight { weight: 1.0 }),
    ];
    let problem = make_fea_problem(objectives, false);

    let mask = FeaVariableMask {
        node_positions: false,
        cross_section_areas: true,
        area_element_indices: vec![0, 1],
        ..Default::default()
    };

    let theta = pack_theta(&problem, &mask);
    fea_fd_gradient_check(&problem, &theta, &mask, 1e-7, 1e-2, 1e-3, "MinWeight / areas");
}

// ─────────────────────────────────────────────────────────────
//  MaxStress — node positions
// ─────────────────────────────────────────────────────────────

#[test]
fn fea_fd_max_stress_node_pos() {
    let ne = 12;
    let objectives: Vec<Box<dyn FeaObjectiveTrait>> = vec![
        Box::new(FeaMaxStress {
            weight: 1.0,
            edge_indices: (0..ne).collect(),
            threshold: vec![100e6; ne],
            sharpness: 10.0,
        }),
    ];
    let problem = make_fea_problem(objectives, false);

    let mask = FeaVariableMask {
        node_positions: true,
        node_indices: vec![4, 5, 6, 7],
        ..Default::default()
    };

    let mut theta = pack_theta(&problem, &mask);
    theta[0] += 0.05;
    theta[3] -= 0.03;

    fea_fd_gradient_check(&problem, &theta, &mask, 1e-7, 1e-2, 1e-3, "MaxStress / node pos");
}

// ─────────────────────────────────────────────────────────────
//  TargetGeometry — node positions
// ─────────────────────────────────────────────────────────────

#[test]
fn fea_fd_target_geometry_node_pos() {
    let objectives: Vec<Box<dyn FeaObjectiveTrait>> = vec![
        Box::new(FeaTargetGeometry {
            weight: 1.0,
            node_indices: vec![4, 5, 6, 7],
            targets: vec![
                [0.1, 0.0, 1.95],
                [1.9, 0.0, 1.95],
                [0.1, 2.0, 1.95],
                [1.9, 2.0, 1.95],
            ],
        }),
    ];
    let problem = make_fea_problem(objectives, false);

    let mask = FeaVariableMask {
        node_positions: true,
        node_indices: vec![4, 5, 6, 7],
        ..Default::default()
    };

    let mut theta = pack_theta(&problem, &mask);
    theta[1] += 0.03;
    theta[10] -= 0.02;

    fea_fd_gradient_check(&problem, &theta, &mask, 1e-7, 1e-2, 1e-3, "TargetGeometry / node pos");
}

// ─────────────────────────────────────────────────────────────
//  Combined objectives — node positions + areas
// ─────────────────────────────────────────────────────────────

#[test]
fn fea_fd_combined_pos_and_areas() {
    let ne = 12;
    let objectives: Vec<Box<dyn FeaObjectiveTrait>> = vec![
        Box::new(FeaCompliance { weight: 1.0 }),
        Box::new(FeaMinWeight { weight: 0.001 }),
        Box::new(FeaMaxStress {
            weight: 0.1,
            edge_indices: (0..ne).collect(),
            threshold: vec![200e6; ne],
            sharpness: 10.0,
        }),
    ];
    let problem = make_fea_problem(objectives, false);

    let mask = FeaVariableMask {
        node_positions: true,
        node_indices: vec![4, 5, 6, 7],
        cross_section_areas: true,
        area_element_indices: vec![0, 1],
        ..Default::default()
    };

    let mut theta = pack_theta(&problem, &mask);
    theta[0] += 0.04;
    theta[5] -= 0.03;
    theta[8] += 0.02;

    fea_fd_gradient_check(&problem, &theta, &mask, 1e-7, 5e-2, 5e-3, "Combined / pos + areas");
}

// ─────────────────────────────────────────────────────────────
//  Self-weight — compliance with gravity
// ─────────────────────────────────────────────────────────────

#[test]
fn fea_fd_compliance_self_weight() {
    let objectives: Vec<Box<dyn FeaObjectiveTrait>> = vec![
        Box::new(FeaCompliance { weight: 1.0 }),
    ];
    let problem = make_fea_problem(objectives, true);

    let mask = FeaVariableMask {
        node_positions: true,
        node_indices: vec![4, 5, 6, 7],
        ..Default::default()
    };

    let mut theta = pack_theta(&problem, &mask);
    theta[2] += 0.03;
    theta[7] -= 0.02;

    fea_fd_gradient_check(&problem, &theta, &mask, 1e-7, 1e-1, 1e-2, "Compliance + self-weight / node pos");
}

// ─────────────────────────────────────────────────────────────
//  Self-weight — min weight with gravity
// ─────────────────────────────────────────────────────────────

#[test]
fn fea_fd_min_weight_self_weight() {
    let objectives: Vec<Box<dyn FeaObjectiveTrait>> = vec![
        Box::new(FeaMinWeight { weight: 1.0 }),
        Box::new(FeaCompliance { weight: 0.001 }),
    ];
    let problem = make_fea_problem(objectives, true);

    let mask = FeaVariableMask {
        node_positions: true,
        node_indices: vec![4, 5, 6, 7],
        cross_section_areas: true,
        area_element_indices: vec![0, 1],
        ..Default::default()
    };

    let mut theta = pack_theta(&problem, &mask);
    theta[0] += 0.03;
    theta[6] -= 0.02;

    fea_fd_gradient_check(&problem, &theta, &mask, 1e-7, 1e-1, 1e-2, "MinWeight + Compliance + self-weight");
}

// ─────────────────────────────────────────────────────────────
//  Kitchen sink: all objectives combined
// ─────────────────────────────────────────────────────────────

#[test]
fn fea_fd_all_objectives_combined() {
    let ne = 12;
    let objectives: Vec<Box<dyn FeaObjectiveTrait>> = vec![
        Box::new(FeaCompliance { weight: 1.0 }),
        Box::new(FeaMaxDisplacement {
            weight: 0.5,
            node_indices: vec![4, 5, 6, 7],
        }),
        Box::new(FeaTargetDisplacement {
            weight: 0.3,
            node_indices: vec![4, 5],
            targets: vec![
                [0.0, 0.0, -0.0005],
                [0.0, 0.0, -0.0005],
            ],
        }),
        Box::new(FeaMinWeight { weight: 0.001 }),
        Box::new(FeaMaxStress {
            weight: 0.1,
            edge_indices: (0..ne).collect(),
            threshold: vec![200e6; ne],
            sharpness: 10.0,
        }),
        Box::new(FeaTargetGeometry {
            weight: 0.2,
            node_indices: vec![4, 5],
            targets: vec![
                [0.05, 0.0, 1.98],
                [1.95, 0.0, 1.98],
            ],
        }),
    ];
    let problem = make_fea_problem(objectives, false);

    let mask = FeaVariableMask {
        node_positions: true,
        node_indices: vec![4, 5, 6, 7],
        cross_section_areas: true,
        area_element_indices: vec![0, 1],
        ..Default::default()
    };

    let mut theta = pack_theta(&problem, &mask);
    theta[0] += 0.04;
    theta[3] -= 0.02;
    theta[9] += 0.03;

    fea_fd_gradient_check(&problem, &theta, &mask, 1e-7, 5e-2, 5e-3, "All objectives combined");
}
