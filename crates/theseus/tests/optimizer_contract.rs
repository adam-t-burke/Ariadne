//! Observable optimizer behavior shared by both parameterization modes.

use ndarray::Array2;
use std::cell::RefCell;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;
use theseus::optimizer::optimize;
use theseus::sparse::SparseColMatOwned;
use theseus::types::*;

type ProgressRecord = (usize, f64, Vec<f64>, Vec<f64>);

thread_local! {
    static PROGRESS: RefCell<Vec<ProgressRecord>> = const { RefCell::new(Vec::new()) };
}

unsafe extern "C" fn record(
    iteration: usize,
    loss: f64,
    xyz: *const f64,
    nodes: usize,
    q: *const f64,
    edges: usize,
) -> u8 {
    PROGRESS.with_borrow_mut(|records| {
        records.push((
            iteration,
            loss,
            unsafe { std::slice::from_raw_parts(xyz, 3 * nodes) }.to_vec(),
            unsafe { std::slice::from_raw_parts(q, edges) }.to_vec(),
        ));
    });
    1
}

unsafe extern "C" fn cancel(
    _: usize,
    _: f64,
    _: *const f64,
    _: usize,
    _: *const f64,
    _: usize,
) -> u8 {
    0
}

fn problem(mode: QParameterizationMode) -> Problem {
    let incidence = SparseColMatOwned::from_coo(1, 2, &[0, 0], &[0, 1], &[-1.0, 1.0]).unwrap();
    let positions = Array2::zeros((1, 3));
    Problem {
        topology: NetworkTopology {
            free_incidence: incidence.extract_columns(&[0]),
            fixed_incidence: incidence.extract_columns(&[1]),
            incidence,
            num_edges: 1,
            num_nodes: 2,
            free_node_indices: vec![0],
            fixed_node_indices: vec![1],
        },
        free_node_loads: Array2::from_shape_vec((1, 3), vec![0.0, 0.0, -1.0]).unwrap(),
        fixed_node_positions: positions.clone(),
        anchors: AnchorInfo::all_fixed(positions),
        objectives: vec![Box::new(TargetXYZ {
            weight: 1.0,
            node_indices: vec![0],
            target: Array2::from_shape_vec((1, 3), vec![0.0, 0.0, -0.25]).unwrap(),
            reduction: TargetGeometryReduction::Sse,
        })],
        bounds: Bounds {
            lower: vec![0.1],
            upper: vec![10.0],
        },
        solver: SolverOptions {
            q_parameterization_mode: mode,
            max_iterations: 100,
            absolute_tolerance: 1e-10,
            relative_tolerance: 0.0,
            barrier_weight: 0.0,
            ..SolverOptions::default()
        },
        self_weight: None,
        pressure: None,
    }
}

const MODES: [QParameterizationMode; 2] = [
    QParameterizationMode::DirectSoftBounds,
    QParameterizationMode::DirectBoxBounds,
];

#[test]
fn progress_contains_accepted_geometry_and_obeys_frequency() {
    for mode in MODES {
        for frequency in [0, 1, 3] {
            PROGRESS.with_borrow_mut(Vec::clear);
            let p = problem(mode);
            let mut state = OptimizationState::new(vec![1.0], Array2::zeros((0, 3)));
            let result = optimize(
                &p,
                &mut state,
                Some(record),
                frequency,
                &AtomicBool::new(false),
            )
            .unwrap();
            assert!(result.converged, "{mode:?}: {}", result.termination_reason);
            assert!((result.q[0] - 4.0).abs() < 1e-6);
            assert!((result.xyz[[0, 2]] + 1.0 / result.q[0]).abs() < 1e-12);
            assert_eq!(state.force_densities, result.q);
            assert_eq!(state.loss_trace, result.loss_trace);
            PROGRESS.with_borrow(|records| {
                let expected: Vec<_> = (1..=result.iterations)
                    .filter(|i| *i == 1 || i % frequency.max(1) == 0)
                    .collect();
                assert_eq!(records.iter().map(|r| r.0).collect::<Vec<_>>(), expected);
                for (_, loss, xyz, q) in records {
                    assert!((xyz[2] + 1.0 / q[0]).abs() < 1e-12);
                    assert!((loss - (xyz[2] + 0.25).powi(2)).abs() < 1e-12);
                }
                assert!(records.windows(2).all(|w| w[1].1 <= w[0].1));
            });
        }
    }
}

#[test]
fn cancellation_preserves_input_state_in_both_modes() {
    for mode in MODES {
        let p = problem(mode);
        for pre_cancelled in [false, true] {
            let mut state = OptimizationState::new(vec![1.0], Array2::zeros((0, 3)));
            let error = optimize(
                &p,
                &mut state,
                Some(cancel),
                1,
                &AtomicBool::new(pre_cancelled),
            )
            .unwrap_err();
            assert!(matches!(error, TheseusError::Cancelled), "{error}");
            assert_eq!(state.force_densities, vec![1.0]);
            assert_eq!(state.iterations, 0);
            assert!(state.loss_trace.is_empty());
        }
    }
}

#[test]
fn iteration_limit_is_not_convergence() {
    for mode in MODES {
        let mut p = problem(mode);
        p.solver.max_iterations = 1;
        p.solver.absolute_tolerance = 0.0;
        let mut state = OptimizationState::new(vec![1.0], Array2::zeros((0, 3)));
        let result = optimize(&p, &mut state, None, 1, &AtomicBool::new(false)).unwrap();
        assert_eq!(result.iterations, 1);
        assert!(!result.converged);
        assert!((result.xyz[[0, 2]] + 1.0 / result.q[0]).abs() < 1e-12);
    }
}

#[test]
fn soft_bounds_accept_one_sided_boxes() {
    let mut p = problem(QParameterizationMode::DirectSoftBounds);
    p.bounds.upper[0] = f64::INFINITY;
    let mut state = OptimizationState::new(vec![1.0], Array2::zeros((0, 3)));
    let result = optimize(&p, &mut state, None, 1, &AtomicBool::new(false)).unwrap();
    assert!((result.q[0] - 4.0).abs() < 1e-6);
}

#[test]
fn non_finite_objective_aborts_without_publishing_state() {
    for mode in MODES {
        let mut p = problem(mode);
        p.objectives = vec![Box::new(TargetXYZ {
            weight: 1.0,
            node_indices: vec![0],
            target: Array2::from_elem((1, 3), f64::NAN),
            reduction: TargetGeometryReduction::Sse,
        })];
        let mut state = OptimizationState::new(vec![1.0], Array2::zeros((0, 3)));
        assert!(optimize(&p, &mut state, None, 1, &AtomicBool::new(false)).is_err());
        assert_eq!(state.force_densities, vec![1.0]);
    }
}

#[derive(Debug)]
struct CancellingObjective(Arc<AtomicBool>);

impl ObjectiveTrait for CancellingObjective {
    fn loss(&self, _: &GeometrySnapshot) -> f64 {
        self.0.store(true, Ordering::Release);
        0.0
    }
    fn accumulate_gradient(&self, _: &mut FdmCache, _: &Problem) {}
    fn weight(&self) -> f64 {
        1.0
    }
}

#[test]
fn cancellation_during_evaluation_is_a_typed_error() {
    for mode in MODES {
        let flag = Arc::new(AtomicBool::new(false));
        let mut p = problem(mode);
        p.objectives = vec![Box::new(CancellingObjective(flag.clone()))];
        let mut state = OptimizationState::new(vec![1.0], Array2::zeros((0, 3)));
        assert!(matches!(
            optimize(&p, &mut state, None, 1, &flag),
            Err(TheseusError::Cancelled)
        ));
        assert_eq!(state.force_densities, vec![1.0]);
        assert!(state.loss_trace.is_empty());
    }
}

#[test]
fn stationary_initial_point_does_not_report_a_step() {
    for mode in MODES {
        PROGRESS.with_borrow_mut(Vec::clear);
        let p = problem(mode);
        let mut state = OptimizationState::new(vec![4.0], Array2::zeros((0, 3)));
        let result = optimize(&p, &mut state, Some(record), 1, &AtomicBool::new(false)).unwrap();
        assert!(result.converged);
        assert_eq!(result.iterations, 0);
        PROGRESS.with_borrow(|records| assert!(records.is_empty()));
    }
}

#[test]
fn variable_support_positions_and_scaled_latents_are_consistent() {
    for mode in MODES {
        for kind in [
            VariableSupportKind::Sphere { radius: 2.0 },
            VariableSupportKind::Rail {
                start: [-1.0, 0.0, 0.0],
                end: [1.0, 0.0, 0.0],
            },
        ] {
            let mut p = problem(mode);
            p.anchors.variable_indices = vec![1];
            p.anchors.fixed_indices = vec![1];
            p.anchors.initial_variable_positions = Array2::zeros((1, 3));
            p.anchors.variable_supports = vec![VariableSupport {
                node_index: 1,
                reference_position: [0.0; 3],
                saturation_lambda: 3.0,
                kind,
            }];
            p.objectives = vec![Box::new(TargetXYZ {
                weight: 1.0,
                node_indices: vec![0, 1],
                target: Array2::from_shape_vec((2, 3), vec![0.5, 0.0, -0.25, 0.5, 0.0, 0.0])
                    .unwrap(),
                reduction: TargetGeometryReduction::Sse,
            })];
            let mut state =
                OptimizationState::new(vec![1.0], p.anchors.initial_variable_positions.clone());
            state.variable_anchor_latents = theseus::variable_supports::initial_parameters(
                &p.anchors,
                p.solver.anchor_saturation_lambda,
                mode,
            )
            .unwrap();
            let result = optimize(&p, &mut state, None, 1, &AtomicBool::new(false)).unwrap();
            let mapped = theseus::variable_supports::map_latents_to_positions(
                &p,
                &state.variable_anchor_latents,
            )
            .unwrap();
            assert_eq!(mapped, result.anchor_positions);
            for (actual, expected) in result.anchor_positions.iter().zip([0.5, 0.0, 0.0]) {
                assert!((actual - expected).abs() < 1e-5, "{mode:?}: {result:?}");
            }
            assert!((result.q[0] - 4.0).abs() < 1e-4);
            assert!((result.xyz[[0, 2]] - result.xyz[[1, 2]] + 1.0 / result.q[0]).abs() < 1e-12);
        }
    }
}

#[derive(Debug)]
struct InconsistentGradient;

impl ObjectiveTrait for InconsistentGradient {
    fn loss(&self, _: &GeometrySnapshot) -> f64 {
        1.0
    }
    fn accumulate_gradient(&self, cache: &mut FdmCache, _: &Problem) {
        // No step can satisfy sufficient decrease for this deliberately invalid derivative.
        cache.grad_q[0] += 1.0;
    }
    fn weight(&self) -> f64 {
        1.0
    }
}

#[test]
fn failed_line_search_returns_consistent_geometry_without_convergence() {
    for mode in MODES {
        PROGRESS.with_borrow_mut(Vec::clear);
        let mut p = problem(mode);
        p.objectives = vec![Box::new(InconsistentGradient)];
        let mut state = OptimizationState::new(vec![1.0], Array2::zeros((0, 3)));
        let result = optimize(&p, &mut state, Some(record), 1, &AtomicBool::new(false)).unwrap();
        assert!(!result.converged);
        assert!(
            result.termination_reason.starts_with("failed:"),
            "{}",
            result.termination_reason
        );
        assert_eq!(result.iterations, 0);
        assert_eq!(result.q, vec![1.0]);
        let mut cache = FdmCache::new(&p).unwrap();
        theseus::fdm::solve_fdm_with_loads(
            &mut cache,
            &result.q,
            &p,
            &result.anchor_positions,
            1e-12,
        )
        .unwrap();
        assert_eq!(result.xyz, cache.nf);
        assert_eq!(result.member_forces, cache.member_forces);
        PROGRESS.with_borrow(|records| assert!(records.is_empty()));
    }
}
