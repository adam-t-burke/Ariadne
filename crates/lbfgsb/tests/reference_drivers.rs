use ariadne_lbfgsb::{
    Backend, Bounds, Control, Convergence, Iteration, Options, SolveAdapter, SolveError, Solver,
    StopReason, Termination,
};
use std::convert::Infallible;

fn reference_last_accepted(csv: &str) -> (usize, usize, f64, f64) {
    let line = csv
        .lines()
        .rev()
        .find(|line| line.starts_with("accepted,"))
        .unwrap();
    let mut fields = line.splitn(6, ',');
    assert_eq!(fields.next(), Some("accepted"));
    (
        fields.next().unwrap().trim().parse().unwrap(),
        fields.next().unwrap().trim().parse().unwrap(),
        fields.next().unwrap().trim().parse().unwrap(),
        fields.next().unwrap().trim().parse().unwrap(),
    )
}

#[derive(Debug)]
struct ReferenceRecord {
    iteration: usize,
    evaluations: usize,
    f: f64,
    pg: f64,
    x: Vec<f64>,
    g: Vec<f64>,
}

fn csv_fields(line: &str) -> Vec<&str> {
    let mut fields = Vec::new();
    let mut quoted = false;
    let mut start = 0;
    for (i, byte) in line.bytes().enumerate() {
        if byte == b'"' {
            quoted = !quoted;
        } else if byte == b',' && !quoted {
            fields.push(line[start..i].trim_matches('"').trim());
            start = i + 1;
        }
    }
    fields.push(line[start..].trim_matches('"').trim());
    fields
}

fn parse_vector(field: &str) -> Vec<f64> {
    if field.is_empty() {
        Vec::new()
    } else {
        field
            .split('|')
            .map(|value| value.trim().parse().unwrap())
            .collect()
    }
}

fn reference_records(csv: &str, kind: &str) -> Vec<ReferenceRecord> {
    csv.lines()
        .skip(1)
        .filter(|line| line.starts_with(kind))
        .map(|line| {
            let fields = csv_fields(line);
            ReferenceRecord {
                iteration: fields[1].parse().unwrap(),
                evaluations: fields[2].parse().unwrap(),
                f: fields[3].parse().unwrap(),
                pg: fields[4].parse().unwrap(),
                x: parse_vector(fields[5]),
                g: parse_vector(fields[6]),
            }
        })
        .collect()
}

fn assert_scaled(actual: f64, expected: f64, tolerance: f64, label: &str) {
    let scale = expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{label}: actual={actual:.17e} expected={expected:.17e}"
    );
}

fn assert_driver_record(iteration: Iteration<'_>, expected: &ReferenceRecord) {
    assert_eq!(iteration.stats.iterations, expected.iteration);
    assert_eq!(iteration.stats.evaluations, expected.evaluations);
    assert_scaled(iteration.value, expected.f, 5e-12, "f");
    assert_scaled(
        iteration.projected_gradient_norm,
        expected.pg,
        5e-11,
        "projected gradient",
    );
    assert_eq!(iteration.x.len(), expected.x.len());
    for (i, (&actual, &reference)) in iteration.x.iter().zip(&expected.x).enumerate() {
        assert_scaled(actual, reference, 5e-11, &format!("x[{i}]"));
    }
    let mut actual_g = vec![0.0; iteration.x.len()];
    driver1_value_gradient(iteration.x, &mut actual_g).unwrap();
    assert_eq!(actual_g.len(), expected.g.len());
    for (i, (&actual, &reference)) in actual_g.iter().zip(&expected.g).enumerate() {
        assert_scaled(actual, reference, 2e-10, &format!("g[{i}]"));
    }
}

fn driver1_value_gradient(x: &[f64], g: &mut [f64]) -> Result<f64, Infallible> {
    let n = x.len();
    let mut f = 0.25 * (x[0] - 1.0).powi(2);
    for i in 1..n {
        let t = x[i] - x[i - 1].powi(2);
        f += t * t;
    }
    f *= 4.0;

    let mut t1 = x[1] - x[0].powi(2);
    g[0] = 2.0 * (x[0] - 1.0) - 16.0 * x[0] * t1;
    for i in 1..n - 1 {
        let t2 = t1;
        t1 = x[i + 1] - x[i].powi(2);
        g[i] = 8.0 * t2 - 16.0 * x[i] * t1;
    }
    g[n - 1] = 8.0 * t1;
    Ok(f)
}

fn audit_value_gradient(x: &[f64], g: &mut [f64]) -> Result<f64, Infallible> {
    const TARGET: [f64; 8] = [1.0, -2.0, 4.0, -1.0, 9.0, 2.0, -1.0, 1.0];
    const SCALE: [f64; 8] = [1.0, 50.0, 0.1, 20.0, 1.0, 0.5, 0.5, 0.5];
    let mut value = 0.0;
    for i in 0..x.len() {
        let residual = x[i] - TARGET[i];
        value += 0.5 * SCALE[i] * residual * residual;
        g[i] = SCALE[i] * residual;
    }
    Ok(value)
}

#[test]
fn driver1_matches_reference_iteration_count() {
    let n = 25;
    let mut x = vec![3.0; n];
    let mut lower = vec![-100.0; n];
    let upper = vec![100.0; n];
    for i in (0..n).step_by(2) {
        lower[i] = 1.0;
    }

    let expected = reference_records(include_str!("reference/fixtures/driver1.csv"), "accepted,");
    let mut record = 0;
    let options = Options::new()
        .with_history_size(5)
        .unwrap()
        .with_relative_function_tolerance(1.0e7 * f64::EPSILON)
        .unwrap()
        .with_projected_gradient_tolerance(1.0e-5)
        .unwrap();
    let report = Solver::new(options)
        .minimize_with_callback(
            &mut x,
            Bounds::new(&lower, &upper, n).unwrap(),
            driver1_value_gradient,
            |iteration| {
                assert_driver_record(iteration, &expected[record]);
                record += 1;
                Control::Continue
            },
        )
        .unwrap();

    let reference = reference_last_accepted(include_str!("reference/fixtures/driver1.csv"));
    assert_eq!(report.stats.iterations, reference.0);
    assert_eq!(report.stats.evaluations, reference.1);
    assert!(report.termination.converged());
    assert!((report.value - reference.2).abs() < 5.0e-9);
    assert_eq!(record, expected.len());
}

#[cfg(feature = "faer-backend")]
#[test]
fn forced_faer_matches_official_driver1_trajectory() {
    let n = 25;
    let mut x = vec![3.0; n];
    let mut lower = vec![-100.0; n];
    let upper = vec![100.0; n];
    for i in (0..n).step_by(2) {
        lower[i] = 1.0;
    }
    let expected = reference_records(include_str!("reference/fixtures/driver1.csv"), "accepted,");
    let options = Options::new()
        .with_history_size(5)
        .unwrap()
        .with_relative_function_tolerance(1.0e7 * f64::EPSILON)
        .unwrap()
        .with_projected_gradient_tolerance(1.0e-5)
        .unwrap()
        .with_backend(Backend::Faer)
        .unwrap();
    let mut record = 0;
    let report = Solver::new(options)
        .minimize_with_callback(
            &mut x,
            Bounds::new(&lower, &upper, n).unwrap(),
            driver1_value_gradient,
            |iteration| {
                assert_driver_record(iteration, &expected[record]);
                record += 1;
                Control::Continue
            },
        )
        .unwrap();

    assert!(report.termination.converged());
    assert_eq!(record, expected.len());
}

#[test]
fn driver2_matches_reference_custom_stop() {
    let n = 25;
    let mut x = vec![3.0; n];
    let mut lower = vec![-100.0; n];
    let upper = vec![100.0; n];
    for i in (0..n).step_by(2) {
        lower[i] = 1.0;
    }
    let expected = reference_records(include_str!("reference/fixtures/driver2.csv"), "accepted,");
    let mut record = 0;
    let options = Options::new()
        .with_history_size(5)
        .unwrap()
        .with_relative_function_tolerance(0.0)
        .unwrap()
        .with_projected_gradient_tolerance(0.0)
        .unwrap();
    let report = Solver::new(options)
        .minimize_with_callback(
            &mut x,
            Bounds::new(&lower, &upper, n).unwrap(),
            driver1_value_gradient,
            |iteration| {
                assert_driver_record(iteration, &expected[record]);
                record += 1;
                if iteration.stats.evaluations >= 99
                    || iteration.projected_gradient_norm <= 1.0e-10 * (1.0 + iteration.value.abs())
                {
                    Control::Stop
                } else {
                    Control::Continue
                }
            },
        )
        .unwrap();
    let reference = reference_last_accepted(include_str!("reference/fixtures/driver2.csv"));
    assert_eq!(report.stats.iterations, reference.0);
    assert_eq!(report.stats.evaluations, reference.1);
    assert_eq!(report.termination, Termination::Stopped(StopReason::User));
    assert!((report.value - reference.2).abs() < 1e-20);
    assert_eq!(record, expected.len());
}

#[test]
fn deterministic_driver3_matches_reference_summary() {
    let n = 1000;
    let mut x = vec![3.0; n];
    let mut lower = vec![-100.0; n];
    let upper = vec![100.0; n];
    for i in (0..n).step_by(2) {
        lower[i] = 1.0;
    }
    let fixture = include_str!("reference/fixtures/driver3-large-n.csv");
    let expected = reference_records(fixture, "accepted,");
    let terminal = reference_records(fixture, "terminal,").pop().unwrap();
    let mut record = 0;
    let options = Options::new()
        .with_history_size(10)
        .unwrap()
        .with_relative_function_tolerance(0.0)
        .unwrap()
        .with_projected_gradient_tolerance(0.0)
        .unwrap()
        .with_backend(Backend::Deterministic)
        .unwrap();
    let mut solver = Solver::new(options);
    let report = solver
        .minimize_with_callback(
            &mut x,
            Bounds::new(&lower, &upper, n).unwrap(),
            driver1_value_gradient,
            |iteration| {
                let expected_iteration = &expected[record];
                assert_eq!(iteration.stats.iterations, expected_iteration.iteration);
                assert_eq!(iteration.stats.evaluations, expected_iteration.evaluations);
                assert_scaled(iteration.value, expected_iteration.f, 5e-11, "driver3 f");
                assert_scaled(
                    iteration.projected_gradient_norm,
                    expected_iteration.pg,
                    5e-10,
                    "driver3 pg",
                );
                assert!(expected_iteration.x.is_empty());
                assert!(expected_iteration.g.is_empty());
                record += 1;
                if iteration.stats.evaluations >= 900
                    || iteration.projected_gradient_norm <= 1.0e-10 * (1.0 + iteration.value.abs())
                {
                    Control::Stop
                } else {
                    Control::Continue
                }
            },
        )
        .unwrap();
    let reference = reference_last_accepted(include_str!("reference/fixtures/driver3-large-n.csv"));
    assert_eq!(
        (report.stats.iterations, report.stats.evaluations),
        (reference.0, reference.1)
    );
    assert!((report.value - reference.2).abs() < 1e-23);
    assert!((report.projected_gradient_norm - reference.3).abs() < 1e-12,);
    assert_eq!(record, expected.len());
    for (i, (&actual, &expected)) in x.iter().zip(&terminal.x).enumerate() {
        assert_scaled(actual, expected, 5e-10, &format!("driver3 x[{i}]"));
    }
    for (i, (&actual, &expected)) in solver
        .workspace()
        .gradient()
        .iter()
        .zip(&terminal.g)
        .enumerate()
    {
        assert_scaled(actual, expected, 2e-9, &format!("driver3 g[{i}]"));
    }
}

#[test]
fn projects_infeasible_start_before_first_evaluation() {
    let mut x = vec![-3.0, 9.0];
    let mut first = None;
    let options = Options::new()
        .with_history_size(5)
        .unwrap()
        .with_projected_gradient_tolerance(1e-10)
        .unwrap();
    let report = Solver::new(options)
        .minimize(
            &mut x,
            Bounds::new(&[0.0, -1.0], &[2.0, 4.0], 2).unwrap(),
            |x, gradient| {
                first.get_or_insert_with(|| x.to_vec());
                gradient[0] = 2.0 * (x[0] - 1.0);
                gradient[1] = 2.0 * (x[1] - 2.0);
                Ok::<_, Infallible>((x[0] - 1.0).powi(2) + (x[1] - 2.0).powi(2))
            },
        )
        .unwrap();

    assert_eq!(first.unwrap(), vec![0.0, 4.0]);
    assert!(report.termination.converged());
}

#[test]
fn mixed_bound_fixture_matches_reference_summary() {
    let mut x = [-3.0, 9.0, 7.0, -8.0];
    let lower = [0.0, -1.0, f64::NEG_INFINITY, 0.0];
    let upper = [2.0, 4.0, f64::INFINITY, f64::INFINITY];
    let bounds = Bounds::new(&lower, &upper, x.len()).unwrap();
    let options = Options::new()
        .with_history_size(5)
        .unwrap()
        .with_projected_gradient_tolerance(1e-12)
        .unwrap();
    let mut solver = Solver::new(options);
    let mut feasible_callbacks = true;
    let report = solver
        .minimize(&mut x, bounds, |x, g| {
            feasible_callbacks &= x
                .iter()
                .zip(lower)
                .zip(upper)
                .all(|((&x, l), u)| x >= l && x <= u);
            let target = [1.0, 2.0, -2.0, 3.0];
            let mut f = 0.0;
            for i in 0..4 {
                let d = x[i] - target[i];
                f += d * d;
                g[i] = 2.0 * d;
            }
            Ok::<_, Infallible>(f)
        })
        .unwrap();
    assert!(feasible_callbacks);
    let reference =
        reference_last_accepted(include_str!("reference/fixtures/edge-mixed-bounds.csv"));
    assert_eq!(
        (report.stats.iterations, report.stats.evaluations),
        (reference.0, reference.1)
    );
    assert!((report.value - reference.2).abs() < 1e-26);
    assert!(x
        .iter()
        .zip([1.0, 2.0, -2.0, 3.0])
        .all(|(&x, target)| (x - target).abs() < 1e-12));
}

#[test]
fn mixed_bound_ties_and_rollover_match_official_trajectory() {
    let mut x = [-3.0, 8.0, 4.0, -6.0, 9.0, 0.0, 0.0, 0.0];
    let lower = [
        0.0,
        0.0,
        f64::NEG_INFINITY,
        f64::NEG_INFINITY,
        1.5,
        -1.0,
        -0.5,
        f64::NEG_INFINITY,
    ];
    let upper = [
        2.0,
        f64::INFINITY,
        2.0,
        f64::INFINITY,
        1.5,
        1.0,
        f64::INFINITY,
        0.5,
    ];
    let expected = reference_records(
        include_str!("reference/fixtures/audit-mixed-rollover.csv"),
        "accepted,",
    );
    let mut first_evaluation = None;
    let mut record = 0;
    let options = Options::new()
        .with_history_size(2)
        .unwrap()
        .with_projected_gradient_tolerance(1.0e-10)
        .unwrap()
        .with_backend(Backend::Deterministic)
        .unwrap();
    let report = Solver::new(options)
        .minimize_with_callback(
            &mut x,
            Bounds::new(&lower, &upper, 8).unwrap(),
            |x, gradient| {
                first_evaluation.get_or_insert_with(|| x.to_vec());
                audit_value_gradient(x, gradient)
            },
            |iteration| {
                let reference = &expected[record];
                assert_eq!(iteration.stats.iterations, reference.iteration);
                assert_eq!(iteration.stats.evaluations, reference.evaluations);
                assert_scaled(iteration.value, reference.f, 5e-12, "audit f");
                assert_scaled(
                    iteration.projected_gradient_norm,
                    reference.pg,
                    5e-11,
                    "audit projected gradient",
                );
                for (i, (&actual, &expected)) in iteration.x.iter().zip(&reference.x).enumerate() {
                    assert_scaled(actual, expected, 5e-11, &format!("audit x[{i}]"));
                }
                let mut gradient = [0.0; 8];
                audit_value_gradient(iteration.x, &mut gradient).unwrap();
                for (i, (&actual, &expected)) in gradient.iter().zip(&reference.g).enumerate() {
                    assert_scaled(actual, expected, 2e-10, &format!("audit g[{i}]"));
                }
                record += 1;
                Control::Continue
            },
        )
        .unwrap();

    assert_eq!(
        first_evaluation.unwrap(),
        [0.0, 8.0, 2.0, -6.0, 1.5, 0.0, 0.0, 0.0]
    );
    assert_eq!(record, expected.len());
    assert_eq!((report.stats.iterations, report.stats.evaluations), (8, 11));
    assert!(
        report.stats.accepted_updates > options.history_size(),
        "the two-pair ring must roll over"
    );
    assert_eq!(x[4], 1.5, "the fixed variable must remain projected");
    assert!(report.termination.converged());
}

#[test]
fn supports_mixed_and_unbounded_variables() {
    let mut x = vec![5.0, -4.0, 3.0];
    let lower = [f64::NEG_INFINITY, 0.0, f64::NEG_INFINITY];
    let upper = [f64::INFINITY, 2.0, 4.0];
    let options = Options::new()
        .with_history_size(5)
        .unwrap()
        .with_projected_gradient_tolerance(1e-10)
        .unwrap();
    let report = Solver::new(options)
        .minimize(
            &mut x,
            Bounds::new(&lower, &upper, 3).unwrap(),
            |x, gradient| {
                let targets = [1.0, 1.0, -2.0];
                for ((gradient, xi), target) in gradient.iter_mut().zip(x).zip(targets) {
                    *gradient = 2.0 * (xi - target);
                }
                let f = x
                    .iter()
                    .zip(targets)
                    .map(|(xi, ti)| (xi - ti).powi(2))
                    .sum();
                Ok::<_, Infallible>(f)
            },
        )
        .unwrap();

    assert!(matches!(report.termination, Termination::Converged(_)));
    assert!((x[0] - 1.0).abs() < 1e-7);
    assert!((x[1] - 1.0).abs() < 1e-7);
    assert!((x[2] + 2.0).abs() < 1e-7);
}

#[test]
fn relative_function_tolerance_is_reported() {
    let mut x = vec![10.0, -5.0];
    let lower = [f64::NEG_INFINITY; 2];
    let upper = [f64::INFINITY; 2];
    let bounds = Bounds::new(&lower, &upper, x.len()).unwrap();
    let options = Options::new()
        .with_history_size(5)
        .unwrap()
        .with_projected_gradient_tolerance(0.0)
        .unwrap()
        .with_relative_function_tolerance(1.0e15 * f64::EPSILON)
        .unwrap();
    let report = Solver::new(options)
        .minimize(&mut x, bounds, |x, gradient| {
            gradient[0] = 4.0 * (x[0] - 1.0).powi(3);
            gradient[1] = 4.0 * (x[1] + 2.0).powi(3);
            Ok::<_, Infallible>((x[0] - 1.0).powi(4) + (x[1] + 2.0).powi(4) + 3.0)
        })
        .unwrap();

    assert_eq!(
        report.termination,
        Termination::Converged(Convergence::RelativeFunction)
    );
    assert!(report.value >= 3.0);
}

#[test]
fn fixed_variables_and_history_rollover_are_supported() {
    let options = Options::new()
        .with_history_size(2)
        .unwrap()
        .with_relative_function_tolerance(0.0)
        .unwrap()
        .with_projected_gradient_tolerance(1e-10)
        .unwrap();
    let mut solver = Solver::new(options);
    let mut x = [5.0, 2.0, -4.0, 8.0];
    let bounds = Bounds::new(&[-10.0, 2.0, -10.0, -10.0], &[10.0, 2.0, 10.0, 10.0], 4).unwrap();
    let report = solver
        .minimize(&mut x, bounds, |x, g| {
            let targets = [1.0, 99.0, -2.0, 3.0];
            let mut f = 0.0;
            for i in 0..x.len() {
                let d = x[i] - targets[i];
                f += (i + 1) as f64 * d * d;
                g[i] = 2.0 * (i + 1) as f64 * d;
            }
            Ok::<_, Infallible>(f)
        })
        .unwrap();
    assert!(report.termination.converged());
    assert_eq!(x[1], 2.0);
    assert!((x[0] - 1.0).abs() < 1e-7);
    assert!((x[2] + 2.0).abs() < 1e-7);
    assert!((x[3] - 3.0).abs() < 1e-7);
}

#[test]
fn line_search_limit_rolls_back_last_accepted_iterate() {
    let options = Options::new()
        .with_max_line_search_evaluations(1)
        .unwrap()
        .with_projected_gradient_tolerance(0.0)
        .unwrap();
    let mut solver = Solver::new(options);
    let mut x = [0.5];
    let bounds = Bounds::new(&[0.0], &[1.0], 1).unwrap();
    let report = solver
        .minimize(&mut x, bounds, |x, g| {
            g[0] = -1.0;
            Ok::<_, Infallible>(if x[0] == 0.5 { 0.0 } else { 1.0 })
        })
        .unwrap();
    assert_eq!(
        report.termination,
        Termination::Stopped(StopReason::MaximumLineSearchEvaluations)
    );
    assert_eq!(report.stats.evaluations, 2);
    assert_eq!(x, [0.5]);
    assert_eq!(solver.workspace().gradient(), [-1.0]);
}

#[test]
fn iteration_limit_above_one_thousand_is_exact() {
    let options = Options::new()
        .with_max_iterations(1001)
        .unwrap()
        .with_max_line_search_evaluations(20)
        .unwrap()
        .with_relative_function_tolerance(0.0)
        .unwrap()
        .with_projected_gradient_tolerance(0.0)
        .unwrap();
    let mut solver = Solver::new(options);
    let mut x = [0.0];
    let bounds = Bounds::new(&[f64::NEG_INFINITY], &[f64::INFINITY], 1).unwrap();
    let report = solver
        .minimize_with_callback(
            &mut x,
            bounds,
            |x, g| {
                g[0] = -1.0;
                Ok::<_, Infallible>(-x[0])
            },
            |_| Control::Continue,
        )
        .unwrap();
    assert_eq!(
        report.termination,
        Termination::Stopped(StopReason::MaximumIterations)
    );
    assert_eq!(report.stats.iterations, 1001);
}

#[test]
fn nonfinite_projected_initial_value_is_rejected() {
    let mut x = [f64::INFINITY];
    let bounds = Bounds::new(&[0.0], &[1.0], 1).unwrap();
    let error = Solver::new(Options::new())
        .minimize(&mut x, bounds, |_x, _g| Ok::<_, Infallible>(0.0))
        .unwrap_err();
    assert!(matches!(
        error,
        SolveError::InvalidInitialValue { index: 0 }
    ));
}

#[test]
fn trial_evaluation_error_rolls_back_iterate_and_gradient() {
    let options = Options::new()
        .with_projected_gradient_tolerance(0.0)
        .unwrap();
    let mut solver = Solver::new(options);
    let mut x = [0.5];
    let bounds = Bounds::new(&[0.0], &[1.0], 1).unwrap();
    solver
        .minimize(&mut x, bounds, |x, gradient| {
            gradient[0] = 2.0 * (x[0] - 0.25);
            Ok::<_, Infallible>((x[0] - 0.25).powi(2))
        })
        .unwrap();
    let committed = solver.workspace().gradient().to_vec();
    x = [0.5];
    let mut evaluations = 0;
    let error = solver
        .minimize(&mut x, bounds, |_x, gradient| {
            evaluations += 1;
            gradient[0] = -1.0;
            if evaluations == 1 {
                Ok(0.0)
            } else {
                Err("trial failed")
            }
        })
        .unwrap_err();

    assert!(matches!(error, SolveError::Objective("trial failed")));
    assert_eq!(x, [0.5]);
    assert_eq!(solver.workspace().gradient(), committed);
}

struct ObserverAdapter {
    fail: bool,
}

impl SolveAdapter for ObserverAdapter {
    type Error = &'static str;

    fn value_and_gradient(&mut self, x: &[f64], gradient: &mut [f64]) -> Result<f64, Self::Error> {
        gradient[0] = 2.0 * (x[0] - 1.0);
        Ok((x[0] - 1.0).powi(2))
    }

    fn accepted_iteration(&mut self, _iteration: Iteration<'_>) -> Result<Control, Self::Error> {
        if self.fail {
            Err("observer failed")
        } else {
            Ok(Control::Continue)
        }
    }
}

#[test]
fn observer_error_preserves_previous_successful_gradient() {
    let mut solver = Solver::new(Options::new());
    let bounds = Bounds::new(&[-10.0], &[10.0], 1).unwrap();
    let mut x = [4.0];
    solver
        .minimize_with_adapter(&mut x, bounds, &mut ObserverAdapter { fail: false })
        .unwrap();
    let committed = solver.workspace().gradient().to_vec();

    x = [4.0];
    let error = solver
        .minimize_with_adapter(&mut x, bounds, &mut ObserverAdapter { fail: true })
        .unwrap_err();

    assert!(matches!(error, SolveError::Callback("observer failed")));
    assert_eq!(solver.workspace().gradient(), committed);
}

#[test]
fn zero_projected_gradient_tolerance_is_disabled() {
    let options = Options::new()
        .with_relative_function_tolerance(0.0)
        .unwrap()
        .with_projected_gradient_tolerance(0.0)
        .unwrap();
    let mut x = [0.0];
    let report = Solver::new(options)
        .minimize(
            &mut x,
            Bounds::new(&[-1.0], &[1.0], 1).unwrap(),
            |_x, gradient| {
                gradient[0] = 0.0;
                Ok::<_, Infallible>(0.0)
            },
        )
        .unwrap();

    assert_eq!(
        report.termination,
        Termination::Failed(ariadne_lbfgsb::Failure::NoDescentDirection)
    );
}
