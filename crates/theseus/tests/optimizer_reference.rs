//! Basin comparison against the independently generated L-BFGS-B 3.0 fixtures.
//! Objective definitions and provenance: ../../lbfgsb/tests/reference/PROVENANCE.md.

use std::time::{Duration, Instant};

struct Case {
    name: &'static str,
    start: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
    history: usize,
    gradient_tolerance: f64,
    cost_tolerance: f64,
    evaluation_limit: usize,
    fixture: &'static str,
    check_feasibility: bool,
}

impl Case {
    fn evaluate(&self, x: &[f64], g: &mut [f64]) -> f64 {
        if self.check_feasibility {
            for ((&value, &lower), &upper) in x.iter().zip(&self.lower).zip(&self.upper) {
                assert!(
                    value >= lower && value <= upper,
                    "infeasible trial evaluation"
                );
            }
        }

        if self.name == "audit" {
            let target = [1.0, -2.0, 4.0, -1.0, 9.0, 2.0, -1.0, 1.0];
            let scale = [1.0, 50.0, 0.1, 20.0, 1.0, 0.5, 0.5, 0.5];
            let mut value = 0.0;
            for i in 0..x.len() {
                let d = x[i] - target[i];
                value += 0.5 * scale[i] * d * d;
                g[i] = scale[i] * d;
            }
            return value;
        }
        if self.name == "mixed" {
            let target = [1.0, 2.0, -2.0, 3.0];
            let mut value = 0.0;
            for i in 0..x.len() {
                let d = x[i] - target[i];
                value += d * d;
                g[i] = 2.0 * d;
            }
            return value;
        }
        let n = x.len();
        let mut value = 0.25 * (x[0] - 1.0).powi(2);
        for i in 1..n {
            value += (x[i] - x[i - 1].powi(2)).powi(2);
        }
        let mut t1 = x[1] - x[0].powi(2);
        g[0] = 2.0 * (x[0] - 1.0) - 16.0 * x[0] * t1;
        for i in 1..n - 1 {
            let t2 = t1;
            t1 = x[i + 1] - x[i].powi(2);
            g[i] = 8.0 * t2 - 16.0 * x[i] * t1;
        }
        g[n - 1] = 8.0 * t1;
        4.0 * value
    }
}

fn cases() -> Vec<Case> {
    let fixtures = [
        include_str!("../../lbfgsb/tests/reference/fixtures/driver1.csv"),
        include_str!("../../lbfgsb/tests/reference/fixtures/driver2.csv"),
        include_str!("../../lbfgsb/tests/reference/fixtures/driver3-large-n.csv"),
    ];
    let mut cases = Vec::new();
    for (index, name) in ["driver1", "driver2", "driver3"].into_iter().enumerate() {
        let n = if index == 2 { 1000 } else { 25 };
        let mut lower = vec![-100.0; n];
        for i in (0..n).step_by(2) {
            lower[i] = 1.0;
        }
        cases.push(Case {
            name,
            start: vec![3.0; n],
            lower,
            upper: vec![100.0; n],
            history: if index == 2 { 10 } else { 5 },
            gradient_tolerance: if index == 0 { 1e-5 } else { 0.0 },
            cost_tolerance: if index == 0 { 1e7 * f64::EPSILON } else { 0.0 },
            evaluation_limit: if index == 2 { 900 } else { 99 },
            fixture: fixtures[index],
            check_feasibility: false,
        });
    }
    cases.push(Case {
        name: "mixed",
        start: vec![-3.0, 9.0, 7.0, -8.0],
        lower: vec![0.0, -1.0, f64::NEG_INFINITY, 0.0],
        upper: vec![2.0, 4.0, f64::INFINITY, f64::INFINITY],
        history: 5,
        gradient_tolerance: 1e-12,
        cost_tolerance: 0.0,
        evaluation_limit: 1000,
        fixture: include_str!("../../lbfgsb/tests/reference/fixtures/edge-mixed-bounds.csv"),
        check_feasibility: false,
    });
    cases.push(Case {
        name: "audit",
        start: vec![-3.0, 8.0, 4.0, -6.0, 9.0, 0.0, 0.0, 0.0],
        lower: vec![
            0.0,
            0.0,
            f64::NEG_INFINITY,
            f64::NEG_INFINITY,
            1.5,
            -1.0,
            -0.5,
            f64::NEG_INFINITY,
        ],
        upper: vec![
            2.0,
            f64::INFINITY,
            2.0,
            f64::INFINITY,
            1.5,
            1.0,
            f64::INFINITY,
            0.5,
        ],
        history: 2,
        gradient_tolerance: 1e-10,
        cost_tolerance: 1e7 * f64::EPSILON,
        evaluation_limit: 1000,
        fixture: include_str!("../../lbfgsb/tests/reference/fixtures/audit-mixed-rollover.csv"),
        check_feasibility: false,
    });
    cases
}

#[derive(Debug)]
struct Outcome {
    x: Vec<f64>,
    value: f64,
    projected_gradient: f64,
    iterations: usize,
    evaluations: usize,
}

#[path = "support/basin_reference.rs"]
mod backend;
use backend::solve;

#[test]
fn bounded_reference_solutions() {
    for mut case in cases() {
        case.check_feasibility = true;
        let result = solve(&case);
        let row = case
            .fixture
            .lines()
            .rev()
            .find(|line| line.starts_with("accepted,"))
            .unwrap();
        let fields: Vec<_> = row.split(',').collect();
        let expected_value: f64 = fields[3].trim().parse().unwrap();
        let tolerance = match case.name {
            "driver1" => 5e-9,
            "audit" => 5e-12 * expected_value.abs().max(1.0),
            _ => 1e-18,
        };
        assert!(
            (result.value - expected_value).abs() <= tolerance,
            "{}: {result:?}",
            case.name
        );
        let pg_limit = if case.name == "driver1" { 1e-3 } else { 1e-8 };
        assert!(
            result.projected_gradient <= pg_limit,
            "{}: {result:?}",
            case.name
        );
        for ((&x, &lower), &upper) in result.x.iter().zip(&case.lower).zip(&case.upper) {
            assert!(x.is_finite() && x >= lower && x <= upper);
        }
        if case.name == "audit" {
            assert_eq!(result.x[4], 1.5);
            assert!(result.iterations > case.history, "history must roll over");
        }
        if case.name == "mixed" {
            for (&x, target) in result.x.iter().zip([1.0, 2.0, -2.0, 3.0]) {
                assert!((x - target).abs() < 1e-12);
            }
        }
    }
}

#[test]
#[ignore = "manual release benchmark"]
fn bench_reference_comparison() {
    for case in cases() {
        std::hint::black_box(solve(&case));
        let mut times = Vec::new();
        for _ in 0..10 {
            let start = Instant::now();
            let mut runs = 0;
            // Batch short solves so scheduling and clock overhead do not dominate.
            while start.elapsed() < Duration::from_millis(50) {
                std::hint::black_box(solve(&case));
                runs += 1;
            }
            times.push(start.elapsed().as_secs_f64() * 1e6 / runs as f64);
        }
        times.sort_by(f64::total_cmp);
        let result = solve(&case);
        eprintln!("{},median_us={:.6},loss={:.12e},projected_gradient={:.6e},iterations={},evaluations={}",
            case.name, (times[4] + times[5]) / 2.0, result.value, result.projected_gradient, result.iterations, result.evaluations);
    }
}
