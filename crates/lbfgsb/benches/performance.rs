use ariadne_lbfgsb::{Backend, Bounds, Options, Solver, Stats, Termination};
use std::convert::Infallible;
use std::hint::black_box;
use std::time::{Duration, Instant};

const DIMENSIONS: [usize; 6] = [25, 180, 512, 1000, 1984, 4096];
const HISTORIES: [usize; 3] = [5, 10, 20];

#[derive(Clone, Copy, Debug)]
enum Objective {
    Diagonal,
    Rotated,
    Rosenbrock,
    ActiveCorner,
    MixedBounds,
    TiedBreakpoints,
}

impl Objective {
    fn name(self) -> &'static str {
        match self {
            Self::Diagonal => "diagonal",
            Self::Rotated => "rotated",
            Self::Rosenbrock => "rosenbrock",
            Self::ActiveCorner => "active-corner",
            Self::MixedBounds => "mixed-bounds",
            Self::TiedBreakpoints => "tied-breakpoints",
        }
    }

    fn from_name(name: &str) -> Option<Self> {
        match name {
            "diagonal" => Some(Self::Diagonal),
            "rotated" => Some(Self::Rotated),
            "rosenbrock" => Some(Self::Rosenbrock),
            "active-corner" => Some(Self::ActiveCorner),
            "mixed-bounds" => Some(Self::MixedBounds),
            "tied-breakpoints" => Some(Self::TiedBreakpoints),
            _ => None,
        }
    }

    fn evaluate(self, x: &[f64], g: &mut [f64]) -> f64 {
        match self {
            Self::Diagonal => {
                let mut f = 0.0;
                for i in 0..x.len() {
                    let d = x[i] - ((i % 11) as f64 - 5.0) * 0.05;
                    let scale = 1.0 + (i % 97) as f64;
                    f += 0.5 * scale * d * d;
                    g[i] = scale * d;
                }
                f
            }
            Self::Rotated => {
                let mut f = 0.0;
                g.fill(0.0);
                for i in 0..x.len() {
                    let j = (i + 1) % x.len();
                    let d = x[i] - 0.25;
                    let coupled = x[i] + 0.35 * x[j] - 0.1;
                    f += 0.5 * d * d + 2.0 * coupled * coupled;
                    g[i] += d + 4.0 * coupled;
                    g[j] += 1.4 * coupled;
                }
                f
            }
            Self::Rosenbrock => {
                let mut f = 0.25 * (x[0] - 1.0).powi(2);
                for i in 1..x.len() {
                    let r = x[i] - x[i - 1].powi(2);
                    f += r * r;
                }
                f *= 4.0;
                let mut r = x[1] - x[0].powi(2);
                g[0] = 2.0 * (x[0] - 1.0) - 16.0 * x[0] * r;
                for i in 1..x.len() - 1 {
                    let previous = r;
                    r = x[i + 1] - x[i].powi(2);
                    g[i] = 8.0 * previous - 16.0 * x[i] * r;
                }
                g[x.len() - 1] = 8.0 * r;
                f
            }
            Self::ActiveCorner => {
                let mut f = 0.0;
                for i in 0..x.len() {
                    let target = if i % 2 == 0 { -3.0 } else { 4.0 };
                    let d = x[i] - target;
                    f += 0.5 * d * d;
                    g[i] = d;
                }
                f
            }
            Self::MixedBounds => {
                let mut f = 0.0;
                for i in 0..x.len() {
                    let target = (i % 7) as f64 * 0.2 - 0.6;
                    let d = x[i] - target;
                    let scale = 1.0 + (i % 13) as f64;
                    f += 0.5 * scale * d * d;
                    g[i] = scale * d;
                }
                f
            }
            Self::TiedBreakpoints => {
                let mut f = 0.0;
                for i in 0..x.len() {
                    let target = if i % 2 == 0 { -2.0 } else { 2.0 };
                    let d = x[i] - target;
                    f += 0.5 * d * d;
                    g[i] = d;
                }
                f
            }
        }
    }
}

struct Problem {
    objective: Objective,
    initial: Vec<f64>,
    lower: Vec<f64>,
    upper: Vec<f64>,
}

fn problem(objective: Objective, n: usize) -> Problem {
    let mut lower = vec![-1.0; n];
    let mut upper = vec![1.0; n];
    let initial = match objective {
        Objective::Rosenbrock => vec![3.0; n],
        Objective::TiedBreakpoints => vec![0.0; n],
        _ => vec![0.75; n],
    };
    match objective {
        Objective::Diagonal | Objective::Rotated => {
            lower.fill(f64::NEG_INFINITY);
            upper.fill(f64::INFINITY);
        }
        Objective::Rosenbrock => {
            lower.fill(-100.0);
            upper.fill(100.0);
            for i in (0..n).step_by(2) {
                lower[i] = 1.0;
            }
        }
        Objective::MixedBounds => {
            for i in 0..n {
                match i % 4 {
                    0 => {
                        lower[i] = f64::NEG_INFINITY;
                        upper[i] = f64::INFINITY;
                    }
                    1 => upper[i] = f64::INFINITY,
                    2 => lower[i] = f64::NEG_INFINITY,
                    _ => {}
                }
            }
        }
        Objective::ActiveCorner | Objective::TiedBreakpoints => {}
    }
    Problem {
        objective,
        initial,
        lower,
        upper,
    }
}

fn options(m: usize, backend: Backend) -> Options {
    Options::new()
        .with_history_size(m)
        .unwrap()
        .with_projected_gradient_tolerance(1.0e-5)
        .unwrap()
        .with_relative_function_tolerance(1.0e7 * f64::EPSILON)
        .unwrap()
        .with_backend(backend)
        .unwrap()
}

#[derive(Clone, Copy, Debug)]
struct Outcome {
    stats: Stats,
    termination: Termination,
    value: f64,
    projected_gradient: f64,
    checksum: f64,
}

fn high_level_once(solver: &mut Solver, problem: &Problem, x: &mut [f64]) -> Outcome {
    x.copy_from_slice(&problem.initial);
    let bounds = Bounds::new(&problem.lower, &problem.upper, x.len()).unwrap();
    let report = solver
        .minimize(x, bounds, |x, g| {
            Ok::<_, Infallible>(problem.objective.evaluate(x, g))
        })
        .unwrap();
    let checksum = x
        .iter()
        .enumerate()
        .map(|(index, value)| (index + 1) as f64 * value)
        .sum();
    black_box(Outcome {
        stats: report.stats,
        termination: report.termination,
        value: report.value,
        projected_gradient: report.projected_gradient_norm,
        checksum,
    })
}

fn assert_equivalent(actual: Outcome, expected: Outcome) {
    assert_eq!(actual.termination, expected.termination);
    for (label, actual, expected) in [
        ("value", actual.value, expected.value),
        (
            "projected gradient",
            actual.projected_gradient,
            expected.projected_gradient,
        ),
        ("x checksum", actual.checksum, expected.checksum),
    ] {
        let tolerance = 1.0e-10 * expected.abs().max(1.0);
        assert!(
            (actual - expected).abs() <= tolerance,
            "{label} differs: actual={actual:.17e}, expected={expected:.17e}"
        );
    }
}

fn measure(objective: Objective, n: usize, m: usize, minimum: Duration, backend: Backend) {
    let problem = problem(objective, n);
    let mut x = problem.initial.clone();
    let mut solver = Solver::new(options(m, backend));
    let mut deterministic = Solver::new(options(m, Backend::Deterministic));
    let expected = high_level_once(&mut deterministic, &problem, &mut x);
    let actual = high_level_once(&mut solver, &problem, &mut x);
    assert_equivalent(actual, expected);
    let resolved = solver.resolved_backend();

    let start = Instant::now();
    let mut runs = 0u64;
    let mut last = actual;
    while runs < 3 || start.elapsed() < minimum {
        last = high_level_once(&mut solver, &problem, &mut x);
        runs += 1;
    }
    let elapsed = start.elapsed();
    println!(
        "{:?},{},{},{},{},{:.3},{},{},{:?},{:.17e},{:.17e},{:.17e}",
        resolved,
        objective.name(),
        n,
        m,
        runs,
        elapsed.as_secs_f64() * 1.0e6 / runs as f64,
        last.stats.iterations,
        last.stats.evaluations,
        last.termination,
        last.value,
        last.projected_gradient,
        last.checksum,
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let full = args.iter().any(|arg| arg == "--full");
    let backend = args
        .iter()
        .find_map(|arg| arg.strip_prefix("--backend="))
        .map(|value| match value {
            "auto" => Backend::Auto,
            "deterministic" | "scalar" => Backend::Deterministic,
            "faer" => Backend::Faer,
            other => panic!("unknown backend {other:?}; use auto, deterministic, or faer"),
        })
        .unwrap_or(Backend::Deterministic);
    let objective_filter = args
        .iter()
        .find_map(|arg| arg.strip_prefix("--objective="))
        .map(|value| {
            Objective::from_name(value).unwrap_or_else(|| {
                panic!(
                    "unknown objective {value:?}; use diagonal, rotated, rosenbrock, \
                     active-corner, mixed-bounds, or tied-breakpoints"
                )
            })
        });
    let dimension_filter = args
        .iter()
        .find_map(|arg| arg.strip_prefix("--n="))
        .map(|value| value.parse::<usize>().expect("--n must be a positive integer"));
    let history_filter = args
        .iter()
        .find_map(|arg| arg.strip_prefix("--m="))
        .map(|value| value.parse::<usize>().expect("--m must be a positive integer"));
    let millis = std::env::var("LBFGSB_BENCH_MILLIS")
        .ok()
        .and_then(|value| value.parse().ok())
        .unwrap_or(250);
    let minimum = Duration::from_millis(millis);
    let objectives = [
        Objective::Diagonal,
        Objective::Rotated,
        Objective::Rosenbrock,
        Objective::ActiveCorner,
        Objective::MixedBounds,
        Objective::TiedBreakpoints,
    ];

    println!(
        "resolved_backend,objective,n,m,runs,microseconds,iterations,evaluations,termination,final_f,final_pg,x_checksum"
    );
    if objective_filter.is_some() || dimension_filter.is_some() || history_filter.is_some() {
        let selected_objectives: Vec<_> = objective_filter
            .map(|objective| vec![objective])
            .unwrap_or_else(|| objectives.to_vec());
        let selected_dimensions: Vec<_> = dimension_filter
            .map(|dimension| vec![dimension])
            .unwrap_or_else(|| DIMENSIONS.to_vec());
        let selected_histories: Vec<_> = history_filter
            .map(|history| vec![history])
            .unwrap_or_else(|| HISTORIES.to_vec());
        for objective in selected_objectives {
            for &n in &selected_dimensions {
                for &m in &selected_histories {
                    measure(objective, n, m, minimum, backend);
                }
            }
        }
    } else if full {
        for objective in objectives {
            for n in DIMENSIONS {
                for m in HISTORIES {
                    measure(objective, n, m, minimum, backend);
                }
            }
        }
    } else {
        let dimensions = [25, 180, 512, 1000, 1984, 4096];
        for (index, objective) in objectives.into_iter().enumerate() {
            measure(
                objective,
                dimensions[index],
                HISTORIES[index % 3],
                minimum,
                backend,
            );
        }
    }
}
