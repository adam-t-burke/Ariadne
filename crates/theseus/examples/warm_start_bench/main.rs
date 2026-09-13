//! Reproducible benchmark for the inverse-FDM warm-start pipeline.
//!
//! Each case builds a synthetic net with a known funicular (see `nets.rs`),
//! perturbs it into a target, computes a warm start with every method in
//! `methods.rs`, clips the result to the L-BFGS-B box, and runs box-constrained
//! L-BFGS-B on the SSE target objective. Tables are printed to stdout; every
//! case has one header line starting with `config/box` and one row per method.
//!
//! Build:  `cargo build --release -p theseus --example warm_start_bench`
//! Run:    `cargo run --release -p theseus --example warm_start_bench -- <subcommand>`
//!
//! Subcommands
//! - `nets`               sanity check of every generator (funicular finite,
//!                        sign counts, depth/extent, max support reaction);
//!                        `BENCH_VERBOSE=1` adds coordinate ranges per net.
//! - `suite [max_iters]`  all nets × {jit2%d, bump10%d} × {loose, snug} × all
//!                        methods. `BENCH_NETS=a,b` and `BENCH_METHODS=x,y`
//!                        restrict the run.
//! - `alt`                `suite` restricted to seven methods and four nets.
//! - `scale`              corner-anchored quads with ≈1k, 4k, 16k, 65k edges:
//!                        warm-start time per method, forward-solve and
//!                        L-BFGS-B evaluation time.
//! - `dense`              sparse vs dense Gram normal equations, ne ≈ 250..4000,
//!                        with the observed time exponent of the dense path.
//! - `reactions`          tied arch and cable truss with and without the
//!                        zero-horizontal-reaction rows, plus the inconsistent
//!                        zero-vertical-reaction request.

mod methods;
mod nets;
mod report;

use methods::{
    clip, dense_gram_solve, diag_note, equilibrium_e, human_bytes, lbfgsb, library_options,
    run_library, run_method, time_eval, time_forward, Method, DENSE_EDGE_CAP, GRAM_RELATIVE_SHIFT,
};
use ndarray::Array2;
use nets::{
    all_nets, cable_truss, net_reaction, quad, safe_err, suite_nets, support_reactions, tied_arch,
    tied_arch_straight, try_forward, Net,
};
use report::{fmt_e, print_case, Row};
use std::time::Instant;
use theseus::types::SolverOptions;

const DEFAULT_MAX_ITERS: usize = 1000;
/// L-BFGS-B is skipped in `scale` above this many edges.
const SCALE_LBFGSB_EDGE_CAP: usize = 16_384;
const SCALE_LBFGSB_ITERS: usize = 200;
const SCALE_TIME_BUDGET_S: f64 = 600.0;

fn env_list(var: &str) -> Option<Vec<String>> {
    std::env::var(var)
        .ok()
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.split(',').map(|v| v.trim().to_string()).collect())
}

/// Methods to run: `BENCH_METHODS` if set, else `default`.
fn selected_methods(default: &[Method]) -> Vec<Method> {
    match env_list("BENCH_METHODS") {
        Some(names) => names
            .iter()
            .filter_map(|n| {
                let m = Method::from_name(n);
                if m.is_none() {
                    eprintln!("unknown method '{n}' in BENCH_METHODS");
                }
                m
            })
            .collect(),
        None => default.to_vec(),
    }
}

/// Nets to run: `BENCH_NETS` if set, else `default` (or all when `None`).
fn selected_nets(nets: Vec<Net>, default: Option<&[&str]>) -> Vec<Net> {
    let filter: Option<Vec<String>> = env_list("BENCH_NETS")
        .or_else(|| default.map(|d| d.iter().map(|s| s.to_string()).collect()));
    match filter {
        Some(names) => nets
            .into_iter()
            .filter(|n| names.iter().any(|f| *f == n.name))
            .collect(),
        None => nets,
    }
}

fn targets(net: &Net) -> Vec<(&'static str, Array2<f64>)> {
    vec![
        ("jit2%d", net.target(0.02, 0.0)),
        ("bump10%d", net.target(0.0, 0.10)),
    ]
}

fn boxes(net: &Net) -> Vec<(&'static str, Vec<f64>, Vec<f64>)> {
    let (l1, h1) = net.box_from_true(100.0, 100.0);
    let (l2, h2) = net.box_from_true(1.0, 1.0);
    vec![("loose", l1, h1), ("snug", l2, h2)]
}

/// Warm start → clip → score → L-BFGS-B for every method; one row each.
fn run_case(
    net: &Net,
    target: &Array2<f64>,
    lo: &[f64],
    hi: &[f64],
    methods: &[Method],
    max_iters: usize,
) -> Vec<Row> {
    let fixed = net.fixed_positions();
    let problem = net.problem(&fixed, Vec::new(), lo, hi, SolverOptions::default());
    methods
        .iter()
        .map(|&m| {
            let label = m.name();
            let w = match run_method(m, net, &problem, target, lo, hi) {
                Ok(w) => w,
                Err(e) => return Row::failed(label, 0.0, e),
            };
            let err_raw = safe_err(&problem, target, &w.q);
            let (q, n_clipped, _) = clip(&w.q, lo, hi);
            let err_clip = safe_err(&problem, target, &q);
            let mut row = Row {
                label: label.to_string(),
                warm_ms: w.ms,
                err_raw,
                err_clip,
                n_clipped,
                run: None,
                note: w.note,
            };
            if !err_clip.is_finite() {
                row.note = join_note(&row.note, "clipped q gives singular forward solve");
                return row;
            }
            match lbfgsb(net, &fixed, target, &q, lo, hi, max_iters) {
                Ok(run) => row.run = Some(run),
                Err(e) => row.note = join_note(&row.note, &e),
            }
            row
        })
        .collect()
}

fn join_note(a: &str, b: &str) -> String {
    if a.is_empty() {
        b.to_string()
    } else {
        format!("{a} | {b}")
    }
}

// ───────────────────────── nets ─────────────────────────

/// Coordinate ranges of the funicular, split into loaded and unloaded free
/// nodes (e.g. ridge vs hoop nodes of the cable dome), plus the anchor box.
fn describe_shape(net: &Net, x: &Array2<f64>) {
    let range = |rows: &[usize], d: usize| -> (f64, f64) {
        rows.iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &i| {
                (lo.min(x[[i, d]]), hi.max(x[[i, d]]))
            })
    };
    let loaded: Vec<usize> = (0..net.free.len())
        .filter(|&i| net.loads[i] != 0.0)
        .collect();
    let unloaded: Vec<usize> = (0..net.free.len())
        .filter(|&i| net.loads[i] == 0.0)
        .collect();
    let fixed = net.fixed_positions();
    let frange = |d: usize| {
        let c = fixed.column(d);
        (
            c.iter().cloned().fold(f64::INFINITY, f64::min),
            c.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        )
    };
    let fmt = |(lo, hi): (f64, f64)| format!("[{lo:.2}, {hi:.2}]");
    println!(
        "    anchors  x {} y {} z {}",
        fmt(frange(0)),
        fmt(frange(1)),
        fmt(frange(2))
    );
    for (label, set) in [("loaded  ", &loaded), ("unloaded", &unloaded)] {
        if set.is_empty() {
            continue;
        }
        println!(
            "    {label} ({:>4} nodes) x {} y {} z {}",
            set.len(),
            fmt(range(set, 0)),
            fmt(range(set, 1)),
            fmt(range(set, 2))
        );
    }
    let lengths = net.edge_lengths(x);
    let (shortest, _) =
        lengths.iter().enumerate().fold(
            (0, f64::INFINITY),
            |acc, (e, &l)| if l < acc.1 { (e, l) } else { acc },
        );
    let (a, b) = net.edges[shortest];
    println!(
        "    edge length [{:.3}, {:.3}] (shortest: edge {shortest} = {a}-{b}, plan {:?} {:?}, q {:.3})  |q| [{:.3}, {:.3}]  load [{:.4}, {:.4}]",
        lengths.iter().cloned().fold(f64::INFINITY, f64::min),
        lengths.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
        net.plan[a],
        net.plan[b],
        net.q_true[shortest],
        net.q_true.iter().map(|q| q.abs()).fold(f64::INFINITY, f64::min),
        net.q_true.iter().map(|q| q.abs()).fold(f64::NEG_INFINITY, f64::max),
        net.loads.iter().cloned().fold(f64::INFINITY, f64::min),
        net.loads.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    );
}

fn cmd_nets() {
    let verbose = std::env::var("BENCH_VERBOSE")
        .map(|v| v != "0")
        .unwrap_or(false);
    println!(
        "{:<18}{:>7}{:>7}{:>7}{:>6}{:>6}{:>10}{:>12}{:>12}  {}",
        "net", "ne", "nfree", "nfix", "q>0", "q<0", "depth/L", "|R|max", "|R|net", "status"
    );
    for net in selected_nets(all_nets(), None) {
        let pos = net.q_true.iter().filter(|&&q| q > 0.0).count();
        let neg = net.q_true.len() - pos;
        match net.try_funicular() {
            Some(x) => {
                let r = support_reactions(&net, &net.q_true, &x);
                let rmax = r
                    .rows()
                    .into_iter()
                    .map(|row| row.iter().map(|v| v * v).sum::<f64>().sqrt())
                    .fold(0.0, f64::max);
                let rn = net_reaction(&r);
                let rnet = (rn[0] * rn[0] + rn[1] * rn[1] + rn[2] * rn[2]).sqrt();
                let plan_extent = plan_span(&x);
                let status = if plan_extent > 3.0 * net.extent {
                    format!("OK (plan span {:.1}L)", plan_extent / net.extent)
                } else {
                    "OK".to_string()
                };
                println!(
                    "{:<18}{:>7}{:>7}{:>7}{:>6}{:>6}{:>10.4}{:>12.4e}{:>12.4e}  {}",
                    net.name,
                    net.edges.len(),
                    net.free.len(),
                    net.fixed.len(),
                    pos,
                    neg,
                    nets::depth_of(&x) / net.extent,
                    rmax,
                    rnet,
                    status
                );
                if verbose {
                    describe_shape(&net, &x);
                }
            }
            None => println!(
                "{:<18}{:>7}{:>7}{:>7}{:>6}{:>6}{:>10}{:>12}{:>12}  FAILED: forward solve at q_true not finite",
                net.name,
                net.edges.len(),
                net.free.len(),
                net.fixed.len(),
                pos,
                neg,
                "-",
                "-",
                "-"
            ),
        }
    }
}

/// Largest horizontal (x or y) span of the free nodes.
fn plan_span(x: &Array2<f64>) -> f64 {
    (0..2)
        .map(|d| {
            let c = x.column(d);
            c.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
                - c.iter().cloned().fold(f64::INFINITY, f64::min)
        })
        .fold(0.0, f64::max)
}

// ───────────────────────── suite / alt ─────────────────────────

fn cmd_suite(max_iters: usize, default_nets: Option<&[&str]>, default_methods: &[Method]) {
    let methods = selected_methods(default_methods);
    let nets = selected_nets(suite_nets(), default_nets);
    if nets.is_empty() || methods.is_empty() {
        eprintln!("nothing to run (check BENCH_NETS / BENCH_METHODS)");
        return;
    }
    println!(
        "# suite: max_iters={max_iters} methods=[{}] nets=[{}]",
        methods
            .iter()
            .map(|m| m.name())
            .collect::<Vec<_>>()
            .join(","),
        nets.iter()
            .map(|n| n.name.as_str())
            .collect::<Vec<_>>()
            .join(",")
    );
    for net in &nets {
        for (case, target) in targets(net) {
            for (box_name, lo, hi) in boxes(net) {
                let rows = run_case(net, &target, &lo, &hi, &methods, max_iters);
                print_case(
                    net,
                    &format!("{case} box={box_name} [{:.3},{:.3}]", lo[0], hi[0]),
                    &rows,
                    net.extent,
                );
            }
        }
    }
}

const ALT_METHODS: [Method; 7] = [
    Method::Uniform,
    Method::S1,
    Method::GramSparse,
    Method::LengthRatio,
    Method::Frozen,
    Method::Pipeline,
    Method::Legacy,
];
const ALT_NETS: [&str; 4] = ["quad21c", "hypar21m", "tiedarch16", "cabledome4x16"];

// ───────────────────────── scale ─────────────────────────

fn cmd_scale() {
    let methods = selected_methods(&Method::ALL);
    println!(
        "{:>5}{:>7}{:>7}  {:<18}{:>10}{:>11}{:>6}{:>7}{:>9}{:>11}  {}",
        "side",
        "ne",
        "nfree",
        "method",
        "warm_ms",
        "err_clip/L",
        "clip",
        "lb_it",
        "lb_ms",
        "final/L",
        "note"
    );
    let sides: Vec<usize> = std::env::var("BENCH_SIDES")
        .ok()
        .map(|v| v.split(',').filter_map(|s| s.trim().parse().ok()).collect())
        .unwrap_or_else(|| vec![23, 46, 91, 181]);
    for side in sides {
        let size_started = Instant::now();
        let net = quad(side, true, 1.0).normalized(0.25);
        let ne = net.edges.len();
        let target = net.target(0.02, 0.0);
        let (lo, hi) = net.box_from_true(1.0, 1.0);
        let fixed = net.fixed_positions();
        let problem = net.problem(&fixed, Vec::new(), &lo, &hi, SolverOptions::default());
        let prefix = format!("{:>5}{:>7}{:>7}  ", side, ne, net.free.len());
        println!(
            "{prefix}{:<18}{:>10.2}",
            "forward_solve",
            time_forward(&problem, &net.q_true)
        );
        println!(
            "{prefix}{:<18}{:>10.2}",
            "lbfgsb_eval",
            time_eval(&net, &fixed, &target, &net.q_true, &lo, &hi)
        );
        for &m in &methods {
            if size_started.elapsed().as_secs_f64() > SCALE_TIME_BUDGET_S {
                println!(
                    "{prefix}{:<18}  skipped (size exceeded {SCALE_TIME_BUDGET_S:.0}s budget)",
                    m.name()
                );
                continue;
            }
            let w = match run_method(m, &net, &problem, &target, &lo, &hi) {
                Ok(w) => w,
                Err(e) => {
                    println!(
                        "{prefix}{:<18}{:>10}{:>11}{:>6}{:>7}{:>9}{:>11}  {e}",
                        m.name(),
                        "-",
                        "-",
                        "-",
                        "-",
                        "-",
                        "-"
                    );
                    continue;
                }
            };
            let (q, n_clipped, _) = clip(&w.q, &lo, &hi);
            let err_clip = safe_err(&problem, &target, &q);
            let (lb_it, lb_ms, final_err, note) =
                if ne > SCALE_LBFGSB_EDGE_CAP || !err_clip.is_finite() {
                    (
                        "-".to_string(),
                        "-".to_string(),
                        f64::NAN,
                        join_note(&w.note, "lbfgsb skipped"),
                    )
                } else {
                    match lbfgsb(&net, &fixed, &target, &q, &lo, &hi, SCALE_LBFGSB_ITERS) {
                        Ok(run) => (
                            run.iters.to_string(),
                            format!("{:.0}", run.ms),
                            run.final_err,
                            w.note.clone(),
                        ),
                        Err(e) => (
                            "-".to_string(),
                            "-".to_string(),
                            f64::NAN,
                            join_note(&w.note, &e),
                        ),
                    }
                };
            println!(
                "{prefix}{:<18}{:>10.1}{:>11}{:>6}{:>7}{:>9}{:>11}  {note}",
                m.name(),
                w.ms,
                fmt_e(err_clip / net.extent),
                n_clipped,
                lb_it,
                lb_ms,
                fmt_e(final_err / net.extent)
            );
        }
        let elapsed = size_started.elapsed().as_secs_f64();
        println!("# side={side} ne={ne} took {elapsed:.1}s");
        if elapsed > SCALE_TIME_BUDGET_S {
            println!("# stopping: last size exceeded {SCALE_TIME_BUDGET_S:.0}s");
            break;
        }
    }
}

// ───────────────────────── dense ─────────────────────────

fn cmd_dense() {
    println!(
        "{:>5}{:>7}{:>10}{:>12}{:>12}{:>9}{:>9}{:>12}",
        "side", "ne", "bytes", "dense_ms", "sparse_ms", "ratio", "slope", "|Δq|/|q|"
    );
    let mut prev: Option<(f64, f64)> = None;
    for side in [12usize, 16, 23, 32, 45] {
        let net = quad(side, true, 1.0).normalized(0.25);
        let ne = net.edges.len();
        let target = net.target(0.02, 0.0);
        let (lo, hi) = net.box_from_true(1.0, 1.0);
        let problem = net.problem(
            &net.fixed_positions(),
            Vec::new(),
            &lo,
            &hi,
            SolverOptions::default(),
        );
        let sparse = match run_method(Method::GramSparse, &net, &problem, &target, &lo, &hi) {
            Ok(w) => w,
            Err(e) => {
                println!("{side:>5}{ne:>7}  gram_sparse failed: {e}");
                continue;
            }
        };
        if ne > DENSE_EDGE_CAP {
            println!(
                "{side:>5}{ne:>7}{:>10}{:>12}{:>12.1}  skipped (ne > cap {DENSE_EDGE_CAP})",
                "-", "-", sparse.ms
            );
            continue;
        }
        let (e, p, scale) = match equilibrium_e(&problem, &target) {
            Ok(v) => v,
            Err(err) => {
                println!("{side:>5}{ne:>7}  assemble failed: {err}");
                continue;
            }
        };
        let started = Instant::now();
        let dense = match dense_gram_solve(&e, &p, GRAM_RELATIVE_SHIFT * scale) {
            Ok(d) => d,
            Err(err) => {
                println!("{side:>5}{ne:>7}  dense failed: {err}");
                continue;
            }
        };
        let dense_ms = started.elapsed().as_secs_f64() * 1e3;
        let slope = prev
            .map(|(ne0, t0)| ((dense_ms / t0).ln() / (ne as f64 / ne0).ln()))
            .map(|s| format!("{s:.2}"))
            .unwrap_or_else(|| "-".into());
        let diff = dense
            .q
            .iter()
            .zip(&sparse.q)
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
            .sqrt()
            / sparse
                .q
                .iter()
                .map(|v| v * v)
                .sum::<f64>()
                .sqrt()
                .max(1e-300);
        println!(
            "{side:>5}{ne:>7}{:>10}{:>12.1}{:>12.1}{:>9.1}{:>9}{:>12.2e}",
            human_bytes(dense.bytes),
            dense_ms,
            sparse.ms,
            dense_ms / sparse.ms.max(1e-9),
            slope,
            diff
        );
        prev = Some((ne as f64, dense_ms));
    }
    println!("# slope = log-log exponent of dense_ms vs ne between consecutive sizes (O(ne³) → 3)");
}

// ───────────────────────── reactions ─────────────────────────

fn print_reactions(net: &Net, q: &[f64], x: &Array2<f64>) {
    let r = support_reactions(net, q, x);
    for (i, &n) in net.fixed.iter().enumerate() {
        println!(
            "    support node {n:>4} at ({:>7.2},{:>7.2},{:>7.2}):  R = ({:>11.4e}, {:>11.4e}, {:>11.4e})",
            net.plan[n][0], net.plan[n][1], net.plan[n][2], r[[i, 0]], r[[i, 1]], r[[i, 2]]
        );
    }
    let rn = net_reaction(&r);
    println!(
        "    net reaction = ({:.4e}, {:.4e}, {:.4e})   |horizontal| = {:.4e}",
        rn[0],
        rn[1],
        rn[2],
        (rn[0] * rn[0] + rn[1] * rn[1]).sqrt()
    );
    if !net.tie_edges.is_empty() {
        let tie: Vec<f64> = net.tie_edges.iter().map(|&e| q[e]).collect();
        let mean = tie.iter().sum::<f64>() / tie.len() as f64;
        println!(
            "    tie / bottom-chord q: mean {:.4}  min {:.4}  max {:.4}",
            mean,
            tie.iter().cloned().fold(f64::INFINITY, f64::min),
            tie.iter().cloned().fold(f64::NEG_INFINITY, f64::max)
        );
    }
}

fn cmd_reactions() {
    for net in [
        tied_arch_straight(16),
        tied_arch(16).normalized(0.25),
        cable_truss(16).normalized(0.25),
    ] {
        let (lo, hi) = net.box_from_true(100.0, 100.0);
        let fixed = net.fixed_positions();
        let problem = net.problem(&fixed, Vec::new(), &lo, &hi, SolverOptions::default());
        println!(
            "\n##### {} (ne={}, nfree={}, span along x)",
            net.name,
            net.edges.len(),
            net.free.len()
        );
        println!("  reference: q_true");
        print_reactions(&net, &net.q_true, &net.funicular());
        for (case, target) in [
            ("jit2%d", net.target(0.02, 0.0)),
            ("exact", net.funicular()),
        ] {
            for (label, rx, rz, weight, outer) in [
                ("rx free", false, false, 1.0, 2),
                ("rx=0", true, false, 1.0, 2),
                ("rx=0, 6 GN steps", true, false, 1.0, 6),
                ("rx=0, weight 10", true, false, 10.0, 2),
                ("rz=0", false, true, 1.0, 2),
            ] {
                let mut opts = library_options(Method::Pipeline, &lo, &hi, 0.0).unwrap();
                opts.enforce_zero_rx = rx;
                opts.enforce_zero_rz = rz;
                opts.reaction_weight = weight;
                opts.max_outer = outer;
                println!("  --- target {case}, pipeline, {label} ---");
                match run_library(&problem, &target, opts) {
                    Ok((r, ms)) => {
                        let (q, n_clipped, _) = clip(&r.q, &lo, &hi);
                        println!(
                            "    {:.1} ms, geometric error/L = {} (library {}), clipped {n_clipped}, {}",
                            ms,
                            fmt_e(safe_err(&problem, &target, &q) / net.extent),
                            fmt_e(r.geometric_error / net.extent),
                            diag_note(&r, true)
                        );
                        match try_forward(&problem, &q) {
                            Some(x) => print_reactions(&net, &q, &x),
                            None => println!("    forward solve at returned q failed"),
                        }
                    }
                    Err(e) => println!("    {e}"),
                }
            }
        }
    }
}

// ───────────────────────── main ─────────────────────────

fn usage() {
    eprintln!(
        "usage: warm_start_bench <nets|suite [max_iters]|alt [max_iters]|scale|dense|reactions>"
    );
    eprintln!("  env: BENCH_NETS=name1,name2  BENCH_METHODS=uniform,s1,...  BENCH_SIDES=23,46");
}

fn main() {
    std::panic::set_hook(Box::new(|info| {
        eprintln!("[caught panic] {}", info);
    }));
    let args: Vec<String> = std::env::args().collect();
    let max_iters = |i: usize| {
        args.get(i)
            .and_then(|s| s.parse().ok())
            .unwrap_or(DEFAULT_MAX_ITERS)
    };
    match args.get(1).map(String::as_str) {
        Some("nets") => cmd_nets(),
        Some("suite") => cmd_suite(max_iters(2), None, &Method::ALL),
        Some("alt") => cmd_suite(max_iters(2), Some(&ALT_NETS), &ALT_METHODS),
        Some("scale") => cmd_scale(),
        Some("dense") => cmd_dense(),
        Some("reactions") => cmd_reactions(),
        _ => usage(),
    }
}
