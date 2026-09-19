//! Gauss–Newton / frozen CWLS versus short L-BFGS-B on the exact geometric
//! objective, starting from the Stage-1 force-residual particular `q*`.
//!
//! The question is whether a few compliance-weighted linearisations (each a
//! sparse saddle LDL) buy more geometric progress than handing `q*` to
//! L-BFGS-B and letting it spend the same budget on the nonlinear `‖x(q)−x*‖²`
//! — in particular after 10 accepted L-BFGS-B steps, which saturates the
//! default compact-history size `m = 10`.

use crate::methods::{base_options, clip, equilibrium_e, lbfgsb, run_library};
use crate::nets::{safe_err, suite_nets, try_forward, Net};
use crate::report::fmt_e;
use ndarray::Array2;
use serde_json::{json, Value};
use std::path::Path;
use std::time::Instant;
use theseus::gradients::value_and_gradient;
use theseus::inverse::{InverseDiagnostics, InverseFdmOptions, InverseMetric};
use theseus::sparse::SparseColMatOwned;
use theseus::types::{
    FdmCache, ObjectiveTrait, Problem, SolverOptions, TargetGeometryReduction, TargetXYZ,
};

/// Cases: easy hanging net plus the presentation showcase (mixed-sign and
/// collapsed Stage-1 seeds). Targets/boxes match `figures.rs`.
const CASES: &[(&str, &str, &str)] = &[
    ("quad21c", "jit2%d", "loose"),
    ("cabledome4x16", "jit2%d", "loose"),
    ("cabledome4x16", "bump10%d", "snug"),
    ("hypar21m", "jit2%d", "loose"),
    ("cabletruss16", "jit2%d", "loose"),
    ("quad21c_d1", "jit2%d", "loose"),
    ("holes21c", "jit2%d", "loose"),
    ("crease21", "jit2%d", "loose"),
    ("barrel16x12", "jit2%d", "loose"),
    ("wheel24a4", "jit2%d", "loose"),
    ("tiedarch16", "jit2%d", "loose"),
    ("oculus21h7", "bump10%d", "loose"),
];

const LBFGS_BUDGETS: [usize; 4] = [1, 3, 10, 20];

#[derive(Clone)]
struct Residuals {
    geom: f64,
    geom_max: f64,
    geom_rms: f64,
    force: f64,
    force_rel: f64,
    sse: f64,
    grad_inf: f64,
    pgrad_inf: f64,
    n_active: usize,
    finite: bool,
}

#[derive(Clone)]
struct Scored {
    label: String,
    family: String,
    ms: f64,
    residuals: Residuals,
    note: String,
    frozen_steps: usize,
    newton_steps: usize,
    factorizations: usize,
    lb_iters: usize,
    lb_evals: usize,
    q: Vec<f64>,
    trace_geom: Vec<f64>,
}

fn target_for(net: &Net, kind: &str) -> Array2<f64> {
    match kind {
        "bump10%d" => net.target(0.0, 0.10),
        "exact" => net.funicular(),
        _ => net.target(0.02, 0.0),
    }
}

fn box_for(net: &Net, kind: &str) -> (Vec<f64>, Vec<f64>) {
    match kind {
        "snug" => net.box_from_true(1.0, 1.0),
        _ => net.box_from_true(100.0, 100.0),
    }
}

fn l2(v: &[f64]) -> f64 {
    v.iter().map(|x| x * x).sum::<f64>().sqrt()
}

fn force_residual(e: &SparseColMatOwned, p: &[f64], q: &[f64]) -> (f64, f64) {
    if q.len() != e.ncols {
        return (f64::NAN, f64::NAN);
    }
    let mut r = e.matvec(q);
    for (ri, &pi) in r.iter_mut().zip(p.iter()) {
        *ri -= pi;
    }
    let force = l2(&r);
    let pnorm = l2(p);
    let rel = if pnorm > 0.0 { force / pnorm } else { force };
    (force, rel)
}

fn target_objective(net: &Net, target: &Array2<f64>) -> Box<dyn ObjectiveTrait> {
    Box::new(TargetXYZ {
        weight: 1.0,
        node_indices: net.free.clone(),
        target: target.clone(),
        reduction: TargetGeometryReduction::Sse,
    })
}

fn projected_grad_inf(q: &[f64], grad: &[f64], lo: &[f64], hi: &[f64]) -> f64 {
    q.iter()
        .zip(grad)
        .zip(lo)
        .zip(hi)
        .map(|(((&qi, &gi), &a), &b)| {
            if qi <= a && gi > 0.0 {
                0.0
            } else if qi >= b && gi < 0.0 {
                0.0
            } else {
                gi.abs()
            }
        })
        .fold(0.0, f64::max)
}

fn score(
    problem: &Problem,
    obj_problem: &Problem,
    target: &Array2<f64>,
    e: &SparseColMatOwned,
    p: &[f64],
    q: &[f64],
    lo: &[f64],
    hi: &[f64],
) -> Residuals {
    let (q, _, n_active) = clip(q, lo, hi);
    let (force, force_rel) = force_residual(e, p, &q);
    let x = try_forward(problem, &q);
    let (geom, geom_max, geom_rms, sse, finite) = match &x {
        Some(x) => {
            let n = x.nrows().max(1) as f64;
            let mut sse = 0.0;
            let mut max = 0.0;
            for i in 0..x.nrows() {
                let mut d2 = 0.0;
                for d in 0..3 {
                    let diff = x[[i, d]] - target[[i, d]];
                    d2 += diff * diff;
                }
                sse += d2;
                max = f64::max(max, d2.sqrt());
            }
            (sse.sqrt(), max, (sse / n).sqrt(), sse, true)
        }
        None => (f64::NAN, f64::NAN, f64::NAN, f64::NAN, false),
    };
    let (grad_inf, pgrad_inf) = if finite {
        let mut cache = match FdmCache::new(obj_problem) {
            Ok(c) => c,
            Err(_) => {
                return Residuals {
                    geom,
                    geom_max,
                    geom_rms,
                    force,
                    force_rel,
                    sse,
                    grad_inf: f64::NAN,
                    pgrad_inf: f64::NAN,
                    n_active,
                    finite,
                }
            }
        };
        let mut grad = vec![0.0; q.len()];
        match value_and_gradient(&mut cache, obj_problem, &q, &mut grad, lo, hi, &[], &[]) {
            Ok(_) => (
                grad.iter().fold(0.0_f64, |a, g| f64::max(a, g.abs())),
                projected_grad_inf(&q, &grad, lo, hi),
            ),
            Err(_) => (f64::NAN, f64::NAN),
        }
    } else {
        (f64::NAN, f64::NAN)
    };
    Residuals {
        geom,
        geom_max,
        geom_rms,
        force,
        force_rel,
        sse,
        grad_inf,
        pgrad_inf,
        n_active,
        finite,
    }
}

fn inverse_opts(lo: &[f64], hi: &[f64], frozen: usize, newton: usize, guard: f64) -> InverseFdmOptions {
    let mut opts = base_options(lo, hi);
    opts.seed_guard_margin = guard;
    if newton == 0 {
        opts.metric = InverseMetric::Geometry;
        opts.max_frozen_outer = 0;
        opts.max_outer = frozen;
    } else {
        opts.metric = InverseMetric::GeometryNewton;
        opts.max_frozen_outer = frozen;
        opts.max_outer = newton;
    }
    opts
}

fn diag_bits(d: &InverseDiagnostics) -> String {
    format!(
        "s1e={:.2e} fr{} gn{} fac={}{}",
        d.stage1_error,
        d.frozen_steps,
        d.newton_steps,
        d.stage2_factorizations,
        if d.used_uniform_seed { " guard→u" } else { "" }
    )
}

fn run_inverse(
    _net: &Net,
    problem: &Problem,
    obj_problem: &Problem,
    target: &Array2<f64>,
    e: &SparseColMatOwned,
    p: &[f64],
    lo: &[f64],
    hi: &[f64],
    label: &str,
    family: &str,
    opts: InverseFdmOptions,
) -> Scored {
    let stage2 = opts.metric.is_geometric();
    match run_library(problem, target, opts) {
        Ok((r, ms)) => {
            let (q, _, _) = clip(&r.q, lo, hi);
            let residuals = score(problem, obj_problem, target, e, p, &q, lo, hi);
            Scored {
                label: label.to_string(),
                family: family.to_string(),
                ms,
                residuals,
                note: if stage2 {
                    diag_bits(&r.diagnostics)
                } else {
                    format!("lib_err={:.2e}", r.geometric_error)
                },
                frozen_steps: r.diagnostics.frozen_steps,
                newton_steps: r.diagnostics.newton_steps,
                factorizations: r.diagnostics.stage2_factorizations,
                lb_iters: 0,
                lb_evals: 0,
                q,
                trace_geom: Vec::new(),
            }
        }
        Err(err) => Scored {
            label: label.to_string(),
            family: family.to_string(),
            ms: 0.0,
            residuals: Residuals {
                geom: f64::NAN,
                geom_max: f64::NAN,
                geom_rms: f64::NAN,
                force: f64::NAN,
                force_rel: f64::NAN,
                sse: f64::NAN,
                grad_inf: f64::NAN,
                pgrad_inf: f64::NAN,
                n_active: 0,
                finite: false,
            },
            note: err,
            frozen_steps: 0,
            newton_steps: 0,
            factorizations: 0,
            lb_iters: 0,
            lb_evals: 0,
            q: Vec::new(),
            trace_geom: Vec::new(),
        },
    }
}

fn run_lbfgs_from(
    net: &Net,
    problem: &Problem,
    obj_problem: &Problem,
    fixed: &Array2<f64>,
    target: &Array2<f64>,
    e: &SparseColMatOwned,
    p: &[f64],
    lo: &[f64],
    hi: &[f64],
    seed_label: &str,
    q0: &[f64],
    max_iters: usize,
) -> Scored {
    let label = format!("{seed_label}+lb{max_iters}");
    if q0.is_empty() {
        return Scored {
            label,
            family: format!("lbfgs:{seed_label}"),
            ms: 0.0,
            residuals: Residuals {
                geom: f64::NAN,
                geom_max: f64::NAN,
                geom_rms: f64::NAN,
                force: f64::NAN,
                force_rel: f64::NAN,
                sse: f64::NAN,
                grad_inf: f64::NAN,
                pgrad_inf: f64::NAN,
                n_active: 0,
                finite: false,
            },
            note: format!("no {seed_label} seed"),
            frozen_steps: 0,
            newton_steps: 0,
            factorizations: 0,
            lb_iters: 0,
            lb_evals: 0,
            q: Vec::new(),
            trace_geom: Vec::new(),
        };
    }
    match lbfgsb(net, fixed, target, q0, lo, hi, max_iters) {
        Ok(run) => {
            let residuals = score(problem, obj_problem, target, e, p, &run.q, lo, hi);
            let trace_geom: Vec<f64> = run
                .trace
                .iter()
                .map(|&sse| if sse >= 0.0 { sse.sqrt() } else { f64::NAN })
                .collect();
            Scored {
                label,
                family: format!("lbfgs:{seed_label}"),
                ms: run.ms,
                residuals,
                note: format!("accepted={} evals≈{}", run.iters, run.trace.len()),
                frozen_steps: 0,
                newton_steps: 0,
                factorizations: 0,
                lb_iters: run.iters,
                lb_evals: run.trace.len(),
                q: run.q,
                trace_geom,
            }
        }
        Err(err) => Scored {
            label,
            family: format!("lbfgs:{seed_label}"),
            ms: 0.0,
            residuals: Residuals {
                geom: f64::NAN,
                geom_max: f64::NAN,
                geom_rms: f64::NAN,
                force: f64::NAN,
                force_rel: f64::NAN,
                sse: f64::NAN,
                grad_inf: f64::NAN,
                pgrad_inf: f64::NAN,
                n_active: 0,
                finite: false,
            },
            note: err,
            frozen_steps: 0,
            newton_steps: 0,
            factorizations: 0,
            lb_iters: 0,
            lb_evals: 0,
            q: Vec::new(),
            trace_geom: Vec::new(),
        },
    }
}

fn cosine(a: &[f64], b: &[f64]) -> f64 {
    if a.len() != b.len() || a.is_empty() {
        return f64::NAN;
    }
    let mut dot = 0.0;
    let mut na = 0.0;
    let mut nb = 0.0;
    for i in 0..a.len() {
        let da = a[i] - 0.0;
        let db = b[i] - 0.0;
        // callers pass deltas already
        let _ = (da, db);
        dot += a[i] * b[i];
        na += a[i] * a[i];
        nb += b[i] * b[i];
    }
    if na <= 0.0 || nb <= 0.0 {
        return f64::NAN;
    }
    (dot / (na.sqrt() * nb.sqrt())).clamp(-1.0, 1.0)
}

fn delta(a: &[f64], b: &[f64]) -> Vec<f64> {
    if a.len() != b.len() {
        return Vec::new();
    }
    a.iter().zip(b).map(|(x, y)| x - y).collect()
}

fn print_residual_header() {
    println!(
        "{:<16}{:>8}{:>6}{:>4}{:>4}{:>6}{:>6}{:>11}{:>10}{:>10}{:>11}{:>10}{:>10}{:>8}  {}",
        "method",
        "ms",
        "fac",
        "fr",
        "gn",
        "lb_it",
        "lb_ev",
        "geom/L",
        "max/L",
        "rms/L",
        "||r||",
        "||r||/||p||",
        "||∇f||∞",
        "active",
        "note"
    );
}

fn print_scored(s: &Scored, l: f64) {
    let r = &s.residuals;
    println!(
        "{:<16}{:>8.1}{:>6}{:>4}{:>4}{:>6}{:>6}{:>11}{:>10}{:>10}{:>11}{:>10}{:>10}{:>8}  {}",
        s.label,
        s.ms,
        if s.factorizations > 0 {
            s.factorizations.to_string()
        } else {
            "-".into()
        },
        if s.family.starts_with("lbfgs") {
            "-".into()
        } else {
            s.frozen_steps.to_string()
        },
        if s.family.starts_with("lbfgs") {
            "-".into()
        } else {
            s.newton_steps.to_string()
        },
        if s.lb_iters > 0 {
            s.lb_iters.to_string()
        } else {
            "-".into()
        },
        if s.lb_evals > 0 {
            s.lb_evals.to_string()
        } else {
            "-".into()
        },
        fmt_e(r.geom / l),
        fmt_e(r.geom_max / l),
        fmt_e(r.geom_rms / l),
        fmt_e(r.force),
        fmt_e(r.force_rel),
        fmt_e(r.pgrad_inf),
        r.n_active,
        s.note
    );
}

fn scored_json(s: &Scored, l: f64) -> Value {
    let r = &s.residuals;
    json!({
        "label": s.label,
        "family": s.family,
        "ms": s.ms,
        "frozen_steps": s.frozen_steps,
        "newton_steps": s.newton_steps,
        "factorizations": s.factorizations,
        "lb_iters": s.lb_iters,
        "lb_evals": s.lb_evals,
        "geom": nan_none(r.geom),
        "geom_over_l": nan_none(r.geom / l),
        "geom_max_over_l": nan_none(r.geom_max / l),
        "geom_rms_over_l": nan_none(r.geom_rms / l),
        "force": nan_none(r.force),
        "force_rel": nan_none(r.force_rel),
        "sse": nan_none(r.sse),
        "grad_inf": nan_none(r.grad_inf),
        "pgrad_inf": nan_none(r.pgrad_inf),
        "n_active": r.n_active,
        "finite": r.finite,
        "note": s.note,
        "trace_geom": s.trace_geom.iter().map(|&v| nan_none(v)).collect::<Vec<_>>(),
        "trace_geom_over_l": s.trace_geom.iter().map(|&v| nan_none(v / l)).collect::<Vec<_>>(),
    })
}

fn nan_none(v: f64) -> Option<f64> {
    v.is_finite().then_some(v)
}

fn find_row<'a>(rows: &'a [Scored], label: &str) -> Option<&'a Scored> {
    rows.iter().find(|s| s.label == label)
}

fn geom_l(rows: &[Scored], label: &str, l: f64) -> f64 {
    find_row(rows, label)
        .map(|s| s.residuals.geom / l)
        .unwrap_or(f64::NAN)
}

fn run_case(net: &Net, target_kind: &str, box_kind: &str) -> Value {
    let target = target_for(net, target_kind);
    let (lo, hi) = box_for(net, box_kind);
    let fixed = net.fixed_positions();
    let problem = net.problem(&fixed, Vec::new(), &lo, &hi, SolverOptions::default());
    let obj_problem = net.problem(
        &fixed,
        vec![target_objective(net, &target)],
        &lo,
        &hi,
        SolverOptions {
            absolute_tolerance: 1e-8,
            relative_tolerance: 1e-12,
            q_parameterization_mode: theseus::types::QParameterizationMode::DirectBoxBounds,
            ..SolverOptions::default()
        },
    );
    let (e, p, _) = match equilibrium_e(&problem, &target) {
        Ok(v) => v,
        Err(err) => {
            println!("  failed to assemble E(x*): {err}");
            return json!({"net": net.name, "error": err});
        }
    };
    let l = net.extent;
    println!(
        "\n=== {} {} box={} (ne={}, nfree={})  L={:.4} ===",
        net.name,
        target_kind,
        box_kind,
        net.edges.len(),
        net.free.len(),
        l
    );

    let mut rows = Vec::new();

    // Stage-1 force-residual particular (no Stage-2).
    let mut s1_opts = base_options(&lo, &hi);
    s1_opts.metric = InverseMetric::Force;
    s1_opts.seed_guard_margin = 0.0;
    rows.push(run_inverse(
        net, &problem, &obj_problem, &target, &e, &p, &lo, &hi, "s1", "inverse", s1_opts,
    ));

    // One / two frozen CWLS steps from q*, guard off.
    rows.push(run_inverse(
        net,
        &problem,
        &obj_problem,
        &target,
        &e,
        &p,
        &lo,
        &hi,
        "frozen1",
        "inverse",
        inverse_opts(&lo, &hi, 1, 0, 0.0),
    ));
    rows.push(run_inverse(
        net,
        &problem,
        &obj_problem,
        &target,
        &e,
        &p,
        &lo,
        &hi,
        "frozen2",
        "inverse",
        inverse_opts(&lo, &hi, 2, 0, 0.0),
    ));

    // Pure Gauss–Newton from q* (no frozen linearisation at x*).
    for n in [1usize, 2, 3] {
        rows.push(run_inverse(
            net,
            &problem,
            &obj_problem,
            &target,
            &e,
            &p,
            &lo,
            &hi,
            &format!("gn{n}"),
            "inverse",
            inverse_opts(&lo, &hi, 0, n, 0.0),
        ));
    }

    // Frozen then GN, guard off (controlled) and on (production).
    rows.push(run_inverse(
        net,
        &problem,
        &obj_problem,
        &target,
        &e,
        &p,
        &lo,
        &hi,
        "pipe",
        "inverse",
        inverse_opts(&lo, &hi, 1, 2, 0.0),
    ));
    rows.push(run_inverse(
        net,
        &problem,
        &obj_problem,
        &target,
        &e,
        &p,
        &lo,
        &hi,
        "pipe_guard",
        "inverse",
        inverse_opts(&lo, &hi, 1, 2, 3.0),
    ));

    let q_s1 = find_row(&rows, "s1").map(|s| s.q.clone()).unwrap_or_default();
    let q_fr = find_row(&rows, "frozen1")
        .map(|s| s.q.clone())
        .unwrap_or_default();
    let q_gn1 = find_row(&rows, "gn1").map(|s| s.q.clone()).unwrap_or_default();

    // L-BFGS-B from each seed at the requested accepted-iteration budgets.
    for (seed, q0) in [("s1", q_s1.as_slice()), ("frozen1", q_fr.as_slice()), ("gn1", q_gn1.as_slice())]
    {
        if q0.is_empty() {
            continue;
        }
        for &k in &LBFGS_BUDGETS {
            rows.push(run_lbfgs_from(
                net,
                &problem,
                &obj_problem,
                &fixed,
                &target,
                &e,
                &p,
                &lo,
                &hi,
                seed,
                q0,
                k,
            ));
        }
    }

    // Direction comparison: frozen CWLS step vs first L-BFGS step vs first GN step.
    if let (Some(s1), Some(fr), Some(gn1), Some(lb1)) = (
        find_row(&rows, "s1"),
        find_row(&rows, "frozen1"),
        find_row(&rows, "gn1"),
        find_row(&rows, "s1+lb1"),
    ) {
        if !s1.q.is_empty() && !fr.q.is_empty() && !gn1.q.is_empty() && !lb1.q.is_empty() {
            let d_fr = delta(&fr.q, &s1.q);
            let d_gn = delta(&gn1.q, &s1.q);
            let d_lb = delta(&lb1.q, &s1.q);
            let cos_fr_lb = cosine(&d_fr, &d_lb);
            let cos_fr_gn = cosine(&d_fr, &d_gn);
            let cos_gn_lb = cosine(&d_gn, &d_lb);
            println!(
                "  q-space step cosine from s1:  <frozen, lb1>={}  <frozen, gn1>={}  <gn1, lb1>={}   ||Δfr||={}  ||Δgn||={}  ||Δlb1||={}",
                fmt_e(cos_fr_lb),
                fmt_e(cos_fr_gn),
                fmt_e(cos_gn_lb),
                fmt_e(l2(&d_fr)),
                fmt_e(l2(&d_gn)),
                fmt_e(l2(&d_lb))
            );
        }
    }

    println!("-- linearised CWLS / Gauss–Newton from q* (seed guard off except pipe_guard) --");
    print_residual_header();
    for s in rows.iter().filter(|s| s.family == "inverse") {
        print_scored(s, l);
    }
    println!("-- L-BFGS-B on ‖x(q)−x*‖² from those seeds (history m=10) --");
    print_residual_header();
    for s in rows.iter().filter(|s| s.family.starts_with("lbfgs")) {
        print_scored(s, l);
    }

    // Per-evaluation traces for the 10-step runs (history saturation).
    for label in ["s1+lb10", "frozen1+lb10", "gn1+lb10"] {
        if let Some(s) = find_row(&rows, label) {
            if !s.trace_geom.is_empty() {
                let vals: Vec<String> = s
                    .trace_geom
                    .iter()
                    .map(|&v| fmt_e(v / l))
                    .collect();
                println!("  trace {label} geom/L per eval: {}", vals.join(" "));
            }
        }
    }

    // Identity check: geometric error from the forward solve vs library probe.
    if let Some(s1) = find_row(&rows, "s1") {
        let fwd = safe_err(&problem, &target, &s1.q);
        println!(
            "  identity check s1: forward ‖x−x*‖={}  scored geom={}  (should match)",
            fmt_e(fwd),
            fmt_e(s1.residuals.geom)
        );
    }

    json!({
        "net": net.name,
        "target": target_kind,
        "box": box_kind,
        "ne": net.edges.len(),
        "nfree": net.free.len(),
        "extent": l,
        "rows": rows.iter().map(|s| scored_json(s, l)).collect::<Vec<_>>(),
        "summary": {
            "s1": nan_none(geom_l(&rows, "s1", l)),
            "frozen1": nan_none(geom_l(&rows, "frozen1", l)),
            "frozen2": nan_none(geom_l(&rows, "frozen2", l)),
            "gn1": nan_none(geom_l(&rows, "gn1", l)),
            "gn2": nan_none(geom_l(&rows, "gn2", l)),
            "gn3": nan_none(geom_l(&rows, "gn3", l)),
            "pipe": nan_none(geom_l(&rows, "pipe", l)),
            "s1_lb10": nan_none(geom_l(&rows, "s1+lb10", l)),
            "s1_lb20": nan_none(geom_l(&rows, "s1+lb20", l)),
            "frozen1_lb10": nan_none(geom_l(&rows, "frozen1+lb10", l)),
            "gn1_lb10": nan_none(geom_l(&rows, "gn1+lb10", l)),
        }
    })
}

fn selected_cases() -> Vec<(Net, &'static str, &'static str)> {
    let filter = std::env::var("BENCH_NETS").ok().filter(|s| !s.trim().is_empty());
    let nets = suite_nets();
    let mut out = Vec::new();
    for &(name, target, box_kind) in CASES {
        if let Some(ref f) = filter {
            if !f.split(',').any(|n| n.trim() == name) {
                continue;
            }
        }
        if let Some(net) = nets.iter().find(|n| n.name == name) {
            out.push((net.clone(), target, box_kind));
        } else {
            eprintln!("unknown net '{name}' in tradeoff case list");
        }
    }
    out
}

fn print_cross_case_summary(cases: &[Value]) {
    println!("\n===== cross-case geometric error / L =====");
    println!(
        "{:<28}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>10}{:>10}{:>10}",
        "case", "s1", "frozen1", "gn2", "pipe", "s1+lb10", "fr+lb10", "s1+lb10/fr", "s1+lb10/gn2", "s1+lb10/pipe"
    );
    for c in cases {
        let Some(s) = c.get("summary") else { continue };
        let g = |k: &str| {
            s.get(k)
                .and_then(|v| v.as_f64())
                .unwrap_or(f64::NAN)
        };
        let ratio = |a: f64, b: f64| {
            if a.is_finite() && b.is_finite() && b > 0.0 {
                format!("{:.2}", a / b)
            } else {
                "-".into()
            }
        };
        let case = format!(
            "{}_{}_{}",
            c.get("net").and_then(|v| v.as_str()).unwrap_or("?"),
            c.get("target").and_then(|v| v.as_str()).unwrap_or("?"),
            c.get("box").and_then(|v| v.as_str()).unwrap_or("?")
        );
        println!(
            "{:<28}{:>9}{:>9}{:>9}{:>9}{:>9}{:>9}{:>10}{:>10}{:>10}",
            case,
            fmt_e(g("s1")),
            fmt_e(g("frozen1")),
            fmt_e(g("gn2")),
            fmt_e(g("pipe")),
            fmt_e(g("s1_lb10")),
            fmt_e(g("frozen1_lb10")),
            ratio(g("s1_lb10"), g("frozen1")),
            ratio(g("s1_lb10"), g("gn2")),
            ratio(g("s1_lb10"), g("pipe")),
        );
    }
}

/// Run the GN-versus-L-BFGS-B residual experiment.
///
/// `out_dir` writes `tradeoff.json` (and is created if needed). Tables always
/// go to stdout.
pub fn cmd_tradeoff(out_dir: Option<&str>) {
    let cases = selected_cases();
    if cases.is_empty() {
        eprintln!("nothing to run (check BENCH_NETS)");
        return;
    }
    println!(
        "# tradeoff: frozen/GN linearisations vs L-BFGS-B on ‖x(q)−x*‖² from q*"
    );
    println!(
        "# history m=10; L-BFGS budgets {:?}; seed guard off except pipe_guard",
        LBFGS_BUDGETS
    );
    println!(
        "# residuals: geom=‖x−x*‖  max=max nodal  rms=RMS nodal  r=E(x*)q−p  ∇f of SSE"
    );
    let started = Instant::now();
    let mut json_cases = Vec::new();
    for (net, target, box_kind) in &cases {
        json_cases.push(run_case(net, target, box_kind));
    }
    print_cross_case_summary(&json_cases);
    println!(
        "\n# finished {} cases in {:.1}s",
        json_cases.len(),
        started.elapsed().as_secs_f64()
    );
    if let Some(dir) = out_dir {
        let path = Path::new(dir);
        std::fs::create_dir_all(path).unwrap_or_else(|e| panic!("mkdir {}: {e}", path.display()));
        let payload = json!({
            "title": "Gauss-Newton / frozen CWLS vs short L-BFGS-B from Stage-1 q*",
            "lbfgs_history": 10,
            "lbfgs_budgets": LBFGS_BUDGETS,
            "cases": json_cases,
        });
        let file = path.join("tradeoff.json");
        std::fs::write(&file, serde_json::to_string_pretty(&payload).unwrap())
            .unwrap_or_else(|e| panic!("write {}: {e}", file.display()));
        println!("wrote {}", file.display());
    }
}
