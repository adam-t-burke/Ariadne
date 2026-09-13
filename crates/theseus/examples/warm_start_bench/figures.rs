//! Data export for the presentation figures (`figures <out_dir>`).
//!
//! Writes one JSON file per showcase case with the full geometry of the
//! target, of every method's clipped warm start and of its L-BFGS-B result,
//! together with the loss traces, plus a gallery file with the funicular of
//! every suite net and a reactions file for the straight-tie arch. The
//! plots themselves are produced by `bench/figures/render.py`.

use crate::methods::{clip, diag_note, lbfgsb, library_options, run_library, run_method, Method};
use crate::nets::{
    net_reaction, safe_err, suite_nets, support_reactions, tied_arch_straight, try_forward, Net,
};
use ndarray::Array2;
use serde_json::{json, Value};
use std::path::Path;
use theseus::types::SolverOptions;

/// Cases shown in the talk: (net, target, box). Targets are the same as the
/// suite (`jit2%d` = 2 % white noise on the depth, `bump10%d` = smooth 10 %
/// push), boxes are `loose` (×/÷100) or `snug` (×/÷1 around `q_true`).
const SHOWCASE: &[(&str, &str, &str)] = &[
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

const FIGURE_METHODS: [Method; 7] = [
    Method::Uniform,
    Method::S1,
    Method::GramSparse,
    Method::LengthRatio,
    Method::Frozen,
    Method::Pipeline,
    Method::Legacy,
];

fn array_json(x: &Array2<f64>) -> Value {
    Value::Array(
        x.rows()
            .into_iter()
            .map(|r| json!([r[0], r[1], r[2]]))
            .collect(),
    )
}

fn nan_none(v: f64) -> Option<f64> {
    v.is_finite().then_some(v)
}

fn finite_vec(v: &[f64]) -> Vec<Option<f64>> {
    v.iter().map(|&x| nan_none(x)).collect()
}

/// Full node positions at `q`, or `None` when the forward solve fails.
fn positions_at(net: &Net, problem: &theseus::types::Problem, q: &[f64]) -> Option<Value> {
    try_forward(problem, q).map(|x| array_json(&net.full_positions(&x)))
}

fn net_json(net: &Net) -> Value {
    json!({
        "name": net.name,
        "edges": net.edges,
        "free": net.free,
        "fixed": net.fixed,
        "extent": net.extent,
        "q_true": net.q_true,
        "tie_edges": net.tie_edges,
        "x_true": array_json(&net.full_positions(&net.funicular())),
    })
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

fn write(path: &Path, value: &Value) {
    let text = serde_json::to_string(value).expect("serialise");
    std::fs::write(path, text).unwrap_or_else(|e| panic!("write {}: {e}", path.display()));
    println!("wrote {}", path.display());
}

fn export_case(net: &Net, target_kind: &str, box_kind: &str, max_iters: usize, out: &Path) {
    let target = target_for(net, target_kind);
    let (lo, hi) = box_for(net, box_kind);
    let fixed = net.fixed_positions();
    let problem = net.problem(&fixed, Vec::new(), &lo, &hi, SolverOptions::default());
    let mut methods = Vec::new();
    for &m in &FIGURE_METHODS {
        let label = m.name();
        let w = match run_method(m, net, &problem, &target, &lo, &hi) {
            Ok(w) => w,
            Err(e) => {
                println!("  {label:<14} warm start failed: {e}");
                methods.push(json!({ "label": label, "failed": e }));
                continue;
            }
        };
        let err_raw = safe_err(&problem, &target, &w.q);
        let (q, n_clipped, _) = clip(&w.q, &lo, &hi);
        let err_clip = safe_err(&problem, &target, &q);
        let mut entry = json!({
            "label": label,
            "warm_ms": w.ms,
            "note": w.note,
            "err_raw": nan_none(err_raw),
            "err_clip": nan_none(err_clip),
            "n_clipped": n_clipped,
            "q_warm": finite_vec(&q),
            "x_warm": positions_at(net, &problem, &q),
        });
        if err_clip.is_finite() {
            match lbfgsb(net, &fixed, &target, &q, &lo, &hi, max_iters) {
                Ok(run) => {
                    entry["trace"] = json!(run.trace);
                    entry["lbfgsb_iters"] = json!(run.iters);
                    entry["lbfgsb_ms"] = json!(run.ms);
                    entry["final_err"] = json!(nan_none(run.final_err));
                    entry["x_final"] = positions_at(net, &problem, &run.q).unwrap_or(Value::Null);
                    entry["q_final"] = json!(finite_vec(&run.q));
                    println!(
                        "  {label:<14} start/L {:.3e}  final/L {:.3e}  ({} it)",
                        err_clip / net.extent,
                        run.final_err / net.extent,
                        run.iters
                    );
                }
                Err(e) => {
                    entry["failed"] = json!(e);
                    println!("  {label:<14} L-BFGS-B failed: {e}");
                }
            }
        } else {
            entry["failed"] = json!("clipped q gives singular forward solve");
            println!("  {label:<14} singular forward solve after clipping");
        }
        methods.push(entry);
    }
    let value = json!({
        "net": net_json(net),
        "target_kind": target_kind,
        "box_kind": box_kind,
        "lo": lo,
        "hi": hi,
        "max_iters": max_iters,
        "x_target": array_json(&net.full_positions(&target)),
        "methods": methods,
    });
    let file = format!(
        "case_{}_{}_{}.json",
        net.name,
        target_kind.replace('%', "pct"),
        box_kind
    );
    write(&out.join(file), &value);
}

fn export_gallery(out: &Path) {
    let nets: Vec<Value> = suite_nets().iter().map(net_json).collect();
    write(&out.join("gallery.json"), &Value::Array(nets));
}

/// Straight-tie arch, exact and jittered targets, with and without the
/// zero-horizontal-reaction rows: geometry, `q`, and per-support reactions.
fn export_reactions(out: &Path) {
    let net = tied_arch_straight(16);
    let (lo, hi) = net.box_from_true(100.0, 100.0);
    let fixed = net.fixed_positions();
    let problem = net.problem(&fixed, Vec::new(), &lo, &hi, SolverOptions::default());
    let reactions_json = |q: &[f64], x: &Array2<f64>| -> Value {
        let r = support_reactions(&net, q, x);
        let rn = net_reaction(&r);
        json!({ "per_support": array_json(&r), "net": [rn[0], rn[1], rn[2]] })
    };
    let x_true = net.funicular();
    let mut runs = vec![json!({
        "label": "reference q_true",
        "target_kind": "exact",
        "x": array_json(&net.full_positions(&x_true)),
        "q": net.q_true,
        "reactions": reactions_json(&net.q_true, &x_true),
    })];
    for target_kind in ["exact", "jit2%d"] {
        let target = target_for(&net, target_kind);
        for (label, rx, weight) in [("rx free", false, 1.0), ("rx=0", true, 1.0)] {
            let mut opts = library_options(Method::Pipeline, &lo, &hi, 0.0).unwrap();
            opts.enforce_zero_rx = rx;
            opts.reaction_weight = weight;
            opts.max_outer = 2;
            match run_library(&problem, &target, opts) {
                Ok((r, ms)) => {
                    let (q, _, _) = clip(&r.q, &lo, &hi);
                    let err = safe_err(&problem, &target, &q);
                    let x = try_forward(&problem, &q);
                    println!(
                        "  reactions {target_kind:<7} {label:<8} err/L {:.3e}  {}",
                        err / net.extent,
                        diag_note(&r, true)
                    );
                    runs.push(json!({
                        "label": label,
                        "target_kind": target_kind,
                        "enforce_zero_rx": rx,
                        "ms": ms,
                        "err": nan_none(err),
                        "x": x.as_ref().map(|x| array_json(&net.full_positions(x))),
                        "q": finite_vec(&q),
                        "reactions": x.as_ref().map(|x| reactions_json(&q, x)),
                        "reaction_residual": r.diagnostics.reaction_residual,
                    }));
                }
                Err(e) => println!("  reactions {target_kind} {label}: {e}"),
            }
        }
    }
    let value = json!({
        "net": net_json(&net),
        "x_target_jit": array_json(&net.full_positions(&target_for(&net, "jit2%d"))),
        "runs": runs,
    });
    write(&out.join("reactions_tiedarch16s.json"), &value);
}

pub fn cmd_figures(out_dir: &str, max_iters: usize) {
    let out = Path::new(out_dir);
    std::fs::create_dir_all(out).expect("create output directory");
    export_gallery(out);
    export_reactions(out);
    let nets = suite_nets();
    let only = std::env::var("BENCH_NETS").ok();
    for &(name, target_kind, box_kind) in SHOWCASE {
        if let Some(filter) = &only {
            if !filter.split(',').any(|f| f.trim() == name) {
                continue;
            }
        }
        let Some(net) = nets.iter().find(|n| n.name == name) else {
            eprintln!("unknown showcase net {name}");
            continue;
        };
        println!("\n=== {name} {target_kind} {box_kind} ===");
        export_case(net, target_kind, box_kind, max_iters, out);
    }
}
