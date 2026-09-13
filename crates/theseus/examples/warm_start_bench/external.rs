//! Case studies exported from other implementations (JAX FDM, compas_cem).
//!
//! Each JSON file in `bench/external/cases/` (written by the Python scripts
//! in `bench/external/`) describes one structure at a reference equilibrium
//! together with the force densities `q_ref` that produce it, the loads, the
//! supports and — where the external tool was run on the same fit — that
//! tool's own result. The case is turned into a [`Net`] with 3-D loads, the
//! export is checked against Theseus' forward solve, and the usual
//! warm-start → clip → L-BFGS-B table is produced for the exact target
//! (`target`) and, when every free node carries one, the designer's target
//! (`target_original`).

use crate::methods::{clip, lbfgsb, Method, Run};
use crate::nets::{lcg, safe_err, try_forward, Net};
use crate::report::{fmt_e, print_case, Row};
use ndarray::Array2;
use serde::Deserialize;
use std::path::{Path, PathBuf};

/// Per-edge bounds; `null` entries mean "unbounded" on that side.
#[derive(Deserialize)]
struct BoundsJson {
    lo: Vec<Option<f64>>,
    hi: Vec<Option<f64>>,
}

#[derive(Deserialize)]
struct Case {
    name: String,
    #[serde(default)]
    source: String,
    #[serde(default)]
    description: String,
    nodes: Vec<[f64; 3]>,
    edges: Vec<[usize; 2]>,
    fixed: Vec<usize>,
    loads: Vec<[f64; 3]>,
    q_ref: Vec<f64>,
    target: Vec<[f64; 3]>,
    #[serde(default)]
    target_original: Option<Vec<Option<[f64; 3]>>>,
    #[serde(default)]
    bounds: Option<BoundsJson>,
    #[serde(default)]
    external: serde_json::Value,
}

/// Methods run on every external case (the dense Gram path and the guard
/// ablation add nothing here).
pub const EXTERNAL_METHODS: [Method; 7] = [
    Method::Uniform,
    Method::S1,
    Method::GramSparse,
    Method::LengthRatio,
    Method::Frozen,
    Method::Pipeline,
    Method::Legacy,
];

/// Loose box factor when the case carries no bounds (see [`Net::box_from_true`]).
const DEFAULT_BOX_FACTOR: f64 = 100.0;

/// Half-width in decades of the randomised uniform seed.
const UNIFORM_RAND_DECADES: f64 = 0.25;

/// Amplitude of the isotropic target perturbation relative to the bounding-box diagonal.
const JITTER_FRACTION: f64 = 0.01;

fn median_abs(v: &[f64]) -> f64 {
    let mut a: Vec<f64> = v.iter().map(|x| x.abs()).collect();
    a.sort_by(|x, y| x.partial_cmp(y).unwrap());
    if a.is_empty() {
        1.0
    } else {
        a[a.len() / 2]
    }
}

fn bbox_diagonal(nodes: &[[f64; 3]]) -> f64 {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for n in nodes {
        for d in 0..3 {
            lo[d] = lo[d].min(n[d]);
            hi[d] = hi[d].max(n[d]);
        }
    }
    (0..3).map(|d| (hi[d] - lo[d]).powi(2)).sum::<f64>().sqrt()
}

impl Case {
    fn validate(&self) -> Result<(), String> {
        let n = self.nodes.len();
        if self.loads.len() != n || self.target.len() != n {
            return Err(format!("loads/target length != {n} nodes"));
        }
        if self.q_ref.len() != self.edges.len() {
            return Err("q_ref length != edges".into());
        }
        if let Some(b) = &self.bounds {
            if b.lo.len() != self.edges.len() || b.hi.len() != self.edges.len() {
                return Err("bounds length != edges".into());
            }
        }
        if let Some(t) = &self.target_original {
            if t.len() != n {
                return Err("target_original length != nodes".into());
            }
        }
        if self.edges.iter().any(|e| e[0] >= n || e[1] >= n) || self.fixed.iter().any(|&f| f >= n) {
            return Err("edge or support index out of range".into());
        }
        Ok(())
    }

    fn into_net(self) -> (Net, Case) {
        let n = self.nodes.len();
        let mut is_fixed = vec![false; n];
        for &f in &self.fixed {
            is_fixed[f] = true;
        }
        let free: Vec<usize> = (0..n).filter(|&i| !is_fixed[i]).collect();
        let mut fixed = self.fixed.clone();
        fixed.sort_unstable();
        fixed.dedup();
        let load_xyz: Vec<[f64; 3]> = free.iter().map(|&i| self.loads[i]).collect();
        let net = Net {
            name: self.name.clone(),
            edges: self.edges.iter().map(|e| (e[0], e[1])).collect(),
            loads: vec![0.0; free.len()],
            free,
            fixed,
            plan: self.nodes.clone(),
            q_true: self.q_ref.clone(),
            extent: bbox_diagonal(&self.nodes),
            tie_edges: Vec::new(),
            load_xyz: Some(load_xyz),
            seed_magnitude: median_abs(&self.q_ref),
        };
        (net, self)
    }

    /// The case's own box, else the same loose box the synthetic suite uses:
    /// per sign group `[min|q|/100, 100·max|q|]` (so `q_ref` is always inside).
    fn box_(&self, net: &Net) -> (Vec<f64>, Vec<f64>, &'static str) {
        let (loose_lo, loose_hi) = net.box_from_true(DEFAULT_BOX_FACTOR, DEFAULT_BOX_FACTOR);
        if let Some(b) = &self.bounds {
            // Unbounded sides take the loose box's value: L-BFGS-B's direct
            // box parameterisation needs a finite interval.
            let lo: Vec<f64> =
                b.lo.iter()
                    .zip(&loose_lo)
                    .map(|(v, l)| v.unwrap_or(*l))
                    .collect();
            let hi: Vec<f64> =
                b.hi.iter()
                    .zip(&loose_hi)
                    .map(|(v, h)| v.unwrap_or(*h))
                    .collect();
            let name = if b.lo.iter().chain(b.hi.iter()).any(|v| v.is_none()) {
                "case+loose"
            } else {
                "case"
            };
            return (lo, hi, name);
        }
        (loose_lo, loose_hi, "loose")
    }

    fn free_target(&self, net: &Net, full: &[[f64; 3]]) -> Array2<f64> {
        Array2::from_shape_fn((net.free.len(), 3), |(i, d)| full[net.free[i]][d])
    }
}

/// External tool's own position-fit run as a table row (`external.position_fit_run`).
fn external_row(case: &Case) -> Option<Row> {
    let run = case.external.get("position_fit_run")?;
    let f = |k: &str| run.get(k).and_then(|v| v.as_f64());
    let trace: Vec<f64> = run
        .get("loss_trace")
        .and_then(|v| v.as_array())
        .map(|a| a.iter().filter_map(|v| v.as_f64()).collect())
        .unwrap_or_default();
    let final_err = f("final_err")?;
    let ms = f("time_s_excl_jit")
        .or_else(|| f("time_s"))
        .unwrap_or(f64::NAN)
        * 1e3;
    let iters = f("iterations").unwrap_or(0.0) as usize;
    let tool = case
        .external
        .get("tool")
        .and_then(|v| v.as_str())
        .unwrap_or(&case.source)
        .to_string();
    let err_raw = trace.first().map(|l| l.sqrt()).unwrap_or(f64::NAN);
    Some(Row {
        label: format!("ext:{}", tool.split_whitespace().next().unwrap_or("tool")),
        warm_ms: 0.0,
        err_raw,
        err_clip: err_raw,
        n_clipped: 0,
        run: Some(Run {
            trace,
            iters,
            ms,
            final_err,
        }),
        note: run
            .get("optimizer")
            .and_then(|v| v.as_str())
            .map(|s| format!("{s}, wall time of the external tool"))
            .unwrap_or_default(),
    })
}

/// Plain L-BFGS-B from a given seed as a table row.
fn seed_row(
    label: &str,
    net: &Net,
    fixed: &Array2<f64>,
    target: &Array2<f64>,
    q0: &[f64],
    lo: &[f64],
    hi: &[f64],
    max_iters: usize,
    problem: &theseus::types::Problem,
) -> Row {
    let err_raw = safe_err(problem, target, q0);
    let (q, n_clipped, _) = clip(q0, lo, hi);
    let err_clip = safe_err(problem, target, &q);
    let mut row = Row {
        label: label.to_string(),
        warm_ms: 0.0,
        err_raw,
        err_clip,
        n_clipped,
        run: None,
        note: String::new(),
    };
    if !err_clip.is_finite() {
        row.note = "seed gives singular forward solve".into();
        return row;
    }
    match lbfgsb(net, fixed, target, &q, lo, hi, max_iters) {
        Ok(run) => row.run = Some(run),
        Err(e) => row.note = e,
    }
    row
}

fn rows_json(rows: &[Row], l: f64) -> serde_json::Value {
    serde_json::Value::Array(
        rows.iter()
            .map(|r| {
                let (iters, ms, final_err, trace_len) = match &r.run {
                    Some(run) => (
                        Some(run.iters),
                        Some(run.ms),
                        Some(run.final_err / l),
                        run.trace.len(),
                    ),
                    None => (None, None, None, 0),
                };
                serde_json::json!({
                    "label": r.label,
                    "warm_ms": r.warm_ms,
                    "err_raw_over_L": nan_none(r.err_raw / l),
                    "err_clip_over_L": nan_none(r.err_clip / l),
                    "n_clipped": r.n_clipped,
                    "lbfgsb_iters": iters,
                    "lbfgsb_ms": ms,
                    "final_err_over_L": final_err.and_then(nan_none),
                    "trace_len": trace_len,
                    "note": r.note,
                })
            })
            .collect(),
    )
}

fn nan_none(v: f64) -> Option<f64> {
    v.is_finite().then_some(v)
}

fn case_files(args: &[String]) -> Vec<PathBuf> {
    let mut files = Vec::new();
    for a in args {
        let p = Path::new(a);
        if p.is_dir() {
            let mut entries: Vec<PathBuf> = std::fs::read_dir(p)
                .map(|it| {
                    it.filter_map(|e| e.ok().map(|e| e.path()))
                        .filter(|p| p.extension().map(|x| x == "json").unwrap_or(false))
                        .collect()
                })
                .unwrap_or_default();
            entries.sort();
            files.extend(entries);
        } else {
            files.push(p.to_path_buf());
        }
    }
    files
}

fn run_target(
    net: &Net,
    case: &Case,
    label: &str,
    target: &Array2<f64>,
    lo: &[f64],
    hi: &[f64],
    box_name: &str,
    methods: &[Method],
    max_iters: usize,
    with_external: bool,
) -> Vec<Row> {
    let fixed = net.fixed_positions();
    let problem = net.problem(
        &fixed,
        Vec::new(),
        lo,
        hi,
        theseus::types::SolverOptions::default(),
    );
    let mut rows = crate::run_case(net, target, lo, hi, methods, max_iters);
    for (tag, factor) in [("uniform_x0.1", 0.1), ("uniform_x10", 10.0)] {
        let q0: Vec<f64> = net.sign_seed().iter().map(|q| q * factor).collect();
        rows.push(seed_row(
            tag, net, &fixed, target, &q0, lo, hi, max_iters, &problem,
        ));
    }
    // Equal-magnitude mixed-sign seeds put exact zeros on the Laplacian
    // diagonal wherever tension and compression balance at a node; a
    // log-uniform ±¼-decade scatter is the practitioner's way around it.
    let mut state = 0x5eed_u64;
    let q0: Vec<f64> = net
        .sign_seed()
        .iter()
        .map(|q| q * 10f64.powf(UNIFORM_RAND_DECADES * (2.0 * lcg(&mut state) - 1.0)))
        .collect();
    rows.push(seed_row(
        "uniform_rand",
        net,
        &fixed,
        target,
        &q0,
        lo,
        hi,
        max_iters,
        &problem,
    ));
    rows.push(seed_row(
        "oracle(q_ref)",
        net,
        &fixed,
        target,
        &net.q_true,
        lo,
        hi,
        max_iters,
        &problem,
    ));
    if with_external {
        if let Some(r) = external_row(case) {
            rows.push(r);
        }
    }
    print_case(
        net,
        &format!("{label} box={box_name} [{:.3},{:.3}]", lo[0], hi[0]),
        &rows,
        net.extent,
    );
    rows
}

pub fn cmd_external(args: &[String], max_iters: usize, results_dir: Option<&Path>) {
    let files = case_files(args);
    if files.is_empty() {
        eprintln!("external: no case files given");
        return;
    }
    let methods = crate::selected_methods(&EXTERNAL_METHODS);
    for file in files {
        let text = match std::fs::read_to_string(&file) {
            Ok(t) => t,
            Err(e) => {
                eprintln!("{}: {e}", file.display());
                continue;
            }
        };
        let case: Case = match serde_json::from_str(&text) {
            Ok(c) => c,
            Err(e) => {
                eprintln!("{}: parse error: {e}", file.display());
                continue;
            }
        };
        if let Err(e) = case.validate() {
            eprintln!("{}: invalid case: {e}", file.display());
            continue;
        }
        let (net, case) = case.into_net();
        let (lo, hi, box_name) = case.box_(&net);
        let n_t = net.q_true.iter().filter(|&&q| q > 0.0).count();
        let n_c = net.q_true.len() - n_t;
        let outside = net
            .q_true
            .iter()
            .zip(lo.iter().zip(hi.iter()))
            .filter(|(q, (l, h))| **q < **l || **q > **h)
            .count();
        println!(
            "\n##### {} [{}] nodes={} free={} fixed={} edges={} (T={} C={}) L={:.3} median|q|={:.3e} box={} q_ref outside box: {}",
            net.name,
            case.source,
            net.n_nodes(),
            net.free.len(),
            net.fixed.len(),
            net.edges.len(),
            n_t,
            n_c,
            net.extent,
            net.seed_magnitude,
            box_name,
            outside
        );
        if !case.description.is_empty() {
            println!("# {}", case.description);
        }
        // Export consistency: Theseus' forward solve at q_ref must reproduce `nodes`.
        let fixed = net.fixed_positions();
        let problem = net.problem(
            &fixed,
            Vec::new(),
            &lo,
            &hi,
            theseus::types::SolverOptions::default(),
        );
        let reference = case.free_target(&net, &case.nodes);
        match try_forward(&problem, &net.q_true) {
            Some(x) => {
                let dev = (&x - &reference).iter().map(|v| v * v).sum::<f64>().sqrt();
                println!(
                    "# forward(q_ref) vs exported nodes: ||dx||/L = {}",
                    fmt_e(dev / net.extent)
                );
                if dev / net.extent > 1e-6 {
                    println!("# WARNING: export is not an FDM equilibrium under Theseus' conventions; results below fit the exported geometry anyway");
                }
            }
            None => println!("# WARNING: forward solve at q_ref failed (singular)"),
        }
        let mut out = serde_json::Map::new();
        out.insert("name".into(), serde_json::json!(net.name));
        out.insert("source".into(), serde_json::json!(case.source));
        out.insert("extent".into(), serde_json::json!(net.extent));
        out.insert("edges".into(), serde_json::json!(net.edges.len()));
        out.insert("free".into(), serde_json::json!(net.free.len()));
        out.insert("fixed".into(), serde_json::json!(net.fixed.len()));
        out.insert("tension".into(), serde_json::json!(n_t));
        out.insert("compression".into(), serde_json::json!(n_c));
        out.insert("box".into(), serde_json::json!(box_name));
        out.insert("q_ref_outside_box".into(), serde_json::json!(outside));

        let target = case.free_target(&net, &case.target);
        let rows = run_target(
            &net, &case, "exact", &target, &lo, &hi, box_name, &methods, max_iters, true,
        );
        out.insert("exact".into(), rows_json(&rows, net.extent));

        // Perturbed target: isotropic white noise on every free node, so that
        // no q reproduces it exactly (the realistic situation; Stage 1 alone
        // cannot solve it). Amplitude JITTER_FRACTION·L, capped at a fifth of
        // the shortest edge so that short auxiliary members are not folded.
        let min_len = net
            .edge_lengths(&reference)
            .into_iter()
            .filter(|l| l.is_finite() && *l > 0.0)
            .fold(f64::INFINITY, f64::min);
        let amplitude = (JITTER_FRACTION * net.extent).min(0.2 * min_len);
        let mut state = 0xc0ffee_u64;
        let mut jittered = target.clone();
        for v in jittered.iter_mut() {
            *v += amplitude * (2.0 * lcg(&mut state) - 1.0);
        }
        let label = format!("jit{:.2}%L", 100.0 * amplitude / net.extent);
        let rows = run_target(
            &net, &case, &label, &jittered, &lo, &hi, box_name, &methods, max_iters, false,
        );
        out.insert("jitter".into(), rows_json(&rows, net.extent));
        out.insert(
            "jitter_amplitude_over_L".into(),
            serde_json::json!(amplitude / net.extent),
        );

        if let Some(t) = &case.target_original {
            let complete = net.free.iter().all(|&i| t[i].is_some());
            if complete {
                let full: Vec<[f64; 3]> = (0..net.n_nodes())
                    .map(|i| t[i].unwrap_or(case.nodes[i]))
                    .collect();
                let target = case.free_target(&net, &full);
                let rows = run_target(
                    &net, &case, "designer", &target, &lo, &hi, box_name, &methods, max_iters,
                    false,
                );
                out.insert("designer".into(), rows_json(&rows, net.extent));
            } else {
                let n_some = net.free.iter().filter(|&&i| t[i].is_some()).count();
                println!(
                    "# target_original covers {n_some}/{} free nodes — partial targets are not run",
                    net.free.len()
                );
            }
        }
        if let Some(dir) = results_dir {
            let _ = std::fs::create_dir_all(dir);
            let path = dir.join(format!("rust_{}.json", net.name));
            match serde_json::to_string_pretty(&serde_json::Value::Object(out)) {
                Ok(s) => {
                    if let Err(e) = std::fs::write(&path, s) {
                        eprintln!("cannot write {}: {e}", path.display());
                    }
                }
                Err(e) => eprintln!("cannot serialise results: {e}"),
            }
        }
    }
}
