//! Warm-start methods under test, plus the downstream L-BFGS-B runner.
//!
//! Every method maps `(net, target, box)` to a force-density vector; the
//! harness then clips it to the box, scores it and hands it to L-BFGS-B.
//! Library calls are wrapped in `catch_unwind` so that a failure in one method
//! is reported in its row and never aborts the run.

use crate::nets::{safe_err, try_forward, Net};
use ndarray::Array2;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::sync::atomic::AtomicBool;
use std::time::Instant;
use theseus::inverse::{
    solve_inverse_fdm, InverseFdmOptions, InverseFdmResult, InverseMetric, LinearAlgebra,
    ParticularMethod, Stage2Method, DEFAULT_LM_DAMPING, DEFAULT_SEED_GUARD_MARGIN,
};
use theseus::nullspace::{EquilibriumSystem, EquilibriumUnknown};
use theseus::sparse::SparseColMatOwned;
use theseus::types::{
    FdmCache, ObjectiveTrait, OptimizationState, Problem, QParameterizationMode, SolverOptions,
    TargetGeometryReduction, TargetXYZ,
};

/// Dense Gram assembly is skipped above this many edges (ne² doubles).
pub const DENSE_EDGE_CAP: usize = 6000;
/// Relative Tikhonov shift for the unboxed Gram solves, scaled by the mean
/// squared column norm of the equilibrium matrix.
pub const GRAM_RELATIVE_SHIFT: f64 = 1e-8;
/// Length-ratio heuristic iteration budget.
pub const LENGTH_RATIO_ITERS: usize = 20;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Method {
    /// Uniform magnitude with the true sign pattern; no inverse solve.
    Uniform,
    /// Stage 1 only: boxed Clarabel force-residual particular (member-force form).
    S1,
    /// Unboxed sparse Gram normal equations `(EᵀE + λI) q = Eᵀp`, then clip.
    GramSparse,
    /// The same normal equations formed densely with an O(ne³) Cholesky.
    GramDense,
    /// Practitioner heuristic `q ← q · ℓ(q)/ℓ*`, one forward solve per sweep.
    LengthRatio,
    /// One compliance-weighted frozen-Jacobian step (`InverseMetric::Geometry`).
    Frozen,
    /// New default pipeline: Stage 1 → guard → 1 frozen + 2 Gauss–Newton steps.
    Pipeline,
    /// Previous branch behaviour: Clarabel Stage 2, no LM damping, no guard, dimensional.
    Legacy,
    /// The pipeline with the Stage-1 collapse guard disabled.
    PipelineNoguard,
}

impl Method {
    pub const ALL: [Method; 9] = [
        Method::Uniform,
        Method::S1,
        Method::GramSparse,
        Method::GramDense,
        Method::LengthRatio,
        Method::Frozen,
        Method::Pipeline,
        Method::Legacy,
        Method::PipelineNoguard,
    ];

    pub fn name(self) -> &'static str {
        match self {
            Method::Uniform => "uniform",
            Method::S1 => "s1",
            Method::GramSparse => "gram_sparse",
            Method::GramDense => "gram_dense",
            Method::LengthRatio => "length_ratio",
            Method::Frozen => "frozen",
            Method::Pipeline => "pipeline",
            Method::Legacy => "legacy",
            Method::PipelineNoguard => "pipeline_noguard",
        }
    }

    pub fn from_name(name: &str) -> Option<Method> {
        Method::ALL.iter().copied().find(|m| m.name() == name)
    }
}

pub struct WarmOut {
    pub q: Vec<f64>,
    pub ms: f64,
    pub note: String,
}

/// Common `solve_inverse_fdm` options for the boxed library methods.
pub fn base_options(lo: &[f64], hi: &[f64]) -> InverseFdmOptions {
    InverseFdmOptions {
        regularization: 0.0,
        use_l2: true,
        max_l1_iter: 1,
        particular_method: ParticularMethod::Clarabel,
        linear_algebra: LinearAlgebra::Direct,
        enforce_zero_rx: false,
        enforce_zero_ry: false,
        enforce_zero_rz: false,
        solve_for_q: false,
        signs: Vec::new(),
        lower: lo.to_vec(),
        upper: hi.to_vec(),
        max_iter: 4000,
        tol: 1e-8,
        metric: InverseMetric::Force,
        q_ref: Vec::new(),
        max_frozen_outer: 0,
        max_outer: 0,
        cwls_damping: 1e-6,
        stage2_method: Stage2Method::ActiveSet,
        lm_damping: DEFAULT_LM_DAMPING,
        seed_guard_margin: DEFAULT_SEED_GUARD_MARGIN,
        nondimensionalize: true,
        reaction_weight: 1.0,
    }
}

fn env_f64(var: &str) -> Option<f64> {
    std::env::var(var).ok().and_then(|s| s.trim().parse().ok())
}

/// Ablation overrides for `pipeline` (only), read from the environment:
/// `BENCH_LM` (lm_damping), `BENCH_GUARD` (seed_guard_margin, 0 disables),
/// `BENCH_S2=clarabel|activeset`, `BENCH_NONDIM=0|1`, `BENCH_CWLS` (cwls_damping).
fn apply_pipeline_overrides(mut opts: InverseFdmOptions) -> InverseFdmOptions {
    if let Some(v) = env_f64("BENCH_LM") {
        opts.lm_damping = v;
    }
    if let Some(v) = env_f64("BENCH_GUARD") {
        opts.seed_guard_margin = v;
    }
    if let Some(v) = env_f64("BENCH_CWLS") {
        opts.cwls_damping = v;
    }
    if let Ok(s) = std::env::var("BENCH_S2") {
        opts.stage2_method = match s.trim() {
            "clarabel" => Stage2Method::Clarabel,
            _ => Stage2Method::ActiveSet,
        };
    }
    if let Ok(s) = std::env::var("BENCH_NONDIM") {
        opts.nondimensionalize = s.trim() != "0";
    }
    opts
}

/// Library options for a method, or `None` for methods that do not call the library.
pub fn library_options(
    method: Method,
    lo: &[f64],
    hi: &[f64],
    gram_shift: f64,
) -> Option<InverseFdmOptions> {
    let base = base_options(lo, hi);
    Some(match method {
        Method::S1 => base,
        Method::GramSparse => InverseFdmOptions {
            regularization: gram_shift,
            particular_method: ParticularMethod::Gram,
            solve_for_q: true,
            lower: Vec::new(),
            upper: Vec::new(),
            ..base
        },
        Method::Frozen => InverseFdmOptions {
            metric: InverseMetric::Geometry,
            max_frozen_outer: 0,
            max_outer: 1,
            ..base
        },
        Method::Pipeline => apply_pipeline_overrides(InverseFdmOptions {
            metric: InverseMetric::GeometryNewton,
            max_frozen_outer: 1,
            max_outer: 2,
            ..base
        }),
        Method::Legacy => InverseFdmOptions {
            metric: InverseMetric::GeometryNewton,
            max_frozen_outer: 1,
            max_outer: 2,
            stage2_method: Stage2Method::Clarabel,
            lm_damping: 0.0,
            seed_guard_margin: 0.0,
            nondimensionalize: false,
            ..base
        },
        Method::PipelineNoguard => InverseFdmOptions {
            metric: InverseMetric::GeometryNewton,
            max_frozen_outer: 1,
            max_outer: 2,
            seed_guard_margin: 0.0,
            ..base
        },
        Method::Uniform | Method::GramDense | Method::LengthRatio => return None,
    })
}

/// Shorten an error message to a single table-friendly line.
pub fn short_error(s: &str) -> String {
    let first = s.split(';').next().unwrap_or("").trim().to_string();
    if first.chars().count() > 90 {
        let cut: String = first.chars().take(87).collect();
        format!("{cut}...")
    } else {
        first
    }
}

/// Run `solve_inverse_fdm`, catching both `Err` and panics. The error text is
/// returned in full; use [`short_error`] for table notes.
pub fn run_library(
    problem: &Problem,
    target: &Array2<f64>,
    opts: InverseFdmOptions,
) -> Result<(InverseFdmResult, f64), String> {
    let started = Instant::now();
    let out = catch_unwind(AssertUnwindSafe(|| {
        solve_inverse_fdm(problem, target, opts)
    }));
    let ms = started.elapsed().as_secs_f64() * 1e3;
    match out {
        Ok(Ok(r)) => Ok((r, ms)),
        Ok(Err(e)) => Err(format!("ERR {e}")),
        Err(p) => Err(format!("PANIC {}", panic_message(&p))),
    }
}

pub fn panic_message(p: &Box<dyn std::any::Any + Send>) -> String {
    if let Some(s) = p.downcast_ref::<&str>() {
        s.to_string()
    } else if let Some(s) = p.downcast_ref::<String>() {
        s.clone()
    } else {
        "unknown panic".into()
    }
}

/// Compact rendering of `InverseDiagnostics` for the note column. Stage-2
/// fields are only shown for geometric solves (`stage2 = true`).
pub fn diag_note(r: &InverseFdmResult, stage2: bool) -> String {
    let d = &r.diagnostics;
    let mut s = format!(
        "it={}{} lib_err={:.2e}",
        r.iterations,
        if r.converged { "" } else { "!" },
        r.geometric_error
    );
    if !stage2 {
        return s;
    }
    s.push_str(&format!(" s1e={:.2e}", d.stage1_error));
    if d.uniform_seed_error.is_finite() {
        s.push_str(&format!(" u={:.2e}", d.uniform_seed_error));
    }
    if d.used_uniform_seed {
        s.push_str(" guard→u");
    }
    s.push_str(&format!(
        " fr{} gn{} fac={} fb={}",
        d.frozen_steps, d.newton_steps, d.stage2_factorizations, d.clarabel_fallbacks
    ));
    if d.active_set_capped > 0 {
        s.push_str(&format!(" cap={}", d.active_set_capped));
    }
    if d.degenerate_linearizations > 0 {
        s.push_str(&format!(" degen={}", d.degenerate_linearizations));
    }
    if d.reaction_residual != 0.0 {
        s.push_str(&format!(" R*={:.2e}", d.reaction_residual));
    }
    s
}

/// Force-density-form equilibrium matrix `E` at the target and the mean
/// squared column norm used to scale the Gram shift.
pub fn equilibrium_e(
    problem: &Problem,
    target: &Array2<f64>,
) -> Result<(SparseColMatOwned, Vec<f64>, f64), String> {
    let sys = EquilibriumSystem::assemble(
        problem,
        target,
        EquilibriumUnknown::ForceDensity,
        false,
        false,
        false,
    )
    .map_err(|e| format!("ERR {e}"))?;
    let ne = sys.a.ncols.max(1) as f64;
    let scale = sys.a.values.iter().map(|v| v * v).sum::<f64>() / ne;
    Ok((sys.a, sys.p, scale))
}

/// Dispatch a warm-start method. `lo`/`hi` is the L-BFGS-B box (also the
/// inverse-solve box for the boxed methods). Errors are shortened to one line.
pub fn run_method(
    method: Method,
    net: &Net,
    problem: &Problem,
    target: &Array2<f64>,
    lo: &[f64],
    hi: &[f64],
) -> Result<WarmOut, String> {
    run_method_inner(method, net, problem, target, lo, hi).map_err(|e| short_error(&e))
}

fn run_method_inner(
    method: Method,
    net: &Net,
    problem: &Problem,
    target: &Array2<f64>,
    lo: &[f64],
    hi: &[f64],
) -> Result<WarmOut, String> {
    match method {
        Method::Uniform => Ok(WarmOut {
            q: net.sign_seed(),
            ms: 0.0,
            note: String::new(),
        }),
        Method::GramDense => gram_dense(problem, target),
        Method::LengthRatio => length_ratio(net, problem, target, lo, hi),
        Method::GramSparse => {
            let started = Instant::now();
            let (_, _, scale) = equilibrium_e(problem, target)?;
            let shift = GRAM_RELATIVE_SHIFT * scale;
            let opts = library_options(method, lo, hi, shift).unwrap();
            let (r, _) = run_library(problem, target, opts)?;
            Ok(WarmOut {
                note: format!("λ={shift:.1e} {}", diag_note(&r, false)),
                q: r.q,
                ms: started.elapsed().as_secs_f64() * 1e3,
            })
        }
        _ => {
            let opts = library_options(method, lo, hi, 0.0).unwrap();
            let stage2 = opts.metric.is_geometric();
            let (r, ms) = run_library(problem, target, opts)?;
            Ok(WarmOut {
                note: diag_note(&r, stage2),
                q: r.q,
                ms,
            })
        }
    }
}

// ───────────────────────── dense Gram ─────────────────────────

pub struct DenseGram {
    pub q: Vec<f64>,
    pub bytes: usize,
    pub assemble_ms: f64,
    pub factor_ms: f64,
}

/// Form `EᵀE` as a dense row-major `ne × ne` array from the CSC `E`, add a
/// Tikhonov shift, and solve with an in-place Cholesky. Returns `Err` above
/// [`DENSE_EDGE_CAP`] edges.
pub fn dense_gram_solve(e: &SparseColMatOwned, p: &[f64], shift: f64) -> Result<DenseGram, String> {
    let ne = e.ncols;
    if ne > DENSE_EDGE_CAP {
        return Err(format!("skipped (ne > cap {DENSE_EDGE_CAP})"));
    }
    let started = Instant::now();
    // row lists of E for the pairwise products
    let mut rows: Vec<Vec<(usize, f64)>> = vec![Vec::new(); e.nrows];
    for col in 0..ne {
        for nz in e.col_ptrs[col] as usize..e.col_ptrs[col + 1] as usize {
            rows[e.row_indices[nz] as usize].push((col, e.values[nz]));
        }
    }
    let mut g = vec![0.0f64; ne * ne];
    let mut rhs = vec![0.0f64; ne];
    for (i, row) in rows.iter().enumerate() {
        for &(a, va) in row {
            rhs[a] += va * p[i];
            for &(b, vb) in row {
                g[a * ne + b] += va * vb;
            }
        }
    }
    for i in 0..ne {
        g[i * ne + i] += shift;
    }
    let assemble_ms = started.elapsed().as_secs_f64() * 1e3;
    let started = Instant::now();
    dense_cholesky_in_place(&mut g, ne)?;
    dense_cholesky_solve(&g, ne, &mut rhs);
    let factor_ms = started.elapsed().as_secs_f64() * 1e3;
    Ok(DenseGram {
        q: rhs,
        bytes: ne * ne * std::mem::size_of::<f64>(),
        assemble_ms,
        factor_ms,
    })
}

/// Row-major Cholesky–Banachiewicz; the lower triangle is overwritten by `L`.
fn dense_cholesky_in_place(g: &mut [f64], n: usize) -> Result<(), String> {
    for i in 0..n {
        for j in 0..=i {
            let (head, tail) = g.split_at_mut(i * n);
            let row_i = &mut tail[..n];
            let row_j: &[f64] = if j == i {
                &row_i[..j]
            } else {
                &head[j * n..j * n + j]
            };
            let dot: f64 = row_i[..j].iter().zip(row_j).map(|(a, b)| a * b).sum();
            let s = row_i[j] - dot;
            if j == i {
                if s <= 0.0 || !s.is_finite() {
                    return Err(format!("dense Cholesky failed at pivot {i} (s={s:.3e})"));
                }
                row_i[i] = s.sqrt();
            } else {
                row_i[j] = s / head[j * n + j];
            }
        }
    }
    Ok(())
}

/// Solve `L Lᵀ x = b` in place given the row-major lower factor.
fn dense_cholesky_solve(l: &[f64], n: usize, b: &mut [f64]) {
    for i in 0..n {
        let row = &l[i * n..i * n + i];
        let dot: f64 = row.iter().zip(&b[..i]).map(|(a, x)| a * x).sum();
        b[i] = (b[i] - dot) / l[i * n + i];
    }
    for i in (0..n).rev() {
        let mut s = b[i];
        for k in i + 1..n {
            s -= l[k * n + i] * b[k];
        }
        b[i] = s / l[i * n + i];
    }
}

fn gram_dense(problem: &Problem, target: &Array2<f64>) -> Result<WarmOut, String> {
    let started = Instant::now();
    let (e, p, scale) = equilibrium_e(problem, target)?;
    let shift = GRAM_RELATIVE_SHIFT * scale;
    let out = catch_unwind(AssertUnwindSafe(|| dense_gram_solve(&e, &p, shift)));
    let dg = match out {
        Ok(r) => r?,
        Err(pnc) => return Err(format!("PANIC {}", panic_message(&pnc))),
    };
    Ok(WarmOut {
        q: dg.q,
        ms: started.elapsed().as_secs_f64() * 1e3,
        note: format!(
            "λ={shift:.1e} bytes={} asm={:.0}ms chol={:.0}ms",
            human_bytes(dg.bytes),
            dg.assemble_ms,
            dg.factor_ms
        ),
    })
}

pub fn human_bytes(b: usize) -> String {
    let b = b as f64;
    if b >= 1e9 {
        format!("{:.2}GB", b / 1e9)
    } else if b >= 1e6 {
        format!("{:.1}MB", b / 1e6)
    } else if b >= 1e3 {
        format!("{:.0}kB", b / 1e3)
    } else {
        format!("{b:.0}B")
    }
}

// ───────────────────────── length-ratio heuristic ─────────────────────────

/// `q_i ← q_i · ℓ_i(q)/ℓ_i*` from the clipped uniform seed, clipping to the
/// box after every sweep and stopping when the error stops improving.
fn length_ratio(
    net: &Net,
    problem: &Problem,
    target: &Array2<f64>,
    lo: &[f64],
    hi: &[f64],
) -> Result<WarmOut, String> {
    let started = Instant::now();
    let target_len = net.edge_lengths(target);
    let (mut q, _, _) = clip(&net.sign_seed(), lo, hi);
    let mut best = q.clone();
    let mut best_err = f64::INFINITY;
    let mut best_k = 0;
    let mut sweeps = 0;
    for k in 0..=LENGTH_RATIO_ITERS {
        let x = match try_forward(problem, &q) {
            Some(x) => x,
            None => break,
        };
        let err = (&x - target).iter().map(|v| v * v).sum::<f64>().sqrt();
        if err < best_err {
            best_err = err;
            best = q.clone();
            best_k = k;
        } else {
            break;
        }
        if k == LENGTH_RATIO_ITERS {
            break;
        }
        let len = net.edge_lengths(&x);
        for i in 0..q.len() {
            if target_len[i] > 0.0 && len[i].is_finite() {
                q[i] *= len[i] / target_len[i];
            }
        }
        q = clip(&q, lo, hi).0;
        sweeps += 1;
    }
    Ok(WarmOut {
        q: best,
        ms: started.elapsed().as_secs_f64() * 1e3,
        note: format!("sweeps={sweeps} best@{best_k}"),
    })
}

// ───────────────────────── clipping ─────────────────────────

/// Clamp `q` into `[lo, hi]`; returns `(clipped q, #changed, #on a bound)`.
pub fn clip(q: &[f64], lo: &[f64], hi: &[f64]) -> (Vec<f64>, usize, usize) {
    let mut clipped = 0;
    let mut active = 0;
    let out: Vec<f64> = q
        .iter()
        .zip(lo)
        .zip(hi)
        .map(|((&v, &a), &b)| {
            let c = if v.is_finite() {
                v.clamp(a, b)
            } else {
                0.5 * (a + b)
            };
            if c != v {
                clipped += 1;
            }
            if (c - a).abs() <= 1e-9 * a.abs().max(1e-12)
                || (c - b).abs() <= 1e-9 * b.abs().max(1e-12)
            {
                active += 1;
            }
            c
        })
        .collect();
    (out, clipped, active)
}

// ───────────────────────── downstream L-BFGS-B ─────────────────────────

pub struct Run {
    pub trace: Vec<f64>,
    pub iters: usize,
    pub ms: f64,
    pub final_err: f64,
}

fn target_objective(net: &Net, target: &Array2<f64>) -> Box<dyn ObjectiveTrait> {
    Box::new(TargetXYZ {
        weight: 1.0,
        node_indices: net.free.clone(),
        target: target.clone(),
        reduction: TargetGeometryReduction::Sse,
    })
}

/// Box-constrained L-BFGS-B on the SSE target objective from `q0`.
pub fn lbfgsb(
    net: &Net,
    fixed: &Array2<f64>,
    target: &Array2<f64>,
    q0: &[f64],
    lo: &[f64],
    hi: &[f64],
    max_iters: usize,
) -> Result<Run, String> {
    let problem = net.problem(
        fixed,
        vec![target_objective(net, target)],
        lo,
        hi,
        SolverOptions {
            absolute_tolerance: 1e-8,
            relative_tolerance: 1e-12,
            max_iterations: max_iters,
            q_parameterization_mode: QParameterizationMode::DirectBoxBounds,
            ..SolverOptions::default()
        },
    );
    let mut state = OptimizationState::new(q0.to_vec(), Array2::zeros((0, 3)));
    let cancel = AtomicBool::new(false);
    let started = Instant::now();
    let out = catch_unwind(AssertUnwindSafe(|| {
        theseus::optimizer::optimize(&problem, &mut state, None, 0, &cancel)
    }));
    let ms = started.elapsed().as_secs_f64() * 1e3;
    let result = match out {
        Ok(Ok(r)) => r,
        Ok(Err(e)) => return Err(short_error(&format!("lbfgsb ERR {e}"))),
        Err(p) => return Err(short_error(&format!("lbfgsb PANIC {}", panic_message(&p)))),
    };
    let final_err = safe_err(&problem, target, &result.q);
    Ok(Run {
        trace: result.loss_trace,
        iters: result.iterations,
        ms,
        final_err,
    })
}

/// Wall time of one forward solve (cache construction excluded).
pub fn time_forward(problem: &Problem, q: &[f64]) -> f64 {
    let mut cache = FdmCache::new(problem).expect("cache");
    let started = Instant::now();
    theseus::fdm::solve_fdm(&mut cache, q, problem, &Array2::zeros((0, 3)), 0.0).expect("forward");
    started.elapsed().as_secs_f64() * 1e3
}

/// Wall time of one L-BFGS-B evaluation (loss + adjoint gradient).
pub fn time_eval(
    net: &Net,
    fixed: &Array2<f64>,
    target: &Array2<f64>,
    q: &[f64],
    lo: &[f64],
    hi: &[f64],
) -> f64 {
    let problem = net.problem(
        fixed,
        vec![target_objective(net, target)],
        lo,
        hi,
        SolverOptions {
            q_parameterization_mode: QParameterizationMode::DirectBoxBounds,
            ..SolverOptions::default()
        },
    );
    let mut cache = FdmCache::new(&problem).expect("cache");
    let mut grad = vec![0.0; q.len()];
    let started = Instant::now();
    theseus::gradients::value_and_gradient(&mut cache, &problem, q, &mut grad, lo, hi, &[], &[])
        .expect("eval");
    started.elapsed().as_secs_f64() * 1e3
}
