//! Table output. Every case prints one header line starting with `config/box`
//! and one row per method, so the output can be parsed line by line.

use crate::methods::Run;
use crate::nets::Net;

pub struct Row {
    pub label: String,
    pub warm_ms: f64,
    pub err_raw: f64,
    pub err_clip: f64,
    pub n_clipped: usize,
    /// `None` when the warm start or L-BFGS-B failed; `note` says why.
    pub run: Option<Run>,
    pub note: String,
}

impl Row {
    pub fn failed(label: &str, warm_ms: f64, note: String) -> Self {
        Row {
            label: label.to_string(),
            warm_ms,
            err_raw: f64::NAN,
            err_clip: f64::NAN,
            n_clipped: 0,
            run: None,
            note,
        }
    }
}

/// `{:.3e}` or `-` for NaN / missing values.
pub fn fmt_e(v: f64) -> String {
    if v.is_finite() {
        format!("{v:.3e}")
    } else {
        "-".into()
    }
}

pub const HEADER: &str = "config/box";

/// Evaluations until the loss trace first drops within `factor` of `reference`.
fn evals_to(run: &Run, reference: f64, factor: f64) -> Option<usize> {
    let loss_t = (factor * reference).powi(2);
    run.trace.iter().position(|&v| v <= loss_t)
}

pub fn print_case(net: &Net, case: &str, rows: &[Row], l: f64) {
    // Reference rows (the oracle started at q_ref, an external tool's run)
    // are excluded from the "best" the eval-count columns are measured against.
    let best = rows
        .iter()
        .filter(|r| !r.label.starts_with("oracle") && !r.label.starts_with("ext:"))
        .filter_map(|r| r.run.as_ref().map(|run| run.final_err))
        .filter(|v| v.is_finite())
        .fold(f64::INFINITY, f64::min);
    let thresholds = [1.5, 1.05, 1.005];
    println!(
        "\n=== {} {} (ne={}, nfree={})  best final err/L = {} ===",
        net.name,
        case,
        net.edges.len(),
        net.free.len(),
        fmt_e(best / l)
    );
    println!(
        "{:<18}{:>9}{:>11}{:>11}{:>6}{:>7}{:>9}{:>11}{:>9}{:>10}{:>11}{:>10}  {}",
        HEADER,
        "warm_ms",
        "err_raw/L",
        "err_clip/L",
        "clip",
        "lb_it",
        "lb_ms",
        "final/L",
        "ev<1.5b",
        "ev<1.05b",
        "ev<1.005b",
        "ms->1.05",
        "note"
    );
    for r in rows {
        let (lb_it, lb_ms, final_err, ev, ms_to) = match &r.run {
            Some(run) => {
                let ev: Vec<String> = thresholds
                    .iter()
                    .map(|&t| {
                        evals_to(run, best, t)
                            .map(|p| p.to_string())
                            .unwrap_or_else(|| "-".into())
                    })
                    .collect();
                let ms_to = evals_to(run, best, 1.05)
                    .map(|p| {
                        format!(
                            "{:.1}",
                            r.warm_ms + run.ms * p as f64 / run.trace.len().max(1) as f64
                        )
                    })
                    .unwrap_or_else(|| "-".into());
                (
                    run.iters.to_string(),
                    format!("{:.0}", run.ms),
                    run.final_err,
                    ev,
                    ms_to,
                )
            }
            None => (
                "-".into(),
                "-".into(),
                f64::NAN,
                vec!["-".into(); 3],
                "-".into(),
            ),
        };
        println!(
            "{:<18}{:>9.1}{:>11}{:>11}{:>6}{:>7}{:>9}{:>11}{:>9}{:>10}{:>11}{:>10}  {}",
            r.label,
            r.warm_ms,
            fmt_e(r.err_raw / l),
            fmt_e(r.err_clip / l),
            r.n_clipped,
            lb_it,
            lb_ms,
            fmt_e(final_err / l),
            ev[0],
            ev[1],
            ev[2],
            ms_to,
            r.note
        );
    }
}
