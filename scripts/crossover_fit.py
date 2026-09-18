#!/usr/bin/env python3
"""Fit the §5.3 Stage-B cost models to a `runs.jsonl` produced by
`bench_sweep.py` / `tests/bench_scale.rs`, and — when both a direct and an
iterative backend are present — solve for the crossover edge count with
bootstrap 90% intervals.

Usage (from the repository root):

    uv run --project scripts scripts/crossover_fit.py \\
        benchmarks/results/<machine-id>/<YYYYMMDD>-<sha7>/runs.jsonl \\
        [--out benchmarks/reports] [--bootstrap 500] [--seed 0] [--no-plot]

Outputs `benchmarks/reports/crossover-<machine-id>.md`, `.json` and a
log-log plot of the data with the fitted curves (`crossover-<machine-id>.svg`).

Models (per fixture class and thread setting; `n` = edges, `t` in ms):

    Direct:     t = a·n^1.5 + b·n·log(n) + c
    Iterative:  t = (d + e·iters)·n + f, with iters = i0 + i1·log(n) per fixture

Both are fitted by least squares on log-scaled data (residuals
`log t_model − log t_obs`, i.e. relative error), parameters constrained to be
non-negative. The direct model is fitted to the per-evaluation time (`cold`
and, for the direct backend, identically `warm_2pct` since a refactor does
not warm-start) and to the full-solve time (`solve<iters>`, the recorded
iteration budget, extrapolated to 40 iterations as `solve40`). Residuals are
reported per point and as RMS / max relative error; the program expects
< 15%.

Crossover = smallest `n` at which the fitted iterative time is below the
fitted direct time. The 90% interval comes from refitting on bootstrap
resamples of the per-run evaluation samples (or of the runs within a cell
for the full-solve metric).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares, nnls

sys.path.insert(0, str(Path(__file__).resolve().parent))
from bench_common import (  # noqa: E402
    eval_median,
    eval_samples,
    fixture_sort_key,
    fmt_ms,
    group_by,
    load_runs,
)

DIRECT_BACKENDS = {"direct", "0"}
SOLVE_TARGET_ITERS = 40


# ── models ─────────────────────────────────────────────────────────────────


def direct_model(p, n):
    a, b, c = p
    return a * n**1.5 + b * n * np.log(n) + c


def iterative_model(p, n, iters):
    d, e, f = p
    return (d + e * iters) * n + f


def _fit_log(model, n, t, p0, extra=()):
    n = np.asarray(n, dtype=float)
    t = np.asarray(t, dtype=float)

    def resid(p):
        pred = model(p, n, *extra)
        return np.log(np.maximum(pred, 1e-12)) - np.log(t)

    res = least_squares(resid, p0, bounds=(0.0, np.inf), x_scale="jac", max_nfev=5000)
    pred = model(res.x, n, *extra)
    rel = pred / t - 1.0
    return {
        "params": [float(v) for v in res.x],
        "rel_residuals": [float(v) for v in rel],
        "rms_rel": float(np.sqrt(np.mean(rel**2))),
        "max_rel": float(np.max(np.abs(rel))),
        "n_points": int(len(n)),
        "converged": bool(res.success),
    }


def _nnls_start(columns, t):
    """Both models are linear in their parameters, so a non-negative linear
    least-squares solve (weighted by 1/t to approximate the log-scale
    objective) gives a robust starting point for the log-scale refinement."""
    A = np.vstack(columns).T / t[:, None]
    p0, _ = nnls(A, np.ones_like(t))
    return [max(float(v), 1e-15) for v in p0]


def fit_direct(n, t):
    n = np.asarray(n, dtype=float)
    t = np.asarray(t, dtype=float)
    p0 = _nnls_start([n**1.5, n * np.log(n), np.ones_like(n)], t)
    return _fit_log(direct_model, n, t, p0)


def fit_iterative(n, t, iters):
    n = np.asarray(n, dtype=float)
    t = np.asarray(t, dtype=float)
    iters = np.asarray(iters, dtype=float)
    p0 = _nnls_start([n, iters * n, np.ones_like(n)], t)
    return _fit_log(iterative_model, n, t, p0, extra=(iters,))


def fit_iters(n, iters):
    """`iters = i0 + i1·log n`, plain least squares (should be ~flat)."""
    n = np.asarray(n, dtype=float)
    A = np.vstack([np.ones_like(n), np.log(n)]).T
    coef, *_ = np.linalg.lstsq(A, np.asarray(iters, dtype=float), rcond=None)
    return [float(coef[0]), float(coef[1])]


def crossover(direct_p, iter_p, iters_coef, n_lo=1e3, n_hi=1e9):
    """Smallest n in [n_lo, n_hi] with iterative(n) < direct(n), or None."""
    grid = np.logspace(math.log10(n_lo), math.log10(n_hi), 400)
    it = iterative_model(iter_p, grid, iters_coef[0] + iters_coef[1] * np.log(grid))
    dr = direct_model(direct_p, grid)
    below = it < dr
    if not below.any():
        return None
    k = int(np.argmax(below))
    if k == 0:
        return float(grid[0])
    lo, hi = grid[k - 1], grid[k]
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        if iterative_model(iter_p, mid, iters_coef[0] + iters_coef[1] * math.log(mid)) < direct_model(direct_p, mid):
            hi = mid
        else:
            lo = mid
    return float(hi)


# ── data extraction ────────────────────────────────────────────────────────


def _is_direct(backend: str) -> bool:
    return str(backend).lower() in DIRECT_BACKENDS


def solve40_ms(rec: dict) -> float:
    """Full-solve time extrapolated to 40 iterations: setup + per-iteration cost × 40."""
    iters = max(int(rec.get("iterations") or 1), 1)
    setup = float(rec.get("setup_ms") or 0.0)
    return setup + (float(rec["total_ms"]) - setup) * SOLVE_TARGET_ITERS / iters


METRICS = {
    "cold": ("per-evaluation time, ms (median of the timed evaluations)", eval_median),
    "warm_2pct": ("per-evaluation time after a 2% q step, ms (direct backend: same as cold)", eval_median),
    "solve40": (f"full solve extrapolated to {SOLVE_TARGET_ITERS} L-BFGS-B iterations, ms", solve40_ms),
}


def cell_series(runs: list[dict], metric: str):
    """Median per (edges) cell of a metric, plus the per-run values for bootstrap."""
    fn = METRICS[metric][1]
    cells = group_by(runs, "edges")
    n, t, iters, per_cell = [], [], [], []
    for (edges,), rs in sorted(cells.items()):
        vals = [fn(r) for r in rs]
        n.append(float(edges))
        t.append(float(np.median(vals)))
        its = [r.get("linear_solver_iterations") for r in rs]
        its = [float(_mean_iters(v)) for v in its if v is not None]
        iters.append(float(np.mean(its)) if its else float("nan"))
        per_cell.append(rs)
    return np.array(n), np.array(t), np.array(iters), per_cell


def _mean_iters(v):
    """`linear_solver_iterations` may be a number or {forward: {mean}, adjoint: {mean}}."""
    if isinstance(v, (int, float)):
        return v
    if isinstance(v, dict):
        vals = []
        for part in ("forward", "adjoint"):
            x = v.get(part)
            if isinstance(x, dict) and "mean" in x:
                vals.append(float(x["mean"]))
            elif isinstance(x, (int, float)):
                vals.append(float(x))
        if vals:
            return sum(vals) / len(vals)
        if "mean" in v:
            return float(v["mean"])
    return float("nan")


def bootstrap_series(per_cell, metric, rng):
    """One bootstrap replicate of the per-cell metric values."""
    fn = METRICS[metric][1]
    out = []
    for rs in per_cell:
        if metric in ("cold", "warm_2pct"):
            # Resample individual evaluation timings within each run, then runs.
            picks = [rs[i] for i in rng.integers(0, len(rs), len(rs))]
            vals = []
            for r in picks:
                s = eval_samples(r)
                vals.append(float(np.median(rng.choice(s, len(s), replace=True))))
            out.append(float(np.median(vals)))
        else:
            picks = [rs[i] for i in rng.integers(0, len(rs), len(rs))]
            out.append(float(np.median([fn(r) for r in picks])))
    return np.array(out)


# ── analysis ───────────────────────────────────────────────────────────────


def analyse(runs: list[dict], n_boot: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    result: dict = {"threads": {}, "fixture_classes": {}}
    by_thread = group_by(runs, "threads")
    for (threads,), t_runs in sorted(by_thread.items(), key=lambda kv: kv[0][0] or 0):
        thread_block: dict = {}
        by_fixture = group_by(t_runs, "fixture")
        for (fixture,), f_runs in sorted(by_fixture.items(), key=lambda kv: fixture_sort_key(kv[0][0])):
            direct_runs = [r for r in f_runs if _is_direct(r["backend"])]
            iter_backends = sorted({r["backend"] for r in f_runs if not _is_direct(r["backend"])})
            fixture_block: dict = {}
            for metric in METRICS:
                entry: dict = {"n_edges": None, "ci90": None, "direct_fit": None, "iterative_fits": {}}
                if len(direct_runs) >= 3:
                    n, t, _, cells = cell_series(direct_runs, metric)
                    fit = fit_direct(n, t)
                    fit["edges"] = [int(v) for v in n]
                    fit["observed_ms"] = [float(v) for v in t]
                    fit["fitted_ms"] = [float(v) for v in direct_model(fit["params"], n)]
                    entry["direct_fit"] = fit
                    direct_cells = cells
                    direct_n = n
                else:
                    direct_cells = None
                    direct_n = None
                for backend in iter_backends:
                    b_runs = [r for r in f_runs if r["backend"] == backend]
                    if len(b_runs) < 3:
                        continue
                    n, t, iters, cells = cell_series(b_runs, metric)
                    if np.isnan(iters).any():
                        iters = np.where(np.isnan(iters), 1.0, iters)
                    iters_coef = fit_iters(n, iters)
                    fit = fit_iterative(n, t, iters)
                    fit["iters_model"] = iters_coef
                    fit["edges"] = [int(v) for v in n]
                    fit["observed_ms"] = [float(v) for v in t]
                    fit["fitted_ms"] = [float(v) for v in iterative_model(fit["params"], n, iters)]
                    if entry["direct_fit"] is not None:
                        x = crossover(entry["direct_fit"]["params"], fit["params"], iters_coef)
                        fit["crossover_edges"] = x
                        if x is not None and n_boot > 0:
                            xs = []
                            for _ in range(n_boot):
                                td = bootstrap_series(direct_cells, metric, rng)
                                ti = bootstrap_series(cells, metric, rng)
                                fd = fit_direct(direct_n, td)["params"]
                                fi = fit_iterative(n, ti, iters)["params"]
                                xb = crossover(fd, fi, iters_coef)
                                if xb is not None:
                                    xs.append(xb)
                            if xs:
                                fit["crossover_ci90"] = [float(np.percentile(xs, 5)), float(np.percentile(xs, 95))]
                        # The headline number is the CPU iterative backend when present.
                        if entry["n_edges"] is None or "cpu" in backend:
                            entry["n_edges"] = x
                            entry["ci90"] = fit.get("crossover_ci90")
                            entry["iterative_backend"] = backend
                    entry["iterative_fits"][backend] = fit
                fixture_block[metric] = entry
            thread_block[fixture] = fixture_block
        result["threads"][str(threads)] = thread_block
    # Headline per fixture class: the largest thread setting (the production one).
    if result["threads"]:
        top = max(result["threads"], key=lambda k: int(k) if k.isdigit() else 0)
        for fixture, block in result["threads"][top].items():
            result["fixture_classes"][fixture] = {
                m: {"n_edges": e["n_edges"], "ci90": e["ci90"], "threads": int(top) if top.isdigit() else top}
                for m, e in block.items()
            }
    return result


# ── rendering ──────────────────────────────────────────────────────────────


def _fmt_params_direct(p):
    a, b, c = p
    return f"a = {a:.3e}, b = {b:.3e}, c = {c:.3e}"


def _fmt_params_iter(fit):
    d, e, f = fit["params"]
    i0, i1 = fit["iters_model"]
    return f"d = {d:.3e}, e = {e:.3e}, f = {f:.3e}; iters = {i0:.2f} + {i1:.3f}·log n"


def render_markdown(result: dict, machine_id: str, source: str, plot_name: str | None) -> str:
    lines = [f"# Crossover analysis — `{machine_id}`", "",
             f"Source: `{source}`. Models of §5.3 (Stage B) fitted by least squares on "
             "log-scaled data; residuals are relative errors of the fit against the "
             "per-cell medians. Direct: `t = a·n^1.5 + b·n·log n + c`; iterative: "
             "`t = (d + e·iters)·n + f`. Times in ms, `n` = edges.", ""]
    have_iterative = any(
        e["iterative_fits"] for tb in result["threads"].values() for fb in tb.values() for e in fb.values())
    if not have_iterative:
        lines += ["Only the `direct` backend is present in this data set, so no crossover can be "
                  "solved; the direct cost model and its residuals are reported for each fixture "
                  "class and thread setting. `warm_2pct` equals `cold` for the direct backend "
                  "(a refactorisation does not warm-start).", ""]
    lines += ["## Crossover (headline)", "", "| fixture | threads | metric | crossover edges | 90% interval | iterative backend |",
              "|---|---:|---|---:|---|---|"]
    for fixture, block in result["fixture_classes"].items():
        for metric, e in block.items():
            n_edges = f"{e['n_edges']:,.0f}" if e["n_edges"] else "n/a (no iterative data)"
            ci = f"[{e['ci90'][0]:,.0f}, {e['ci90'][1]:,.0f}]" if e.get("ci90") else "—"
            backend = result["threads"][str(e["threads"])][fixture][metric].get("iterative_backend", "—")
            lines.append(f"| {fixture} | {e['threads']} | {metric} | {n_edges} | {ci} | {backend} |")
    lines.append("")
    if plot_name:
        lines += [f"![fits]({plot_name})", ""]
    for threads, tb in result["threads"].items():
        lines += [f"## Fits at {threads} thread(s)", ""]
        for fixture, fb in tb.items():
            lines += [f"### {fixture}", ""]
            for metric, e in fb.items():
                desc = METRICS[metric][0]
                lines += [f"**{metric}** — {desc}", ""]
                if metric == "warm_2pct" and not e["iterative_fits"]:
                    lines += ["Direct only: identical to `cold` (no warm start in a refactorisation).", ""]
                    continue
                fit = e["direct_fit"]
                if fit is None:
                    lines += ["direct: fewer than three sizes, not fitted.", ""]
                else:
                    ok = "ok" if fit["max_rel"] < 0.15 else "exceeds 15%"
                    lines += [f"direct fit: {_fmt_params_direct(fit['params'])}; "
                              f"RMS relative residual {100 * fit['rms_rel']:.1f}%, max {100 * fit['max_rel']:.1f}% ({ok}).", "",
                              "| edges | observed ms | fitted ms | residual |", "|---:|---:|---:|---:|"]
                    for n, o, f_, r in zip(fit["edges"], fit["observed_ms"], fit["fitted_ms"], fit["rel_residuals"]):
                        lines.append(f"| {n:,} | {fmt_ms(o)} | {fmt_ms(f_)} | {100 * r:+.1f}% |")
                    lines.append("")
                for backend, ifit in e["iterative_fits"].items():
                    x = ifit.get("crossover_edges")
                    ci = ifit.get("crossover_ci90")
                    xs = f"{x:,.0f}" if x else "none below 1e9 edges"
                    cis = f" (90%: [{ci[0]:,.0f}, {ci[1]:,.0f}])" if ci else ""
                    lines += [f"`{backend}` fit: {_fmt_params_iter(ifit)}; RMS relative residual "
                              f"{100 * ifit['rms_rel']:.1f}%, max {100 * ifit['max_rel']:.1f}%. Crossover: {xs}{cis}.", ""]
    return "\n".join(lines) + "\n"


def plot_fits(result: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    metrics = ["cold", "solve40"]
    threads = list(result["threads"].keys())
    fig, axes = plt.subplots(1, len(metrics), figsize=(6 * len(metrics), 4.5), squeeze=False)
    colours = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for ax, metric in zip(axes[0], metrics):
        k = 0
        for t in threads:
            for fixture, fb in result["threads"][t].items():
                e = fb[metric]
                fit = e["direct_fit"]
                if fit is None:
                    continue
                c = colours[k % len(colours)]
                k += 1
                n = np.array(fit["edges"], dtype=float)
                ax.loglog(n, fit["observed_ms"], "o", color=c, label=f"direct {fixture} {t}t")
                grid = np.logspace(math.log10(n.min() / 1.5), math.log10(n.max() * 1.5), 100)
                ax.loglog(grid, direct_model(fit["params"], grid), "-", color=c, alpha=0.8)
                for backend, ifit in e["iterative_fits"].items():
                    c = colours[k % len(colours)]
                    k += 1
                    ni = np.array(ifit["edges"], dtype=float)
                    ax.loglog(ni, ifit["observed_ms"], "s", color=c, label=f"{backend} {fixture} {t}t")
                    i0, i1 = ifit["iters_model"]
                    ax.loglog(grid, iterative_model(ifit["params"], grid, i0 + i1 * np.log(grid)), "--", color=c, alpha=0.8)
        ax.set_xlabel("edges")
        ax.set_ylabel("ms")
        ax.set_title(f"{metric}: data and fitted cost models")
        ax.grid(True, which="both", alpha=0.3)
        ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, format="svg")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("runs", type=Path, help="runs.jsonl produced by bench_sweep.py")
    parser.add_argument("--out", type=Path, default=None,
                        help="report directory (default: <repo>/benchmarks/reports)")
    parser.add_argument("--bootstrap", type=int, default=500, help="bootstrap replicates for the 90%% interval")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("--machine-id", default=None, help="override the id recorded in the runs")
    args = parser.parse_args()

    runs = load_runs(args.runs)
    if not runs:
        print(f"no successful bench_scale records in {args.runs}", file=sys.stderr)
        return 1
    machine_id = args.machine_id or runs[0].get("machine_id") or "unknown-machine"
    root = Path(__file__).resolve().parent.parent
    out_dir = args.out or root / "benchmarks" / "reports"
    out_dir.mkdir(parents=True, exist_ok=True)

    result = analyse(runs, args.bootstrap, args.seed)
    result["machine_id"] = machine_id
    result["source"] = str(args.runs)
    result["models"] = {"direct": "a*n^1.5 + b*n*log(n) + c", "iterative": "(d + e*iters)*n + f, iters = i0 + i1*log(n)"}
    result["metrics"] = {k: v[0] for k, v in METRICS.items()}

    plot_name = None
    if not args.no_plot:
        plot_name = f"crossover-{machine_id}.svg"
        plot_fits(result, out_dir / plot_name)
    try:
        source = str(args.runs.resolve().relative_to(root))
    except ValueError:
        source = str(args.runs)
    md = render_markdown(result, machine_id, source, plot_name)
    (out_dir / f"crossover-{machine_id}.md").write_text(md, encoding="utf-8")
    (out_dir / f"crossover-{machine_id}.json").write_text(json.dumps(result, indent=1), encoding="utf-8")
    print(md)
    print(f"wrote {out_dir / f'crossover-{machine_id}.md'} and .json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
