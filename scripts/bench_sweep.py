#!/usr/bin/env python3
"""Drive `tests/bench_scale.rs` over a parameter grid, collect its JSON lines
and render a Markdown report with SVG plots (§5.1–5.2 of
`crates/theseus/ITERATIVE_SOLVER_PROGRAM.md`).

Usage (from the repository root; `RUSTUP_TOOLCHAIN=1.89.0` in the environment
is honoured by cargo):

    uv run --project scripts scripts/bench_sweep.py \\
        --sizes 72,160,224,320,448 --fixtures grid,irregular,dome \\
        --threads 1,4 --reps 5 --iters 10 --solvers direct

    uv run --project scripts scripts/bench_sweep.py --config sweep.json
    uv run --project scripts scripts/bench_sweep.py --report-only \\
        benchmarks/results/<machine-id>/<YYYYMMDD>-<sha7>

Config file (JSON; every key optional, CLI flags override):

    {"sizes": [72, 160, 224], "fixtures": ["grid", "irregular", "dome"],
     "threads": [1, 4], "reps": 5, "iters": 10, "solvers": ["direct"],
     "wall_cap_s": 600, "warmup": false}

Each cell (solver × fixture × threads × size) runs in its own process so peak
RSS and thread pools are per cell. The harness discards one evaluation before
timing `reps` evaluations (median + IQR); `--warmup` additionally discards a
whole run per cell (§5.2). Cells exceeding `wall_cap_s` or failing are
recorded in `runs.jsonl` with a `status` field.

Outputs:

    benchmarks/results/<machine-id>/<YYYYMMDD>-<sha7>/runs.jsonl   one JSON object per cell
    benchmarks/results/<machine-id>/<YYYYMMDD>-<sha7>/machine.json OS, CPU, cores, RAM, GPUs, Rust
    benchmarks/results/<machine-id>/<YYYYMMDD>-<sha7>/config.json  the sweep configuration
    benchmarks/results/<machine-id>/<YYYYMMDD>-<sha7>/harness-output.txt  harness tables
    benchmarks/reports/<machine-id>-<YYYYMMDD>.md                  tables (median ± IQR) and fit
    benchmarks/reports/<machine-id>-<YYYYMMDD>-{eval,total}.svg    log-log time vs edges
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import crossover_fit  # noqa: E402
from bench_common import (  # noqa: E402
    eval_iqr,
    eval_median,
    fixture_sort_key,
    fmt_int,
    fmt_ms,
    fmt_pm,
    git_sha,
    group_by,
    load_failures,
    load_runs,
    machine_info,
)

DEFAULTS = {
    "sizes": [72, 160, 224],
    "fixtures": ["grid"],
    "threads": [1],
    "reps": 5,
    "iters": 10,
    "solvers": ["direct"],
    "wall_cap_s": 600,
    "warmup": False,
}
KNOWN_FIXTURES = ["grid", "irregular", "dome", "few-supports", "anisotropic"]


def parse_list(text: str, cast):
    return [cast(x.strip()) for x in text.split(",") if x.strip()]


def build_harness(root: Path, env: dict) -> Path:
    """Compiles the bench_scale test binary once and returns its path."""
    cmd = ["cargo", "test", "--locked", "--release", "-p", "theseus", "--test", "bench_scale",
           "--no-run", "--message-format=json"]
    print("building", " ".join(cmd), flush=True)
    out = subprocess.run(cmd, cwd=root, env=env, text=True, capture_output=True)
    if out.returncode:
        sys.stderr.write(out.stderr)
        raise SystemExit("cargo build failed")
    for line in out.stdout.splitlines():
        try:
            msg = json.loads(line)
        except json.JSONDecodeError:
            continue
        if msg.get("executable") and msg.get("target", {}).get("name") == "bench_scale":
            return Path(msg["executable"])
    raise SystemExit("bench_scale executable not found in cargo output")


def run_cell(binary: Path, root: Path, base_env: dict, cell: dict, runs_path: Path, log, wall_cap: float,
             machine_id: str) -> str:
    env = dict(base_env)
    env.update({
        "THESEUS_FIXTURE": cell["fixture"],
        "THESEUS_SCALE_GRIDS": str(cell["size"]),
        "THESEUS_SCALE_ITERS": str(cell["iters"]),
        "THESEUS_BENCH_REPS": str(cell["reps"]),
        "THESEUS_LINEAR_SOLVER": cell["solver"],
        "THESEUS_BENCH_JSON": str(runs_path),
        "THESEUS_MACHINE_ID": machine_id,
        "RAYON_NUM_THREADS": str(cell["threads"]),
    })
    cmd = [str(binary), "--ignored", "--nocapture", "--test-threads=1", "bench_scale_grid_solves"]
    t0 = time.time()
    status = "ok"
    try:
        out = subprocess.run(cmd, cwd=root, env=env, text=True, capture_output=True, timeout=wall_cap)
        text = out.stdout + out.stderr
        if out.returncode:
            status = "failed"
        elif "not yet available" in text:
            status = "unavailable"
    except subprocess.TimeoutExpired as exc:
        parts = [p.decode(errors="replace") if isinstance(p, bytes) else (p or "") for p in (exc.stdout, exc.stderr)]
        text = "".join(parts)
        status = "timeout"
    elapsed = time.time() - t0
    log.write(f"### {cell['solver']} {cell['fixture']} threads={cell['threads']} size={cell['size']} "
              f"status={status} wall={elapsed:.1f}s\n{text}\n")
    log.flush()
    if status != "ok":
        record = {
            "schema": 1, "harness": "bench_scale", "status": status,
            "timestamp_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            "fixture": cell["fixture"], "grid_side": cell["size"], "backend": cell["solver"],
            "threads": cell["threads"], "wall_s": elapsed, "wall_cap_s": wall_cap,
            "parameters": {"iterations_budget": cell["iters"], "eval_reps": cell["reps"]},
            "machine_id": machine_id,
        }
        with open(runs_path, "a", encoding="utf-8") as fh:
            fh.write(json.dumps(record) + "\n")
    return status


# ── report ─────────────────────────────────────────────────────────────────


def render_tables(runs: list[dict], failures: list[dict]) -> list[str]:
    lines: list[str] = []
    for (backend, threads), group in sorted(group_by(runs, "backend", "threads").items(),
                                            key=lambda kv: (kv[0][0], kv[0][1] or 0)):
        lines += [f"### `{backend}`, {threads} thread(s)", "",
                  "| fixture | grid | edges | free nodes | setup ms | eval ms (median ± IQR) | evals | iters | total ms | ms/iter | peak RSS MB | final loss |",
                  "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
        for r in sorted(group, key=lambda r: (fixture_sort_key(r["fixture"]), r["edges"])):
            rss = r.get("peak_rss_bytes")
            rss_mb = f"{rss / 2**20:.0f}" if rss is not None else "n/a"
            ms_per_iter = r.get("ms_per_iteration", r["total_ms"] / max(r["iterations"], 1))
            lines.append(
                f"| {r['fixture']} | {r.get('grid_side', '')} | {fmt_int(r['edges'])} | {fmt_int(r['free_nodes'])} | "
                f"{fmt_ms(r['setup_ms'])} | {fmt_pm(eval_median(r), eval_iqr(r))} | {r['evaluations']} | "
                f"{r['iterations']} | {fmt_ms(r['total_ms'])} | {fmt_ms(ms_per_iter)} | {rss_mb} | {r['final_loss']:.4e} |")
        lines.append("")
    if failures:
        lines += ["### Cells that did not complete", "", "| backend | fixture | grid | threads | status | wall s |",
                  "|---|---|---:|---:|---|---:|"]
        for f in failures:
            lines.append(f"| {f.get('backend')} | {f.get('fixture')} | {f.get('grid_side')} | {f.get('threads')} | "
                         f"{f.get('status')} | {f.get('wall_s', 0):.0f} |")
        lines.append("")
    return lines


def plot_metric(runs: list[dict], path: Path, metric: str, title: str) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.8))
    series = group_by(runs, "backend", "fixture", "threads")
    markers = {1: "o", 2: "^", 4: "s", 8: "D"}
    for (backend, fixture, threads), group in sorted(series.items(),
                                                     key=lambda kv: (kv[0][0], fixture_sort_key(kv[0][1]), kv[0][2] or 0)):
        group = sorted(group, key=lambda r: r["edges"])
        x = [r["edges"] for r in group]
        if metric == "eval":
            y = [eval_median(r) for r in group]
            err = [[max(m - r["eval_ms"].get("q25", m), 0) for m, r in zip(y, group)],
                   [max(r["eval_ms"].get("q75", m) - m, 0) for m, r in zip(y, group)]]
            ax.errorbar(x, y, yerr=err, marker=markers.get(threads, "x"), capsize=2, linestyle="-",
                        label=f"{backend} {fixture} {threads}t")
        else:
            y = [r["total_ms"] for r in group]
            ax.plot(x, y, marker=markers.get(threads, "x"), linestyle="-", label=f"{backend} {fixture} {threads}t")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("edges")
    ax.set_ylabel("ms")
    ax.set_title(title)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(path, format="svg")
    plt.close(fig)


def render_report(results_dir: Path, reports_dir: Path, root: Path, fit: bool = True) -> Path:
    runs_path = results_dir / "runs.jsonl"
    runs = load_runs(runs_path)
    failures = load_failures(runs_path)
    machine = json.loads((results_dir / "machine.json").read_text(encoding="utf-8"))
    config = json.loads((results_dir / "config.json").read_text(encoding="utf-8")) if (results_dir / "config.json").exists() else {}
    machine_id = machine["machine_id"]
    date = results_dir.name.split("-")[0]
    stem = f"{machine_id}-{date}"
    reports_dir.mkdir(parents=True, exist_ok=True)
    try:
        rel_results = results_dir.resolve().relative_to(root)
    except ValueError:
        rel_results = results_dir

    lines = [f"# Benchmark sweep — `{machine_id}` — {date}", "",
             f"Raw data: `{rel_results}/runs.jsonl` ({len(runs)} completed cells, {len(failures)} not completed). "
             f"Git: `{machine.get('git_sha', 'unknown')}`.", "",
             "## Machine", "",
             f"* OS: {machine.get('os_release')}",
             f"* CPU: {machine.get('cpu_model')} — {machine.get('logical_cores')} logical / {machine.get('physical_cores')} physical cores",
             f"* RAM: {machine.get('ram_gib')} GiB",
             f"* GPU adapters: {', '.join(machine.get('gpu_adapters') or []) or 'none discovered'}",
             f"* Rust: {machine.get('rust_version')}", ""]
    if config:
        lines += ["## Configuration", "", "```json", json.dumps(config, indent=1), "```", ""]
    lines += ["## Protocol", "",
              "Release build (`--locked`). One process per cell. The harness performs the cache setup and first "
              "factorisation (`setup ms`), discards one fused evaluation, times `reps` evaluations "
              "(`eval ms`: median ± interquartile range; ` !` flags cells whose IQR exceeds 10% of the median), "
              "then runs one full L-BFGS-B solve with a fixed iteration budget and tolerances set to zero "
              "(`total ms`, box bounds `[0.1, 10]`, start `q = 1`). Non-grid fixtures are generated at the "
              "edge count of the grid of the same `grid` column. `threads` is `RAYON_NUM_THREADS`. "
              "Peak RSS is `VmHWM` of the cell's process.", ""]
    if runs:
        lines += ["## Results", ""] + render_tables(runs, failures)
        eval_svg = f"{stem}-eval.svg"
        total_svg = f"{stem}-total.svg"
        plot_metric(runs, reports_dir / eval_svg, "eval", "Fused evaluation time vs edges")
        plot_metric(runs, reports_dir / total_svg, "total", "Full solve time vs edges (fixed iteration budget)")
        lines += ["## Plots", "", f"![evaluation time]({eval_svg})", "", f"![solve time]({total_svg})", ""]
        if fit:
            result = crossover_fit.analyse(runs, n_boot=0, seed=0)
            lines += ["## Cost-model fit", "",
                      "Direct model `t = a·n^1.5 + b·n·log n + c` fitted by least squares on log-scaled data, "
                      "per fixture and thread setting, to the evaluation time (`cold`) and to the full solve "
                      "extrapolated to 40 iterations (`solve40`). Full details, the iterative models and the "
                      f"crossover search are in `crossover-{machine_id}.md` (from `scripts/crossover_fit.py`).", "",
                      "| threads | fixture | metric | a | b | c | RMS residual | max residual |",
                      "|---:|---|---|---:|---:|---:|---:|---:|"]
            for threads, tb in result["threads"].items():
                for fixture, fb in tb.items():
                    for metric in ("cold", "solve40"):
                        f = fb[metric]["direct_fit"]
                        if f is None:
                            continue
                        a, b, c = f["params"]
                        flag = "" if f["max_rel"] < 0.15 else " (> 15%)"
                        lines.append(f"| {threads} | {fixture} | {metric} | {a:.3e} | {b:.3e} | {c:.3e} | "
                                     f"{100 * f['rms_rel']:.1f}% | {100 * f['max_rel']:.1f}%{flag} |")
            lines.append("")
    else:
        lines += ["No completed cells.", ""]
    report = reports_dir / f"{stem}.md"
    report.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report


# ── main ───────────────────────────────────────────────────────────────────


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", type=Path, help="JSON config file")
    parser.add_argument("--sizes", help="comma-separated grid sides")
    parser.add_argument("--fixtures", help=f"comma-separated subset of {','.join(KNOWN_FIXTURES)}")
    parser.add_argument("--threads", help="comma-separated RAYON_NUM_THREADS values")
    parser.add_argument("--reps", type=int, help="timed evaluations per cell (>= 5 recommended)")
    parser.add_argument("--iters", type=int, help="L-BFGS-B iteration budget")
    parser.add_argument("--solvers", help="comma-separated THESEUS_LINEAR_SOLVER values")
    parser.add_argument("--wall-cap", type=float, dest="wall_cap_s", help="per-cell wall-clock cap in seconds")
    parser.add_argument("--warmup", action="store_true", default=None, help="discard one full run per cell first")
    parser.add_argument("--out", type=Path, default=None, help="benchmarks directory (default <repo>/benchmarks)")
    parser.add_argument("--label", default=None, help="results sub-directory name (default <YYYYMMDD>-<sha7>)")
    parser.add_argument("--report-only", type=Path, default=None, help="re-render the report for a results directory")
    parser.add_argument("--no-fit", action="store_true", help="skip the cost-model section of the report")
    parser.add_argument("--dry-run", action="store_true", help="print the cells and exit")
    args = parser.parse_args()

    root = Path(__file__).resolve().parent.parent
    bench_dir = args.out or root / "benchmarks"
    reports_dir = bench_dir / "reports"

    if args.report_only:
        report = render_report(args.report_only, reports_dir, root, fit=not args.no_fit)
        print(f"wrote {report}")
        return 0

    cfg = dict(DEFAULTS)
    if args.config:
        cfg.update(json.loads(args.config.read_text(encoding="utf-8")))
    if args.sizes:
        cfg["sizes"] = parse_list(args.sizes, int)
    if args.fixtures:
        cfg["fixtures"] = parse_list(args.fixtures, str)
    if args.threads:
        cfg["threads"] = parse_list(args.threads, int)
    if args.reps is not None:
        cfg["reps"] = args.reps
    if args.iters is not None:
        cfg["iters"] = args.iters
    if args.solvers:
        cfg["solvers"] = parse_list(args.solvers, str)
    if args.wall_cap_s is not None:
        cfg["wall_cap_s"] = args.wall_cap_s
    if args.warmup is not None:
        cfg["warmup"] = args.warmup
    for f in cfg["fixtures"]:
        if f not in KNOWN_FIXTURES:
            raise SystemExit(f"unknown fixture {f!r}; known: {KNOWN_FIXTURES}")

    cells = [{"solver": s, "fixture": f, "threads": t, "size": n, "iters": cfg["iters"], "reps": cfg["reps"]}
             for s in cfg["solvers"] for f in cfg["fixtures"] for t in cfg["threads"] for n in sorted(cfg["sizes"])]
    print(f"{len(cells)} cells: solvers={cfg['solvers']} fixtures={cfg['fixtures']} threads={cfg['threads']} "
          f"sizes={sorted(cfg['sizes'])} reps={cfg['reps']} iters={cfg['iters']} warmup={cfg['warmup']}")
    if args.dry_run:
        for c in cells:
            print(c)
        return 0

    env = dict(os.environ)
    info = machine_info(root)
    machine_id = info["machine_id"]
    sha = git_sha(root)
    label = args.label or f"{dt.datetime.now(dt.timezone.utc).strftime('%Y%m%d')}-{sha[:7]}"
    results_dir = bench_dir / "results" / machine_id / label
    results_dir.mkdir(parents=True, exist_ok=True)
    info["cpu_governor"] = _cpu_governor()
    info["swept_at_utc"] = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    (results_dir / "machine.json").write_text(json.dumps(info, indent=1) + "\n", encoding="utf-8")
    (results_dir / "config.json").write_text(json.dumps(cfg, indent=1) + "\n", encoding="utf-8")
    print(f"machine id {machine_id}; results in {results_dir}")

    binary = build_harness(root, env)
    runs_path = results_dir / "runs.jsonl"
    t_start = time.time()
    unavailable: set[str] = set()
    with open(results_dir / "harness-output.txt", "a", encoding="utf-8") as log:
        for i, cell in enumerate(cells, 1):
            desc = f"[{i}/{len(cells)}] {cell['solver']} {cell['fixture']} threads={cell['threads']} size={cell['size']}"
            if cell["solver"] in unavailable:
                print(f"{desc}: skipped (backend not yet available)", flush=True)
                continue
            if cfg["warmup"]:
                run_cell(binary, root, env, cell, results_dir / "warmup.jsonl", log, cfg["wall_cap_s"], machine_id)
            t0 = time.time()
            status = run_cell(binary, root, env, cell, runs_path, log, cfg["wall_cap_s"], machine_id)
            print(f"{desc}: {status} in {time.time() - t0:.1f}s (elapsed {time.time() - t_start:.0f}s)", flush=True)
            if status == "unavailable":
                unavailable.add(cell["solver"])
    report = render_report(results_dir, reports_dir, root, fit=not args.no_fit)
    print(f"wrote {report}")
    return 0


def _cpu_governor() -> str | None:
    path = Path("/sys/devices/system/cpu/cpu0/cpufreq/scaling_governor")
    try:
        return path.read_text().strip()
    except OSError:
        return None


if __name__ == "__main__":
    sys.exit(main())
