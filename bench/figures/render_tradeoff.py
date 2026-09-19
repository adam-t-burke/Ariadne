"""Render GN-versus-L-BFGS-B tradeoff figures from `tradeoff.json`.

    cargo run --release -p theseus --example warm_start_bench -- tradeoff bench/figures/data
    cd bench/figures && uv run render_tradeoff.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

INVERSE = {
    "s1": ("Stage 1 (q*)", "#ff7f0e"),
    "frozen1": ("1 frozen CWLS", "#2ca02c"),
    "frozen2": ("2 frozen CWLS", "#98df8a"),
    "gn1": ("1 GN (no frozen)", "#9467bd"),
    "gn2": ("2 GN (no frozen)", "#c5b0d5"),
    "gn3": ("3 GN (no frozen)", "#7b4173"),
    "pipe": ("1 frozen + 2 GN", "#1f5fbf"),
}

LBFGS = {
    "s1+lb10": ("q* + 10 L-BFGS-B", "#d62728"),
    "frozen1+lb10": ("frozen + 10 L-BFGS-B", "#17becf"),
    "gn1+lb10": ("GN1 + 10 L-BFGS-B", "#e377c2"),
}


def load(path: Path) -> dict:
    return json.loads(path.read_text())


def case_label(c: dict) -> str:
    return f"{c['net']}\n{c['target']} {c['box']}"


def row_map(c: dict) -> dict:
    return {r["label"]: r for r in c.get("rows", [])}


def geom(row: dict | None) -> float:
    if not row:
        return np.nan
    v = row.get("geom_over_l")
    return float(v) if v is not None else np.nan


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "figure.dpi": 140,
            "savefig.dpi": 160,
            "savefig.bbox": "tight",
        }
    )


def fig_summary_bars(data: dict, out: Path) -> None:
    cases = data["cases"]
    labels = [case_label(c) for c in cases]
    series = [
        ("s1", "Stage 1", "#ff7f0e"),
        ("frozen1", "1 frozen", "#2ca02c"),
        ("gn2", "2 GN, no frozen", "#9467bd"),
        ("pipe", "frozen+2 GN", "#1f5fbf"),
        ("s1+lb10", "q* + 10 L-BFGS", "#d62728"),
        ("frozen1+lb10", "frozen + 10 L-BFGS", "#17becf"),
    ]
    x = np.arange(len(cases))
    width = 0.14
    fig, ax = plt.subplots(figsize=(12.5, 5.2))
    for i, (key, name, color) in enumerate(series):
        vals = []
        for c in cases:
            rows = row_map(c)
            vals.append(geom(rows.get(key)))
        ax.bar(x + (i - 2.5) * width, vals, width, label=name, color=color, log=True)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=0, fontsize=7)
    ax.set_ylabel("geometric error / L")
    ax.set_title("Warm-start geometric error: linearised CWLS / GN vs 10 L-BFGS-B steps from q*")
    ax.legend(ncol=3, fontsize=8, loc="upper right")
    ax.set_yscale("log")
    fig.savefig(out)
    plt.close(fig)


def fig_pareto(data: dict, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    markers = {
        "s1": ("o", "#ff7f0e"),
        "frozen1": ("s", "#2ca02c"),
        "gn2": ("D", "#9467bd"),
        "pipe": ("^", "#1f5fbf"),
        "s1+lb10": ("x", "#d62728"),
        "frozen1+lb10": ("+", "#17becf"),
    }
    for key, (mk, color) in markers.items():
        xs, ys = [], []
        for c in data["cases"]:
            r = row_map(c).get(key)
            if not r:
                continue
            xs.append(r.get("ms") or np.nan)
            ys.append(geom(r))
        ax.scatter(xs, ys, marker=mk, c=color, label=key, s=36)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("wall time (ms)")
    ax.set_ylabel("geometric error / L")
    ax.set_title("Cost vs geometric error (one marker per case)")
    ax.legend(fontsize=8)
    fig.savefig(out)
    plt.close(fig)


def fig_traces(data: dict, out: Path) -> None:
    interesting = [
        "quad21c",
        "cabledome4x16",
        "hypar21m",
        "cabletruss16",
        "quad21c_d1",
        "barrel16x12",
    ]
    picked = []
    for name in interesting:
        for c in data["cases"]:
            if c["net"] == name and c not in picked:
                picked.append(c)
                break
    n = len(picked)
    cols = 3
    rows = (n + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(11.5, 3.2 * rows), squeeze=False)
    for ax, c in zip(axes.ravel(), picked):
        rm = row_map(c)
        l = c["extent"]
        for key, color, ls in [
            ("s1+lb10", "#d62728", "-"),
            ("s1+lb20", "#d62728", ":"),
            ("frozen1+lb10", "#17becf", "-"),
            ("gn1+lb10", "#e377c2", "--"),
        ]:
            r = rm.get(key)
            if not r:
                continue
            tr = [v for v in (r.get("trace_geom_over_l") or []) if v is not None]
            if tr:
                ax.plot(np.arange(len(tr)), tr, color=color, ls=ls, lw=1.2, label=key)
        for key, color in [("s1", "#ff7f0e"), ("frozen1", "#2ca02c"), ("gn2", "#9467bd"), ("pipe", "#1f5fbf")]:
            r = rm.get(key)
            if r and r.get("geom_over_l") is not None:
                ax.axhline(r["geom_over_l"], color=color, ls="--", lw=0.9, alpha=0.85)
        ax.set_yscale("log")
        ax.set_title(f"{c['net']} {c['target']} {c['box']}", fontsize=8)
        ax.set_xlabel("L-BFGS-B evaluation")
        ax.set_ylabel("geom / L")
    axes[0, 0].legend(fontsize=6, loc="upper right")
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.suptitle("L-BFGS-B evaluation traces vs linearised warm-start levels", y=1.01)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    _ = l


def fig_force_vs_geom(data: dict, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    for key, (name, color) in {**INVERSE, **LBFGS}.items():
        xs, ys = [], []
        for c in data["cases"]:
            r = row_map(c).get(key)
            if not r:
                continue
            if r.get("force_rel") is None or r.get("geom_over_l") is None:
                continue
            xs.append(r["force_rel"])
            ys.append(r["geom_over_l"])
        ax.scatter(xs, ys, c=color, label=name, s=28, alpha=0.85)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("force residual  ‖E(x*)q − p‖ / ‖p‖")
    ax.set_ylabel("geometric error / L")
    ax.set_title("Metric mismatch: force residual vs geometric error")
    ax.legend(fontsize=7, loc="best")
    fig.savefig(out)
    plt.close(fig)


def fig_ratio(data: dict, out: Path) -> None:
    labels, r_fr, r_gn, r_pipe = [], [], [], []
    for c in data["cases"]:
        s = c.get("summary") or {}
        def ratio(a, b):
            va, vb = s.get(a), s.get(b)
            if va and vb and vb > 0:
                return va / vb
            return np.nan
        labels.append(f"{c['net']} {c['box']}")
        r_fr.append(ratio("s1_lb10", "frozen1"))
        r_gn.append(ratio("s1_lb10", "gn2"))
        r_pipe.append(ratio("s1_lb10", "pipe"))
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(11.5, 4.4))
    ax.bar(x - 0.25, r_fr, 0.24, label="(q*+10 L-BFGS) / frozen1", color="#d62728")
    ax.bar(x, r_gn, 0.24, label="(q*+10 L-BFGS) / 2 GN", color="#9467bd")
    ax.bar(x + 0.25, r_pipe, 0.24, label="(q*+10 L-BFGS) / pipeline", color="#1f5fbf")
    ax.axhline(1.0, color="k", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=7)
    ax.set_ylabel("geometric-error ratio  (>1 ⇒ L-BFGS worse)")
    ax.set_title("Can 10 L-BFGS-B steps from q* match a linearised warm start?")
    ax.set_yscale("log")
    ax.legend(fontsize=8)
    fig.savefig(out)
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=Path("data/tradeoff.json"))
    p.add_argument("--out", type=Path, default=Path("../../docs/figures"))
    args = p.parse_args()
    data = load(args.data)
    args.out.mkdir(parents=True, exist_ok=True)
    setup_style()
    fig_summary_bars(data, args.out / "tradeoff_geom_bars.png")
    fig_pareto(data, args.out / "tradeoff_cost_vs_error.png")
    fig_traces(data, args.out / "tradeoff_lbfgs_traces.png")
    fig_force_vs_geom(data, args.out / "tradeoff_force_vs_geom.png")
    fig_ratio(data, args.out / "tradeoff_lb10_ratio.png")
    print(f"wrote figures to {args.out}")


if __name__ == "__main__":
    main()
