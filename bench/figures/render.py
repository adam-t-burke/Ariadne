"""Render the presentation figures from the `warm_start_bench figures` export.

    cargo run --release -p theseus --example warm_start_bench -- figures bench/figures/data
    cargo run --release -p theseus --example warm_start_bench -- dense > bench/figures/data/dense.txt
    cargo run --release -p theseus --example warm_start_bench -- scale > bench/figures/data/scale.txt
    cd bench/figures && uv run render.py [--data data] [--out ../../docs/figures]

Every figure is written as PNG; `docs/presentation_outline.md` links them.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402
from mpl_toolkits.mplot3d.art3d import Line3DCollection  # noqa: E402

TENSION = "#1f5fbf"
COMPRESSION = "#c8321e"
TARGET = "#9a9a9a"
FIXED = "#222222"

METHOD_STYLE = {
    "uniform": ("uniform seed", "#7f7f7f", "-"),
    "s1": ("Stage 1 only (force residual)", "#ff7f0e", "-"),
    "gram_sparse": ("sparse Gram", "#bcbd22", "--"),
    "length_ratio": ("length-ratio heuristic", "#8c564b", ":"),
    "frozen": ("one frozen CWLS step", "#2ca02c", "-"),
    "pipeline": ("pipeline (this work)", "#1f5fbf", "-"),
    "legacy": ("previous branch (Clarabel Stage 2)", "#9467bd", "--"),
}

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "legend.fontsize": 8.5,
        "figure.dpi": 100,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
    }
)


# ───────────────────────── helpers ─────────────────────────


def load(path: Path):
    with path.open() as f:
        return json.load(f)


def to_array(rows, fill=np.nan):
    """List of rows that may contain None → float array."""
    return np.array([[fill if v is None else v for v in r] for r in rows], dtype=float)


def forces(x, edges, q):
    a = x[edges[:, 0]]
    b = x[edges[:, 1]]
    length = np.linalg.norm(b - a, axis=1)
    return q * length


def planar_axes(x):
    """Indices of the two axes to draw when the net is (nearly) planar, else None."""
    ext = x.max(0) - x.min(0)
    thin = np.argmin(ext)
    if ext[thin] < 0.05 * ext.max():
        return [i for i in range(3) if i != thin]
    return None


def draw_net(ax, x, edges, q, fixed, planar, *, target=None, width=(0.4, 3.2), alpha=1.0):
    """Edges coloured by sign of q, width by |force|; target as thin grey overlay."""
    segs = np.stack([x[edges[:, 0]], x[edges[:, 1]]], axis=1)
    if q is None:
        colors = [TARGET] * len(edges)
        lw = np.full(len(edges), 0.9)
    else:
        colors = [TENSION if v >= 0 else COMPRESSION for v in q]
        f = np.abs(forces(x, edges, q))
        fmax = np.nanmax(f) if np.isfinite(f).any() and np.nanmax(f) > 0 else 1.0
        lw = width[0] + (width[1] - width[0]) * f / fmax
    if planar is not None:
        segs2 = segs[:, :, planar]
        if target is not None:
            tsegs = np.stack([target[edges[:, 0]], target[edges[:, 1]]], axis=1)[:, :, planar]
            ax.add_collection(LineCollection(tsegs, colors=TARGET, linewidths=0.7, alpha=0.7, zorder=1))
        ax.add_collection(LineCollection(segs2, colors=colors, linewidths=lw, alpha=alpha, zorder=2))
        fx = x[fixed][:, planar]
        ax.scatter(fx[:, 0], fx[:, 1], s=14, marker="s", color=FIXED, zorder=3)
        ax.set_aspect("equal")
        lo = np.minimum(x.min(0), target.min(0) if target is not None else x.min(0))[planar]
        hi = np.maximum(x.max(0), target.max(0) if target is not None else x.max(0))[planar]
        pad = 0.05 * (hi - lo).max()
        ax.set_xlim(lo[0] - pad, hi[0] + pad)
        ax.set_ylim(lo[1] - pad, hi[1] + pad)
        ax.set_axis_off()
    else:
        if target is not None:
            tsegs = np.stack([target[edges[:, 0]], target[edges[:, 1]]], axis=1)
            ax.add_collection3d(Line3DCollection(tsegs, colors=TARGET, linewidths=0.6, alpha=0.6))
        ax.add_collection3d(Line3DCollection(segs, colors=colors, linewidths=lw, alpha=alpha))
        fx = x[fixed]
        ax.scatter(fx[:, 0], fx[:, 1], fx[:, 2], s=10, marker="s", color=FIXED, depthshade=False)
        ref = x if target is None else np.vstack([x, target])
        lo, hi = ref.min(0), ref.max(0)
        ext = np.maximum(hi - lo, 1e-9)
        ax.set_xlim(lo[0], hi[0])
        ax.set_ylim(lo[1], hi[1])
        ax.set_zlim(lo[2], hi[2])
        ax.set_box_aspect(ext / ext.max(), zoom=1.35)
        ax.view_init(elev=24, azim=-58)
        ax.set_axis_off()


def new_axes(fig, spec, planar):
    rows, cols, index = spec
    if planar is not None:
        return fig.add_subplot(rows, cols, index)
    return fig.add_subplot(rows, cols, index, projection="3d")


def err_label(err, L):
    return "–" if err is None or not np.isfinite(err) else f"{err / L:.2e}"


def sign_legend(fig, y=0.02):
    from matplotlib.lines import Line2D

    handles = [
        Line2D([], [], color=TENSION, lw=2, label="tension (q > 0)"),
        Line2D([], [], color=COMPRESSION, lw=2, label="compression (q < 0)"),
        Line2D([], [], color=TARGET, lw=1, label="target geometry"),
        Line2D([], [], color=FIXED, lw=0, marker="s", markersize=5, label="support"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, y))


# ───────────────────────── gallery ─────────────────────────


def fig_gallery(data: Path, out: Path):
    nets = load(data / "gallery.json")
    n = len(nets)
    cols = 4
    rows = int(np.ceil(n / cols))
    fig = plt.figure(figsize=(13, 2.9 * rows))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.94, bottom=0.04, wspace=0.02, hspace=0.12)
    for i, net in enumerate(nets):
        x = np.array(net["x_true"])
        edges = np.array(net["edges"])
        q = np.array(net["q_true"])
        planar = planar_axes(x)
        ax = new_axes(fig, (rows, cols, i + 1), planar)
        draw_net(ax, x, edges, q, np.array(net["fixed"]), planar)
        pos = int((q > 0).sum())
        neg = int((q < 0).sum())
        ax.set_title(f"{net['name']}  (ne={len(edges)}, {pos}+ / {neg}−)", fontsize=10, y=0.97 if planar is None else 1.0)
    sign_legend(fig, y=-0.005)
    fig.suptitle("Synthetic benchmark nets at their funicular (edge width ∝ |force|)", y=0.995)
    fig.savefig(out / "gallery.png")
    plt.close(fig)


# ───────────────────────── per-case figures ─────────────────────────


def case_title(case):
    net = case["net"]
    return f"{net['name']}  ·  target {case['target_kind']}  ·  box {case['box_kind']}"


def method(case, label):
    for m in case["methods"]:
        if m["label"] == label:
            return m
    return None


def fig_case_geometry(case, out: Path, stem: str):
    net = case["net"]
    L = net["extent"]
    edges = np.array(net["edges"])
    fixed = np.array(net["fixed"])
    target = np.array(case["x_target"])
    planar = planar_axes(target)
    iters = case["max_iters"]

    panels = [("target", "target geometry", None, None, None)]
    uni = method(case, "uniform")
    if uni is not None:
        if uni.get("x_final") is not None:
            panels.append(
                (
                    "uniform",
                    f"uniform seed → L-BFGS-B ({uni.get('lbfgsb_iters', iters)} it)",
                    to_array(uni["x_final"]),
                    to_array([uni["q_final"]])[0],
                    uni.get("final_err"),
                )
            )
        else:
            panels.append(("uniform", "uniform seed", None, None, uni.get("failed", "failed")))
    for label, title in [
        ("s1", "Stage 1 seed (force residual), clipped"),
        ("pipeline", "pipeline warm start, before L-BFGS-B"),
    ]:
        m = method(case, label)
        if m is None:
            continue
        if m.get("x_warm") is not None:
            panels.append((label, title, to_array(m["x_warm"]), to_array([m["q_warm"]])[0], m.get("err_clip")))
        else:
            panels.append((label, title, None, None, m.get("failed", "failed")))

    fig = plt.figure(figsize=(3.9 * len(panels), 3.9))
    fig.subplots_adjust(left=0.01, right=0.99, top=0.82, bottom=0.06, wspace=0.03)
    for i, (label, title, x, q, err) in enumerate(panels):
        ax = new_axes(fig, (1, len(panels), i + 1), planar)
        if label == "target":
            draw_net(ax, target, edges, None, fixed, planar)
            ax.set_title(title)
        elif x is None:
            draw_net(ax, target, edges, None, fixed, planar)
            msg = err if isinstance(err, str) else "failed"
            if "singular" in msg:
                msg = "singular Laplacian:\nforward solve fails"
            ax.text2D(0.5, 0.5, msg, transform=ax.transAxes, ha="center", va="center", fontsize=10, color=COMPRESSION) if planar is None else ax.text(
                0.5, 0.5, msg, transform=ax.transAxes, ha="center", va="center", fontsize=10, color=COMPRESSION
            )
            ax.set_title(title)
        else:
            draw_net(ax, x, edges, q, fixed, planar, target=target)
            ax.set_title(f"{title}\nerr/L = {err_label(err, L)}")
    sign_legend(fig, y=-0.03)
    fig.suptitle(case_title(case), y=0.995)
    fig.savefig(out / f"{stem}_geometry.png")
    plt.close(fig)


def fig_case_convergence(case, out: Path, stem: str, ax=None, legend=True):
    net = case["net"]
    L = net["extent"]
    own = ax is None
    if own:
        fig, ax = plt.subplots(figsize=(7.2, 5.6))
    finals = [m["final_err"] for m in case["methods"] if m.get("final_err") is not None]
    best = min(finals) if finals else None
    for m in case["methods"]:
        label = m["label"]
        if label not in METHOD_STYLE:
            continue
        name, color, ls = METHOD_STYLE[label]
        if m.get("trace") is None:
            if own:
                ax.plot([], [], color=color, ls=ls, label=f"{name}: {m.get('failed', 'failed')}")
            continue
        tr = np.sqrt(np.maximum(np.minimum.accumulate(np.array(m["trace"], dtype=float)), 0.0)) / L
        ev = np.arange(len(tr))
        warm = m.get("warm_ms", 0.0)
        lab = f"{name} (warm {warm:.0f} ms, final {m['final_err'] / L:.2e})"
        ax.plot(ev, tr, color=color, ls=ls, lw=1.5 if label == "pipeline" else 1.1, label=lab)
        ax.plot([0], [tr[0]], "o", color=color, ms=4)
    if best is not None:
        ax.axhline(best / L, color="k", lw=0.7, ls=":", label=f"best final ({best / L:.2e})")
    ax.set_yscale("log")
    ax.set_xlabel("L-BFGS-B evaluations (loss + adjoint gradient)")
    ax.set_ylabel("best geometric error so far, ‖x(q) − x*‖ / L")
    ax.grid(True, which="both", alpha=0.25)
    if legend:
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2, frameon=False)
    if own:
        ax.set_title(case_title(case))
        fig.savefig(out / f"{stem}_convergence.png")
        plt.close(fig)


def fig_convergence_grid(cases, out: Path):
    n = len(cases)
    cols = 3
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(5.0 * cols, 3.4 * rows), squeeze=False)
    for ax, (stem, case) in zip(axes.flat, cases):
        fig_case_convergence(case, out, stem, ax=ax, legend=False)
        ax.set_title(case_title(case), fontsize=9.5)
        ax.set_xlabel("")
        ax.set_ylabel("")
    for ax in axes.flat[n:]:
        ax.set_axis_off()
    from matplotlib.lines import Line2D

    handles = [Line2D([], [], color=c, ls=ls, lw=1.6, label=name) for name, c, ls in METHOD_STYLE.values()]
    handles.append(Line2D([], [], color="k", ls=":", lw=0.8, label="best final error"))
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, -0.02))
    fig.supxlabel("L-BFGS-B evaluations", y=0.03)
    fig.supylabel("geometric error / L")
    fig.tight_layout(rect=(0.02, 0.06, 1, 1))
    fig.savefig(out / "convergence_grid.png")
    plt.close(fig)


def fig_start_vs_final(cases, out: Path):
    """Clipped start error and final error after L-BFGS-B, relative to the best final, per case."""
    labels = [m for m in METHOD_STYLE if m != "legacy"]
    names = [stem.replace("case_", "").replace("_jit2pctd", " jit").replace("_bump10pctd", " bump") for stem, _ in cases]
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.2), sharey=True)
    width = 0.8 / len(labels)
    xs = np.arange(len(cases))
    for k, label in enumerate(labels):
        name, color, _ = METHOD_STYLE[label]
        start, final = [], []
        for _, case in cases:
            finals = [m["final_err"] for m in case["methods"] if m.get("final_err") is not None]
            best = min(finals) if finals else np.nan
            m = method(case, label)
            if m is None or m.get("err_clip") is None:
                start.append(np.nan)
                final.append(np.nan)
                continue
            start.append(m["err_clip"] / best)
            final.append(m["final_err"] / best if m.get("final_err") is not None else np.nan)
        for ax, vals in zip(axes, (start, final)):
            vals = np.array(vals)
            ax.bar(xs + (k - len(labels) / 2 + 0.5) * width, vals - 0.5, width, bottom=0.5, color=color, label=name)
            for xi, v in zip(xs, vals):
                if np.isnan(v):
                    ax.plot(xi + (k - len(labels) / 2 + 0.5) * width, 1.0, "x", color=color, ms=5)
    for ax, title in zip(axes, ("clipped warm start / best final", "after L-BFGS-B (1000 it) / best final")):
        ax.set_yscale("log")
        ax.axhline(1.0, color="k", lw=0.7)
        ax.axhline(1.05, color="k", lw=0.6, ls=":")
        ax.set_xticks(xs)
        ax.set_xticklabels(names, rotation=35, ha="right", fontsize=8.5)
        ax.set_title(title)
        ax.grid(True, axis="y", which="both", alpha=0.25)
    axes[0].set_ylabel("error relative to best final (log)")
    axes[0].set_ylim(0.5, None)
    handles, labels_ = axes[1].get_legend_handles_labels()
    fig.legend(handles, labels_, loc="upper center", bbox_to_anchor=(0.5, 1.08), ncol=6, frameon=False)
    fig.text(0.5, -0.22, "× = seed fails (singular Laplacian, L-BFGS-B cannot start); dotted line = 5 % above the best final", ha="center", fontsize=8.5)
    fig.savefig(out / "start_vs_final.png")
    plt.close(fig)


# ───────────────────────── reactions ─────────────────────────


def fig_reactions(data: Path, out: Path):
    d = load(data / "reactions_tiedarch16s.json")
    net = d["net"]
    edges = np.array(net["edges"])
    fixed = np.array(net["fixed"])
    tie = np.array(net["tie_edges"], dtype=int)
    L = net["extent"]
    runs = {(r["target_kind"], r["label"]): r for r in d["runs"] if r.get("x") is not None}
    order = [("exact", "rx free"), ("exact", "rx=0"), ("jit2%d", "rx free"), ("jit2%d", "rx=0")]
    fig, axes = plt.subplots(2, 2, figsize=(12, 6.2))
    for ax, key in zip(axes.flat, order):
        r = runs.get(key)
        if r is None:
            ax.set_axis_off()
            continue
        x = to_array(r["x"])
        q = to_array([r["q"]])[0]
        target = np.array(d["x_target_jit"]) if key[0] != "exact" else np.array(net["x_true"])
        planar = planar_axes(x) or [0, 2]
        draw_net(ax, x, edges, q, fixed, planar, target=target, width=(0.5, 4.0))
        R = to_array(r["reactions"]["per_support"])
        span = x[:, 0].max() - x[:, 0].min()
        scale = 0.25 * span / max(np.abs(R).max(), 1e-12)
        for i, n in enumerate(fixed):
            p = x[n][planar]
            v = -R[i][planar] * scale  # support reaction = −(member resultant)
            ax.add_patch(
                FancyArrowPatch(p, p + v, arrowstyle="-|>", mutation_scale=12, color="#111111", lw=1.4, zorder=4, clip_on=False)
            )
        q_tie = np.nanmean(q[tie]) if len(tie) else np.nan
        rx = np.abs(R[:, 0]).max()
        title = f"target {key[0]}, {key[1]}\n|Rx| = {rx:.2e}, mean q_tie = {q_tie:.2f}, err/L = {err_label(r.get('err'), L)}"
        ax.set_title(title, fontsize=9.5)
        x0, x1 = ax.get_xlim()
        ax.set_xlim(x0 - 0.28 * span, x1 + 0.28 * span)
        ax.set_ylim(ax.get_ylim()[0] - 0.12 * span, ax.get_ylim()[1] + 0.05 * span)
    fig.subplots_adjust(wspace=0.02, hspace=0.35)
    sign_legend(fig, y=-0.02)
    fig.suptitle("Straight-tie arch (tiedarch16s): support reactions with and without `enforce_zero_rx`  (arrows: support reaction, common scale per panel)", fontsize=10.5)
    fig.savefig(out / "reactions_tiedarch.png")
    plt.close(fig)


# ───────────────────────── dense / scale ─────────────────────────


def parse_dense(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 8 and parts[0].isdigit():
            rows.append((int(parts[1]), float(parts[3]), float(parts[4])))
    return np.array(rows)


def fig_dense(data: Path, out: Path):
    rows = parse_dense(data / "dense.txt")
    ne, dense, sparse = rows[:, 0], rows[:, 1], rows[:, 2]
    fig, ax = plt.subplots(figsize=(6.6, 4.3))
    ax.loglog(ne, dense, "o-", color=COMPRESSION, label="dense Gram: form EᵀE, Cholesky (O(ne³))")
    ax.loglog(ne, sparse, "s-", color=TENSION, label="sparse Gram: same equations, sparse LDLᵀ")
    grid = np.array([ne[0], 150_000.0])
    for t, x0, col, p, name in [
        (dense[-1], ne[-1], COMPRESSION, 3.0, "cubic"),
        (sparse[-1], ne[-1], TENSION, np.polyfit(np.log(ne), np.log(sparse), 1)[0], "fit"),
    ]:
        ax.loglog(grid, t * (grid / x0) ** p, ls=":", color=col, lw=1)
        yend = t * (150_000.0 / x0) ** p
        if yend < 1e3:
            txt = f"150 k edges: {yend:.0f} ms"
        elif yend < 3.6e6:
            txt = f"150 k edges: {yend / 1e3:.0f} s"
        else:
            txt = f"150 k edges: {yend / 3.6e6:.0f} h  (matrix {8 * 150_000.0 ** 2 / 1e9:.0f} GB)"
        ax.annotate(
            txt,
            (150_000, yend),
            textcoords="offset points",
            xytext=(-6, 6),
            ha="right",
            fontsize=8.5,
            color=col,
        )
    for n, t in zip(ne, dense):
        ax.annotate(f"{8 * n * n / 1e6:.0f} MB", (n, t), textcoords="offset points", xytext=(-7, 6), ha="right", fontsize=7.5, color=COMPRESSION)
    ax.axvline(150_000, color="k", lw=0.6, ls="--")
    ax.set_xlabel("edges ne")
    ax.set_ylabel("solve time (ms)")
    ax.set_title("The dense trap: Gram normal equations, dense vs sparse")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="upper left", frameon=False)
    fig.savefig(out / "dense_vs_sparse.png")
    plt.close(fig)


def parse_scale(path: Path):
    """Return {method: [(ne, warm_ms, err_clip/L, final/L)]} plus forward/eval timings."""
    out: dict[str, list] = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) < 5 or not parts[0].isdigit():
            continue
        ne = int(parts[1])
        m = parts[3]
        if m in ("forward_solve", "lbfgsb_eval"):
            out.setdefault(m, []).append((ne, float(parts[4]), np.nan, np.nan))
            continue
        if len(parts) < 10 or parts[4] == "skipped" or parts[4] == "-":
            continue
        warm = float(parts[4])
        err = float(parts[5]) if parts[5] != "-" else np.nan
        fin = float(parts[9]) if parts[9] != "-" else np.nan
        out.setdefault(m, []).append((ne, warm, err, fin))
    return {k: np.array(v) for k, v in out.items()}


def fig_scale(data: Path, out: Path):
    path = data / "scale.txt"
    if not path.exists():
        print("scale.txt missing, skipping")
        return
    d = parse_scale(path)
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.3))
    ax = axes[0]
    for label, (name, color, ls) in METHOD_STYLE.items():
        if label not in d or label == "uniform":
            continue
        r = d[label]
        ax.loglog(r[:, 0], r[:, 1] / 1e3, marker="o", color=color, ls=ls, label=name)
    if "gram_dense" in d:
        r = d["gram_dense"]
        ax.loglog(r[:, 0], r[:, 1] / 1e3, "x--", color=COMPRESSION, label="dense Gram (capped at 6 k edges)")
    if "lbfgsb_eval" in d:
        r = d["lbfgsb_eval"]
        ax.loglog(r[:, 0], r[:, 1], "k^:", lw=1, label="1 000 L-BFGS-B evaluations")
    ax.axvline(150_000, color="k", lw=0.6, ls="--")
    ax.set_xlabel("edges ne")
    ax.set_ylabel("warm-start wall time (s)")
    ax.set_title("Warm-start cost against size (corner-anchored quads)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=8)

    ax = axes[1]
    for label, (name, color, ls) in METHOD_STYLE.items():
        if label not in d:
            continue
        r = d[label]
        ax.loglog(r[:, 0], r[:, 2], marker="o", color=color, ls=ls, label=name)
    ax.set_xlabel("edges ne")
    ax.set_ylabel("clipped warm-start error / L")
    ax.set_title("Quality of the handed-over start (jit2%d, snug box)")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    fig.savefig(out / "scale.png")
    plt.close(fig)


# ───────────────────────── suite summary (from the write-up table) ─────────────────────────


def parse_summary_table(doc: Path):
    text = doc.read_text()
    start = text.index("### 4.2 Summary over the 64 cases")
    block = text[start:]
    rows = []
    for line in block.splitlines()[1:]:
        if line.startswith("|") and not line.startswith("|---") and "method" not in line:
            cells = [c.strip().strip("*") for c in line.strip().strip("|").split("|")]
            rows.append(cells)
        elif rows and not line.startswith("|"):
            break
    return rows


def fig_summary(doc: Path, out: Path):
    rows = parse_summary_table(doc)
    keep = {
        "uniform": "uniform",
        "s1 (Stage 1 only)": "s1",
        "gram_sparse": "gram_sparse",
        "length_ratio": "length_ratio",
        "frozen (1 CWLS step)": "frozen",
        "pipeline": "pipeline",
        "legacy (Clarabel Stage 2)": "legacy",
    }
    names, within2, evals, fails, start_ratio = [], [], [], [], []
    for cells in rows:
        if cells[0] not in keep:
            continue
        label = keep[cells[0]]
        names.append(METHOD_STYLE[label][0])
        fails.append(int(cells[1]))
        within2.append(int(cells[3]))
        evals.append(float(cells[6]))
        start_ratio.append(float(cells[8]))
    colors = [METHOD_STYLE[keep[c[0]]][1] for c in rows if c[0] in keep]
    y = np.arange(len(names))
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.9))
    axes[0].barh(y, within2, color=colors)
    axes[0].set_xlim(0, 64)
    axes[0].set_title("cases within 2 % of the best final error\n(of 64, after 1000 L-BFGS-B iterations)")
    for yi, v, f in zip(y, within2, fails):
        axes[0].text(v + 0.8, yi, f"{v}" + (f"  ({f} fail)" if f else ""), va="center", fontsize=8.5)
    axes[1].barh(y, start_ratio, color=colors)
    axes[1].set_xscale("log")
    axes[1].axvline(1, color="k", lw=0.7)
    axes[1].set_xlim(0.7, max(start_ratio) * 4)
    axes[1].set_title("median clipped start / best final\n(1 = optimal before L-BFGS-B moves)")
    for yi, v in zip(y, start_ratio):
        axes[1].text(v * 1.1, yi, f"{v:g}×", va="center", fontsize=8.5)
    axes[2].barh(y, evals, color=colors)
    axes[2].set_title("median evaluations to reach 5 % of the best")
    for yi, v in zip(y, evals):
        axes[2].text(v + 8, yi, f"{v:g}", va="center", fontsize=8.5)
    axes[2].set_xlim(0, max(evals) * 1.18)
    for ax in axes:
        ax.set_yticks(y)
        ax.invert_yaxis()
        ax.grid(True, axis="x", alpha=0.25)
    axes[0].set_yticklabels(names)
    axes[1].set_yticklabels([])
    axes[2].set_yticklabels([])
    fig.suptitle("Suite summary: 16 nets × {jit2%d, bump10%d} × {loose, snug}", y=1.02)
    fig.tight_layout()
    fig.savefig(out / "suite_summary.png")
    plt.close(fig)


# ───────────────────────── pipeline diagram ─────────────────────────


def fig_pipeline(out: Path):
    fig, ax = plt.subplots(figsize=(13, 3.6))
    ax.set_xlim(0, 13)
    ax.set_ylim(0, 3.6)
    ax.set_axis_off()
    boxes = [
        (0.2, "Stage 1\nforce-residual particular\nmin ‖Eq − p‖ in the box\n(Clarabel, t-form, sparse)", "#f2f2f2"),
        (2.9, "Seed guard\nrace Stage 1 vs scaled\nuniform through Stage 2,\nkeep lower geometric error", "#f2f2f2"),
        (5.6, "Frozen CWLS step\nweights D(x*)⁻¹, Jacobian at\nthe target; bounds by sparse\nactive-set BVLS", "#dce8f7"),
        (8.3, "Gauss–Newton ×2\nre-linearise D at x(q),\nLM safeguard, degenerate\nedges fall back to frozen J", "#dce8f7"),
        (11.0, "clip to box\n→ L-BFGS-B\n(target objective,\nadjoint gradient)", "#f2f2f2"),
    ]
    w, h, y0 = 2.4, 2.1, 0.9
    for x0, text, col in boxes:
        ax.add_patch(FancyBboxPatch((x0, y0), w, h, boxstyle="round,pad=0.05", fc=col, ec="#444444", lw=1))
        ax.text(x0 + w / 2, y0 + h / 2, text, ha="center", va="center", fontsize=9)
    for (x0, _, _), (x1, _, _) in zip(boxes, boxes[1:]):
        ax.add_patch(FancyArrowPatch((x0 + w + 0.02, y0 + h / 2), (x1 - 0.02, y0 + h / 2), arrowstyle="-|>", mutation_scale=14, color="#444444"))
    ax.text(
        6.5,
        0.35,
        "identity behind the weighting:  x(q) − x* = −D(q)⁻¹ r(q),   r(q) = E(x*) q − p  (affine in q);   "
        "all steps are sparse LDLᵀ factorisations of a saddle system — no EᵀE, no dense matrix",
        ha="center",
        va="center",
        fontsize=9.2,
        style="italic",
    )
    ax.text(6.5, 3.35, "Warm-start pipeline (default `InverseFdmOptions`, `Stage2Method::ActiveSet`)", ha="center", fontsize=11)
    fig.savefig(out / "pipeline.png")
    plt.close(fig)


# ───────────────────────── main ─────────────────────────


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default="../../docs/figures")
    ap.add_argument("--doc", default="../../docs/warm_start.md")
    args = ap.parse_args()
    data, out, doc = Path(args.data), Path(args.out), Path(args.doc)
    out.mkdir(parents=True, exist_ok=True)

    fig_pipeline(out)
    if (data / "gallery.json").exists():
        fig_gallery(data, out)
    if (data / "reactions_tiedarch16s.json").exists():
        fig_reactions(data, out)
    if (data / "dense.txt").exists():
        fig_dense(data, out)
    fig_scale(data, out)
    if doc.exists():
        fig_summary(doc, out)

    cases = []
    for path in sorted(data.glob("case_*.json")):
        case = load(path)
        stem = path.stem
        fig_case_geometry(case, out, stem)
        fig_case_convergence(case, out, stem)
        cases.append((stem, case))
        print("rendered", stem)
    if cases:
        fig_convergence_grid(cases, out)
        fig_start_vs_final(cases, out)
    print("figures written to", out.resolve())


if __name__ == "__main__":
    main()
