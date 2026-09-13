"""
Export jax_fdm's bundled examples as portable inverse-FDM benchmark cases.

Usage
-----
    uv run python export_jax_fdm.py            # all cases
    uv run python export_jax_fdm.py truss dome # a subset

For every example this script

1. re-runs the ORIGINAL optimisation of the upstream example (same goals,
   optimizer, bounds, iteration budget and tolerance) and takes the optimised
   force densities q* as the REFERENCE equilibrium;
2. verifies with an independent NumPy FDM solve that (q*, fixed nodes, loads)
   reproduce the reference geometry;
3. runs a SECOND "position fit" with jax_fdm's L-BFGS-B, starting from a
   uniform sign-consistent force density, minimising the sum of squared
   distances between the free nodes and the reference geometry (NodePointGoal
   on every free node, weight 1). This is the external baseline the Rust
   pipeline is compared against, together with a 0.1x / 10x seed sweep.

One JSON file per case is written to ``cases/``.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
from math import radians
from math import sqrt
from time import perf_counter
from typing import Any

import numpy as np

import jax

jax.config.update("jax_enable_x64", True)

import jax_fdm  # noqa: E402
from compas.geometry import Polygon  # noqa: E402
from compas.geometry import Translation  # noqa: E402
from compas.geometry import add_vectors  # noqa: E402
from compas.geometry import cross_vectors  # noqa: E402
from compas.geometry import offset_polygon  # noqa: E402
from compas.geometry import rotate_points  # noqa: E402
from compas.geometry import subtract_vectors  # noqa: E402
from compas.itertools import pairwise  # noqa: E402
from compas.topology import dijkstra_path  # noqa: E402
from jax_fdm.constraints import EdgeForceConstraint  # noqa: E402
from jax_fdm.constraints import EdgeLengthConstraint  # noqa: E402
from jax_fdm.constraints import VertexCurvatureConstraint  # noqa: E402
from jax_fdm.datastructures import FDMesh  # noqa: E402
from jax_fdm.datastructures import FDNetwork  # noqa: E402
from jax_fdm.equilibrium import fdm  # noqa: E402
from jax_fdm.equilibrium.fdm import _fdm  # noqa: E402
from jax_fdm.equilibrium.fdm import datastructure_validate  # noqa: E402
from jax_fdm.equilibrium.fdm import model_from_sparsity  # noqa: E402
from jax_fdm.equilibrium.fdm import structure_from_datastructure  # noqa: E402
from jax_fdm.goals import EdgeDirectionGoal  # noqa: E402
from jax_fdm.goals import EdgeForceGoal  # noqa: E402
from jax_fdm.goals import EdgeLengthGoal  # noqa: E402
from jax_fdm.goals import MeshLoadPathGoal  # noqa: E402
from jax_fdm.goals import MeshPlanarityGoal  # noqa: E402
from jax_fdm.goals import MeshSmoothGoal  # noqa: E402
from jax_fdm.goals import NetworkLoadPathGoal  # noqa: E402
from jax_fdm.goals import NodePlaneGoal  # noqa: E402
from jax_fdm.goals import NodePointGoal  # noqa: E402
from jax_fdm.goals import NodeResidualDirectionGoal  # noqa: E402
from jax_fdm.goals import NodeResidualForceGoal  # noqa: E402
from jax_fdm.goals import NodesColinearGoal  # noqa: E402
from jax_fdm.goals import VertexLineGoal  # noqa: E402
from jax_fdm.goals import VertexPointGoal  # noqa: E402
from jax_fdm.goals import VertexResidualForceGoal  # noqa: E402
from jax_fdm.losses import L2Regularizer  # noqa: E402
from jax_fdm.losses import Loss  # noqa: E402
from jax_fdm.losses import MeanSquaredError  # noqa: E402
from jax_fdm.losses import PredictionError  # noqa: E402
from jax_fdm.losses import RootMeanSquaredError  # noqa: E402
from jax_fdm.losses import SquaredError  # noqa: E402
from jax_fdm.optimization import LBFGSB  # noqa: E402
from jax_fdm.optimization import SLSQP  # noqa: E402
from jax_fdm.parameters import EdgeForceDensityParameter  # noqa: E402
from jax_fdm.parameters import VertexSupportXParameter  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
CASES_DIR = os.path.join(HERE, "cases")
UPSTREAM_DATA = os.path.join(HERE, "upstream_data")

# The upstream commit the examples were read from. ``upstream_data/`` holds
# verbatim copies of the two ``data/json`` files the examples need.
JAX_FDM_REPO = "https://github.com/arpastrana/jax_fdm"
JAX_FDM_COMMIT = "9a17fd379071fb3511ca914f46e8d1e8c286460e"

# Position-fit configuration (the external baseline).
FIT_MAXITER = 1000
FIT_TOL = 1e-10
SEED_FACTORS = (10.0, 0.1)
VERIFY_RTOL = 1e-8


# ==========================================================================
# Small helpers
# ==========================================================================


def source_ref(path: str) -> str:
    return f"{path} (commit {JAX_FDM_COMMIT})"


def flist(a: Any) -> list:
    return np.asarray(a, dtype=np.float64).tolist()


# ==========================================================================
# Optimisation driver: constrained_fdm with instrumentation
# ==========================================================================


class RunResult:
    def __init__(self, **kw: Any) -> None:
        self.__dict__.update(kw)


def run_constrained(
    datastructure,
    optimizer,
    loss,
    parameters=None,
    constraints=None,
    maxiter: int = 100,
    tol: float = 1e-6,
    sparse: bool = True,
) -> RunResult:
    """
    Reproduce ``jax_fdm.equilibrium.constrained_fdm`` step by step so that the
    timings, the SciPy result and the per-iteration loss trace can be recorded.

    The optimisation itself is byte-for-byte the library's: ``optimizer.problem``
    builds (and JIT-warms) the objective and ``optimizer.solve`` runs SciPy. The
    only addition is a callback that stores the iterates; the loss trace is
    evaluated afterwards, outside the timed region.
    """
    datastructure_validate(datastructure)

    if constraints and sparse:
        sparse = False  # as constrained_fdm does

    model = model_from_sparsity(sparse=sparse, tmax=1, eta=1e-6)
    structure = structure_from_datastructure(datastructure, sparse)

    t0 = perf_counter()
    problem = optimizer.problem(
        model,
        structure,
        datastructure,
        loss,
        parameters,
        constraints,
        maxiter,
        tol,
        None,
        True,
    )
    t_setup = perf_counter() - t0  # includes JIT compile + warm-up evaluation

    iterates: list[np.ndarray] = []

    def callback(xk, *args, **kwargs):
        iterates.append(np.array(xk, dtype=np.float64))

    problem.callback = callback

    t0 = perf_counter()
    x_opt = optimizer.solve(problem)
    t_solve = perf_counter() - t0

    fun = problem.fun
    trace = [float(fun(x)[0]) for x in iterates]
    loss_end = float(fun(x_opt)[0])

    params = optimizer.parameters_fdm(x_opt)
    out = _fdm(model, params, structure, datastructure)

    res = optimizer.result
    return RunResult(
        datastructure=out,
        x=np.asarray(x_opt, dtype=np.float64),
        trace=trace,
        loss_start=trace[0],
        loss_end=loss_end,
        nit=int(res.nit),
        nfev=int(res.nfev),
        message=str(res.message),
        success=bool(res.success),
        time_s=t_setup + t_solve,
        time_s_excl_jit=t_solve,
    )


# ==========================================================================
# Extraction of a case from a form-found datastructure
# ==========================================================================


def extract_case(ds) -> dict:
    """
    Read nodes, edges, supports, nodal loads and force densities off a
    form-found FDNetwork / FDMesh, re-indexed to 0-based contiguous integers.
    """
    if isinstance(ds, FDMesh):
        keys = list(ds.vertices())
        coords = ds.vertex_coordinates
        is_fixed = ds.is_vertex_support
        load = ds.vertex_load
    else:
        keys = list(ds.nodes())
        coords = ds.node_coordinates
        is_fixed = ds.is_node_support
        load = ds.node_load

    key2idx = {k: i for i, k in enumerate(keys)}
    nodes = [list(map(float, coords(k))) for k in keys]
    fixed = [key2idx[k] for k in keys if is_fixed(k)]
    fixed_set = set(fixed)
    loads = []
    for k in keys:
        i = key2idx[k]
        loads.append([0.0, 0.0, 0.0] if i in fixed_set else list(map(float, load(k))))

    edges = [[key2idx[u], key2idx[v]] for u, v in ds.edges()]
    q = [float(ds.edge_forcedensity(e)) for e in ds.edges()]

    return {
        "nodes": nodes,
        "edges": edges,
        "fixed": fixed,
        "loads": loads,
        "q_ref": q,
        "key2idx": key2idx,
    }


def fdm_forward(nodes, edges, fixed, loads, q) -> np.ndarray:
    """Independent dense NumPy force-density solve for the free node positions."""
    X = np.asarray(nodes, dtype=np.float64)
    P = np.asarray(loads, dtype=np.float64)
    q = np.asarray(q, dtype=np.float64)
    n = len(X)
    E = np.asarray(edges, dtype=int)
    C = np.zeros((len(E), n))
    C[np.arange(len(E)), E[:, 0]] = -1.0
    C[np.arange(len(E)), E[:, 1]] = 1.0
    free = np.array([i for i in range(n) if i not in set(fixed)], dtype=int)
    fix = np.asarray(fixed, dtype=int)
    D = C.T @ (q[:, None] * C)
    rhs = P[free] - D[np.ix_(free, fix)] @ X[fix]
    X_free = np.linalg.solve(D[np.ix_(free, free)], rhs)
    out = X.copy()
    out[free] = X_free
    return out


def stiffness_free(nodes, edges, fixed, q) -> np.ndarray:
    """The free-free block of the force-density matrix C_f^T Q C_f."""
    n = len(nodes)
    E = np.asarray(edges, dtype=int)
    q = np.asarray(q, dtype=np.float64)
    C = np.zeros((len(E), n))
    C[np.arange(len(E)), E[:, 0]] = -1.0
    C[np.arange(len(E)), E[:, 1]] = 1.0
    free = np.array([i for i in range(n) if i not in set(fixed)], dtype=int)
    D = C.T @ (q[:, None] * C)
    return D[np.ix_(free, free)]


def verify_reference(case: dict) -> dict:
    X = np.asarray(case["nodes"])
    Xs = fdm_forward(case["nodes"], case["edges"], case["fixed"], case["loads"], case["q_ref"])
    free = [i for i in range(len(X)) if i not in set(case["fixed"])]
    err = np.linalg.norm(Xs[free] - X[free]) / max(np.linalg.norm(X[free]), 1e-300)
    print(f"  forward-solve check: relative error {err:.3e}")
    if not err < VERIFY_RTOL:
        raise RuntimeError(f"reference equilibrium not reproduced (rel err {err:.3e})")

    eig = np.linalg.eigvalsh(
        stiffness_free(case["nodes"], case["edges"], case["fixed"], case["q_ref"])
    )
    return {
        "rel_err": float(err),
        "Dff_eig_min": float(eig.min()),
        "Dff_eig_max": float(eig.max()),
        "Dff_cond": float(np.abs(eig).max() / np.abs(eig).min()),
    }


def network_from_case(case: dict, q) -> FDNetwork:
    """Build a plain FDNetwork (nodal loads only) from an exported case."""
    net = FDNetwork.from_nodes_and_edges(case["nodes"], [tuple(e) for e in case["edges"]])
    for i in case["fixed"]:
        net.node_support(i)
    for i, p in enumerate(case["loads"]):
        net.node_load(i, list(p))
    for e, qe in zip(net.edges(), q):
        net.edge_forcedensity(e, float(qe))
    return net


# ==========================================================================
# Position fit (external baseline)
# ==========================================================================


def default_bounds(signs: np.ndarray, scale: float) -> tuple[list, list]:
    """Sign-consistent box [|q|/100, 100|q|] scaled by median|q_ref|."""
    lo, hi = [], []
    for s in signs:
        if s > 0:
            lo.append(scale / 100.0)
            hi.append(scale * 100.0)
        else:
            lo.append(-scale * 100.0)
            hi.append(-scale / 100.0)
    return lo, hi


def position_fit(case: dict, q0: np.ndarray, bounds: dict) -> RunResult:
    net = network_from_case(case, q0)
    target = case["target"]
    free = [i for i in range(len(case["nodes"])) if i not in set(case["fixed"])]

    goals = [NodePointGoal(i, target=target[i], weight=1.0) for i in free]
    loss = Loss(SquaredError(goals, alpha=1.0, name="NodePointGoal"))

    lo, hi = bounds["lo"], bounds["hi"]
    parameters = [
        EdgeForceDensityParameter(e, lo[k], hi[k]) for k, e in enumerate(net.edges())
    ]

    optimizer = LBFGSB()
    res = run_constrained(
        net,
        optimizer,
        loss,
        parameters=parameters,
        maxiter=FIT_MAXITER,
        tol=FIT_TOL,
    )

    # final geometry error, recomputed independently from q_final
    Xs = fdm_forward(case["nodes"], case["edges"], case["fixed"], case["loads"], res.x)
    res.final_err = float(np.linalg.norm(Xs[free] - np.asarray(target)[free]))
    return res


def run_position_fits(case: dict, bounds: dict | None) -> tuple[dict, dict]:
    q_ref = np.asarray(case["q_ref"])
    signs = np.sign(q_ref)
    n_zero = int(np.sum(signs == 0))
    if n_zero:
        # An edge whose optimised q sits exactly on a zero bound carries no force.
        # Its sign is taken from the (one-signed) box of the original example.
        if bounds is not None and all(h is not None and h <= 0.0 for h in bounds["hi"]):
            fill = -1.0
        elif bounds is not None and all(lo is not None and lo >= 0.0 for lo in bounds["lo"]):
            fill = 1.0
        else:
            raise RuntimeError("zero force density at the reference; sign pattern undefined")
        print(f"  {n_zero} edges have q_ref == 0 exactly; sign set to {int(fill)} from the example bounds")
        signs = np.where(signs == 0, fill, signs)
    case["signs"] = [int(s) for s in signs]
    case["n_zero_q_ref"] = n_zero

    scale = float(np.median(np.abs(q_ref)))
    bounds_source = "example"
    if bounds is None:
        lo, hi = default_bounds(signs, scale)
        bounds = {"lo": lo, "hi": hi}
        bounds_source = "default sign-consistent box [|q|/100, 100|q|] x median|q_ref|"

    lo_b = np.array([-np.inf if v is None else v for v in bounds["lo"]])
    hi_b = np.array([np.inf if v is None else v for v in bounds["hi"]])

    def seed(factor: float) -> np.ndarray:
        q0 = signs * scale * factor
        return np.clip(q0, lo_b, hi_b)  # SciPy would project x0 anyway

    print(f"  position fit: q0 = sign * {scale:.6g} (median|q_ref|)")
    main = position_fit(case, seed(1.0), bounds)
    print(
        f"    -> {main.nit} it, loss {main.loss_start:.3e} -> {main.loss_end:.3e}, "
        f"final_err {main.final_err:.3e}, {main.time_s_excl_jit:.2f}s"
    )

    sweep = []
    for f in SEED_FACTORS:
        r = position_fit(case, seed(f), bounds)
        print(
            f"  seed x{f:g}: {r.nit} it, loss {r.loss_start:.3e} -> {r.loss_end:.3e}, "
            f"final_err {r.final_err:.3e}, {r.time_s_excl_jit:.2f}s"
        )
        sweep.append(
            {
                "q0": f"uniform sign * {scale * f:.10g} ({f:g} x median|q_ref|, clipped to bounds)",
                "iterations": r.nit,
                "nfev": r.nfev,
                "time_s": r.time_s,
                "time_s_excl_jit": r.time_s_excl_jit,
                "loss_start": r.loss_start,
                "loss_end": r.loss_end,
                "final_err": r.final_err,
                "message": r.message,
            }
        )

    fit = {
        "optimizer": (
            "jax_fdm.optimization.LBFGSB (scipy.optimize.minimize method='L-BFGS-B', "
            f"tol={FIT_TOL:g} -> ftol=gtol={FIT_TOL:g}, maxiter={FIT_MAXITER}, "
            "maxfun/maxls/maxcor = scipy defaults 15000/20/10, jac from jax.value_and_grad, jit)"
        ),
        "q0": f"uniform sign * {scale:.10g} (median|q_ref|)",
        "bounds_source": bounds_source,
        "time_s": main.time_s,
        "time_s_excl_jit": main.time_s_excl_jit,
        "iterations": main.nit,
        "nfev": main.nfev,
        "loss_start": main.loss_start,
        "loss_end": main.loss_end,
        "loss_trace": main.trace,
        "final_err": main.final_err,
        "q_final": flist(main.x),
        "message": main.message,
        "notes": (
            "NodePointGoal on every free node (weight 1) with SquaredError -> loss = sum of "
            "squared distances to `target`; final_err = sqrt(loss) recomputed with an "
            "independent NumPy FDM solve at q_final. loss_trace[k] is the loss at the k-th "
            "SciPy callback (index 0 = q0). time_s includes the JIT compile/warm-up done by "
            "Optimizer.problem(); time_s_excl_jit is the SciPy solve only."
        ),
        "seed_sweep": sweep,
    }
    return fit, bounds


# ==========================================================================
# Case assembly
# ==========================================================================


def original_run_record(res: RunResult, optimizer: str, notes: str, **extra: Any) -> dict:
    rec = {
        "optimizer": optimizer,
        "time_s": res.time_s,
        "time_s_excl_jit": res.time_s_excl_jit,
        "iterations": res.nit,
        "nfev": res.nfev,
        "loss_start": res.loss_start,
        "loss_end": res.loss_end,
        "message": res.message,
        "notes": notes,
    }
    rec.update(extra)
    return rec


def finish_case(
    name: str,
    ref_ds,
    source: str,
    description: str,
    original_run: dict,
    bounds_from_example,
    target_original: list | None = None,
    caveats: list[str] | None = None,
) -> dict:
    """Extract, verify, run the position fits and write the JSON."""
    case = extract_case(ref_ds)
    case.pop("key2idx")
    n_free = len(case["nodes"]) - len(case["fixed"])
    print(
        f"  {len(case['nodes'])} nodes ({len(case['fixed'])} fixed, {n_free} free), "
        f"{len(case['edges'])} edges"
    )
    check = verify_reference(case)
    case["target"] = [list(p) for p in case["nodes"]]

    bounds = None
    if bounds_from_example is not None:
        lo, hi = bounds_from_example
        n_e = len(case["edges"])
        bounds = {
            "lo": [lo] * n_e if not isinstance(lo, list) else lo,
            "hi": [hi] * n_e if not isinstance(hi, list) else hi,
        }

    fit, bounds = run_position_fits(case, bounds)

    out = {
        "name": name,
        "source": "jax_fdm",
        "source_ref": source_ref(source),
        "description": description,
        "nodes": case["nodes"],
        "edges": case["edges"],
        "fixed": case["fixed"],
        "loads": case["loads"],
        "q_ref": case["q_ref"],
        "signs": case["signs"],
        "target": case["target"],
        "bounds": bounds,
        "external": {
            "tool": f"jax_fdm {jax_fdm.__version__} (jax {jax.__version__}, numpy {np.__version__})",
            "reference_check_rel_err": check["rel_err"],
            "reference_Dff": {k: v for k, v in check.items() if k != "rel_err"},
            "original_run": original_run,
            "position_fit_run": fit,
        },
    }
    if target_original is not None:
        out["target_original"] = target_original

    n_t = sum(1 for s in case["signs"] if s > 0)
    n_c = len(case["signs"]) - n_t
    out["sign_mix"] = {"tension": n_t, "compression": n_c}
    out["n_zero_q_ref"] = case["n_zero_q_ref"]

    caveats = list(caveats or [])
    if not np.any(np.asarray(case["loads"])):
        caveats.append(
            "All nodal loads are zero (pure prestress): the geometry is invariant to a "
            "uniform scaling of q, so the exact-target inverse problem has a one-parameter "
            "family of solutions and the seed sweep starts from identical losses."
        )
    if case["n_zero_q_ref"]:
        caveats.append(
            f"{case['n_zero_q_ref']} edges have q_ref == 0 (on the zero bound of the original "
            "example); their sign was set from the example's one-signed box."
        )
    if check["Dff_eig_min"] < 0.0 < check["Dff_eig_max"]:
        caveats.append(
            "Mixed sign: the free-free force-density matrix is indefinite at the reference "
            f"(eigenvalues {check['Dff_eig_min']:.3g} .. {check['Dff_eig_max']:.3g}); the "
            "forward problem is still uniquely solvable but the position-fit landscape is "
            "non-convex."
        )
    runs = [("median seed", fit)] + [
        (f"seed {f:g}x", s) for f, s in zip(SEED_FACTORS, fit["seed_sweep"])
    ]
    for label, r in runs:
        if r["final_err"] > 1e-2:
            caveats.append(
                f"{label}: L-BFGS-B did not reach the exact target (final_err "
                f"{r['final_err']:.3g} after {r['iterations']} iterations, '{r['message']}')."
            )
    out["caveats"] = caveats

    path = os.path.join(CASES_DIR, f"{name}.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=None, separators=(",", ":"))
    print(f"  wrote {os.path.relpath(path, HERE)}")
    return out


# ==========================================================================
# Adapters (one per upstream example)
# ==========================================================================


def case_truss_equal_force() -> dict:
    """examples/truss/truss_equal_force.py — MIXED SIGN."""
    length, num_segments, depth, target_force = 5.0, 11, -2.0, 1.5
    step = length / num_segments
    xs = [-length / 2.0 + i * step for i in range(num_segments + 1)]
    top_points = [[x, 0.0, 0.0] for x in xs]
    bottom_points = [[x, depth, 0.0] for x in xs[1:-1]]

    network = FDNetwork()
    top_nodes = [network.add_node(x=x, y=y, z=z) for x, y, z in top_points]
    top_nodes_free = top_nodes[1:-1]
    bottom_nodes = [network.add_node(x=x, y=y, z=z) for x, y, z in bottom_points]
    bottom_chord = [top_nodes[0]] + bottom_nodes + [top_nodes[-1]]
    top_edges = [network.add_edge(u, v) for u, v in pairwise(top_nodes)]
    bottom_edges = [network.add_edge(u, v) for u, v in pairwise(bottom_chord)]
    strut_edges = [network.add_edge(u, v) for u, v in zip(bottom_nodes, top_nodes_free)]

    network.node_support(top_nodes[0])
    network.node_support(top_nodes[-1])
    for node in top_nodes_free:
        network.node_load(node, [0.0, 0.0, -0.2])
    for edge in top_edges:
        network.edge_forcedensity(edge, -1.0)
    for edge in bottom_edges:
        network.edge_forcedensity(edge, 2.0)
    for edge in strut_edges:
        network.edge_forcedensity(edge, -0.5)

    goals_force = [EdgeForceGoal(edge, target=target_force) for edge in bottom_edges]
    goals_colinear = [NodesColinearGoal(key=top_nodes)]
    loss = Loss(
        MeanSquaredError(goals_force, name="ForceTarget"),
        PredictionError(goals_colinear, name="ColinearTopChord"),
    )
    res = run_constrained(network, LBFGSB(), loss, maxiter=10000, tol=1e-9)

    return finish_case(
        "jaxfdm_truss_equal_force",
        res.datastructure,
        "examples/truss/truss_equal_force.py",
        "Planar 'truss' network in the xy plane: a straight top chord in compression, a "
        "hanging bottom chord in tension and vertical struts in compression, loaded on "
        "the free top nodes. The original example makes every bottom-chord force equal "
        "to 1.5 while keeping the top chord colinear.",
        original_run_record(
            res,
            "LBFGSB (maxiter=10000, tol=1e-9, unbounded force densities)",
            "MeanSquaredError(EdgeForceGoal bottom chord = 1.5) + "
            "PredictionError(NodesColinearGoal top chord); parameters: all edge force "
            "densities, no bounds.",
        ),
        bounds_from_example=None,
    )


def case_arches(vertical_comp: float = 0.8) -> dict:
    """examples/arch/arches.py — one of the six reaction-direction variants."""
    length_arch, num_segments, q_init, pz = 5.0, 10, -1, -0.3
    length_segment = length_arch / num_segments
    xs = [-length_arch / 2.0 + i * length_segment for i in range(num_segments + 1)]
    nodes = [[x, 0.0, 0.0] for x in xs]
    edges = [(i, i + 1) for i in range(num_segments)]
    network = FDNetwork.from_nodes_and_edges(nodes, edges)
    network.node_support(0)
    network.node_support(num_segments)
    for node in network.nodes_free():
        network.node_load(node, load=[0.0, 0.0, pz])
    for edge in network.edges():
        network.edge_forcedensity(edge, q_init)
    parameters = [EdgeForceDensityParameter(edge) for edge in network.edges()]

    goals = [
        NodeResidualDirectionGoal(0, target=[-1.0, 0.0, -vertical_comp]),
        NodeResidualDirectionGoal(num_segments, target=[1.0, 0.0, -vertical_comp]),
    ]
    res = run_constrained(
        network, SLSQP(), Loss(SquaredError(goals)), parameters=parameters, maxiter=200, tol=1e-9
    )
    tag = f"{vertical_comp:g}".replace(".", "p")
    return finish_case(
        f"jaxfdm_arch_reaction_{tag}",
        res.datastructure,
        "examples/arch/arches.py",
        "A 10-segment planar arch (11 nodes, 2 supports) under uniform vertical nodal "
        f"loads. The original example steers the support reaction directions to "
        f"[+-1, 0, -{vertical_comp:g}] (one of six variants swept in the script).",
        original_run_record(
            res,
            "SLSQP (maxiter=200, tol=1e-9, unbounded force densities)",
            "SquaredError(NodeResidualDirectionGoal at both supports, target direction "
            f"[+-1,0,-{vertical_comp:g}]); parameters: all edge force densities, no bounds.",
        ),
        bounds_from_example=None,
        caveats=[
            "The optimum of the original problem is a uniform q (a parabolic arch), so the "
            "median-seeded position fit starts at the exact solution (0 iterations); only the "
            "seed sweep is informative."
        ],
    )


def case_arch_loadpath() -> dict:
    """examples/arch/arch_loadpath.py — load-path minimisation with length constraints."""
    arch_length, num_segments = 5.0, 10
    segment_length = arch_length / num_segments
    xs = [-arch_length / 2.0 + i * segment_length for i in range(num_segments + 1)]
    nodes = [[x, 0.0, 0.0] for x in xs]
    edges = [(i, i + 1) for i in range(num_segments)]
    network = FDNetwork.from_nodes_and_edges(nodes, edges)
    network.edges_forcedensities(q=-1.0)
    network.nodes_anchors(keys=[n for n in network.nodes() if network.is_leaf(n)])
    network.nodes_loads([0.0, 0.0, -0.3])

    loss = Loss(PredictionError(goals=[NetworkLoadPathGoal()]))
    constraints = [EdgeLengthConstraint(edge, 0.75, 1.0) for edge in network.edges()]
    res = run_constrained(network, SLSQP(), loss, constraints=constraints)  # defaults: 100 / 1e-6

    return finish_case(
        "jaxfdm_arch_loadpath",
        res.datastructure,
        "examples/arch/arch_loadpath.py",
        "A 10-segment planar arch under uniform vertical nodal loads whose total load "
        "path is minimised subject to hard edge-length constraints 0.75 <= L <= 1.0.",
        original_run_record(
            res,
            "SLSQP (maxiter=100, tol=1e-6 defaults, unbounded force densities, hard constraints)",
            "PredictionError(NetworkLoadPathGoal) with EdgeLengthConstraint(0.75, 1.0) on all "
            "edges (dense solve because constraints are present). The upstream script also "
            "loads the anchors with [0,0,-0.3]; those loads only enter the reactions and are "
            "stored as zeros here.",
        ),
        bounds_from_example=None,
    )


def case_pringle() -> dict:
    """examples/pringle/pringle.py."""
    length_vault, width_vault, num_u, num_v = 6.0, 3.0, 10, 9
    q_init, pz, rz_min, rz_max = -0.25, -0.1, 0.45, 2.0

    network = FDNetwork()
    xyz_origin = [0.0, 0.0, 0.0]
    length_u = length_vault / (num_u - 1)
    length_v = width_vault / (num_v - 1)
    arches, cross_edges = [], []
    for i in range(num_v):
        arch = []
        start = add_vectors(xyz_origin, [0.0, i * length_v, 0.0])
        for j in range(num_u):
            x, y, z = add_vectors(start, [j * length_u, 0.0, 0.0])
            arch.append(network.add_node(x=x, y=y, z=z))
        arches.append(arch)
        a = i * num_u
        b = a + num_u - 1
        for u, v in zip(range(a, b), range(a + 1, b + 1)):
            network.add_edge(u, v)
    for i in range(num_u):
        seq = [arch[i] for arch in arches]
        for u, v in zip(seq[:-1], seq[1:]):
            cross_edges.append(network.add_edge(u, v))

    for arch in arches:
        network.node_anchor(arch[0])
        network.node_anchor(arch[-1])
    for node in network.nodes_free():
        network.node_load(node, load=[0.0, 0.0, pz])
    for edge in network.edges():
        network.edge_forcedensity(edge, q_init)
    network.transform(Translation.from_vector([-length_vault / 2.0, -width_vault / 2.0, 0.0]))

    num_steps = (num_v - 1) / 2.0
    step_size = (rz_max - rz_min) / num_steps
    rzs = [rz_min + i * step_size for i in range(int(num_steps) + 1)]
    rzs = rzs + rzs[0:-1][::-1]

    parameters = [EdgeForceDensityParameter(edge, -5.0, -0.1) for edge in network.edges()]

    goals_a = []
    for rz, arch in zip(rzs, arches):
        goals_a.append(NodeResidualForceGoal(arch[0], target=rz, weight=100.0))
        goals_a.append(NodeResidualForceGoal(arch[-1], target=rz, weight=100.0))
    goals_b = []
    for node in network.nodes_free():
        origin = network.node_coordinates(node)
        goals_b.append(NodePlaneGoal(node, target=(origin, [1.0, 0.0, 0.0]), weight=10.0))
    goals_c = [EdgeLengthGoal(e, target=network.edge_length(e), weight=1.0) for e in cross_edges]

    loss = Loss(
        SquaredError(goals_a, alpha=1.0, name="ResidualForceGoal"),
        SquaredError(goals_b, alpha=1.0, name="NodePlaneGoal"),
        SquaredError(goals_c, alpha=1.0, name="EdgeLengthGoal"),
    )
    res = run_constrained(network, SLSQP(), loss, parameters=parameters, maxiter=1000, tol=1e-9)

    return finish_case(
        "jaxfdm_pringle",
        res.datastructure,
        "examples/pringle/pringle.py",
        "A 10 x 9 compression grid (a 'pringle' vault, 6 m x 3 m) anchored along the two "
        "short ends. The original example targets a graded distribution of support "
        "reaction magnitudes, keeps free nodes on their transverse planes and preserves the "
        "transverse edge lengths.",
        original_run_record(
            res,
            "SLSQP (maxiter=1000, tol=1e-9, bounds [-5, -0.1])",
            "SquaredError(NodeResidualForceGoal at 18 anchors, w=100) + "
            "SquaredError(NodePlaneGoal x-planes on free nodes, w=10) + "
            "SquaredError(EdgeLengthGoal on transverse edges, w=1).",
        ),
        bounds_from_example=(-5.0, -0.1),
    )


def case_dome() -> dict:
    """examples/dome/dome.py."""
    diameter, num_sides, num_rings, offset_distance = 1.0, 8, 30, 0.01
    q0_ring, q0_cross, pz = -2.0, -0.5, -0.1
    qmin, qmax = None, -1e-6
    maxiter, tol = 10000, 1e-8
    length_target = 0.03
    angle_vector, angle_base, angle_top = [0.0, 0.0, 1.0], 15.0, 45.0

    network = FDNetwork()
    polygon = Polygon.from_sides_and_radius_xy(num_sides, diameter / 2.0).points
    rings = []
    for _ in range(num_rings + 1):
        polygon = offset_polygon(polygon, offset_distance)
        rings.append([network.add_node(x=x, y=y, z=z) for x, y, z in polygon])

    edges_rings = []
    for ring in rings[1:]:
        for u, v in pairwise(ring + ring[:1]):
            edges_rings.append(network.add_edge(u, v))
    for i in range(num_sides):
        radial = [ring[i] for ring in rings]
        for u, v in pairwise(radial):
            network.add_edge(u, v)
    edges_cross_rings = [list(zip(*rings_pair)) for rings_pair in pairwise(rings)]

    for node in rings[0]:
        network.node_anchor(node)
    for node in network.nodes_free():
        network.node_load(node, load=[0.0, 0.0, pz])
    q0_scale = sqrt(network.number_of_edges()) / 2.0
    for edge in edges_rings:
        network.edge_forcedensity(edge, q0_ring)
    for i, cross_ring in enumerate(edges_cross_rings):
        for edge in cross_ring:
            network.edge_forcedensity(edge, q0_cross * q0_scale * (num_rings - i))

    parameters = [EdgeForceDensityParameter(edge, qmin, qmax) for edge in network.edges()]

    goals = []
    for cross_ring in edges_cross_rings:
        for edge in cross_ring:
            goals.append(EdgeLengthGoal(edge, target=length_target, weight=1.0))
    for i, cross_ring in enumerate(edges_cross_rings):
        angle = angle_base + (angle_top - angle_base) * (i / (num_rings - 1))
        for u, v in cross_ring:
            xyz = network.node_coordinates(u)
            normal = cross_vectors(network.edge_vector((u, v)), angle_vector)
            point = add_vectors(xyz, angle_vector)
            end = rotate_points([point], -radians(angle), axis=normal, origin=xyz).pop()
            vector = subtract_vectors(end, xyz)
            goals.append(EdgeDirectionGoal((u, v), target=vector, weight=1.0))

    loss = Loss(SquaredError(goals=goals))
    network = fdm(network)  # the script form-finds once before optimising
    res = run_constrained(network, LBFGSB(), loss, parameters=parameters, maxiter=maxiter, tol=tol)

    return finish_case(
        "jaxfdm_dome",
        res.datastructure,
        "examples/dome/dome.py",
        "An 8-sided, 30-ring compression dome network (radial + hoop edges) anchored at "
        "the outer ring under uniform vertical nodal loads. The original example gives the "
        "radial edges a target length of 0.03 and target inclinations graded from 15 to 45 "
        "degrees.",
        original_run_record(
            res,
            "LBFGSB (maxiter=10000, tol=1e-8, bounds (-inf, -1e-6])",
            "SquaredError(EdgeLengthGoal 0.03 on radial edges + EdgeDirectionGoal graded "
            "15..45 deg on radial edges), started from the fdm() solution of the initial q.",
        ),
        bounds_from_example=(None, -1e-6),
    )


def case_creased_shell() -> dict:
    """examples/creased_shell/creased_shell.py — shape approximation with an explicit target."""
    network = FDNetwork.from_json(os.path.join(UPSTREAM_DATA, "creased_shell.json"))
    targets = {node: network.node_coordinates(node) for node in network.nodes()}

    supports = [node for node in network.nodes() if network.is_leaf(node)]
    network.nodes_supports(supports)
    network.nodes_loads([0.0, 0.0, -0.2], keys=network.nodes_free())
    network.edges_forcedensities(q=-1.0)

    parameters = [EdgeForceDensityParameter(edge, -20.0, 0.0) for edge in network.edges()]
    goals = [NodePointGoal(node, target=targets[node]) for node in network.nodes_free()]
    loss = Loss(RootMeanSquaredError(goals))
    res = run_constrained(network, LBFGSB(), loss, parameters=parameters, maxiter=1000, tol=1e-6)

    keys = list(res.datastructure.nodes())
    target_original = [list(map(float, targets[k])) for k in keys]

    return finish_case(
        "jaxfdm_creased_shell",
        res.datastructure,
        "examples/creased_shell/creased_shell.py (data/json/creased_shell.json; the same "
        "structure appeared earlier as examples/butt, examples/vault and examples/creased_vault)",
        "A creased compression shell network (193 nodes, 37 leaf supports) under uniform "
        "vertical loads. The original example is a best-fit problem: it drives every free "
        "node onto the designer's target surface stored in the input JSON (kept here as "
        "target_original).",
        original_run_record(
            res,
            "LBFGSB (maxiter=1000, tol=1e-6, bounds [-20, 0])",
            "RootMeanSquaredError(NodePointGoal on every free node -> target_original).",
        ),
        bounds_from_example=(-20.0, 0.0),
        target_original=target_original,
    )


def case_cablenet() -> dict:
    """examples/cablenet/cablenet.py — tension only."""
    length, nx, support_height = 10.0, 10, 5.0
    force_boundary, force_interior, length_interior = 20.0, 1.0, 1.0
    qmin, qmax, maxiter, tol = 0.1, None, 1000, 1e-8

    mesh = FDMesh.from_meshgrid(length, nx=nx)
    edges_boundary = [e for e in mesh.edges() if mesh.is_edge_on_boundary(e)]
    edges_interior = [e for e in mesh.edges() if not mesh.is_edge_on_boundary(e)]

    corners = list(mesh.vertices_where(vertex_degree=2))
    mesh.vertices_supports(corners)
    mesh.vertex_attribute(corners[0], "z", support_height)
    mesh.vertex_attribute(corners[-1], "z", support_height)
    for edge in mesh.edges():
        force = force_boundary if mesh.is_edge_on_boundary(edge) else force_interior
        mesh.edge_forcedensity(edge, force / mesh.edge_length(edge))

    parameters = [EdgeForceDensityParameter(edge, qmin, qmax) for edge in mesh.edges()]
    goals_force = [EdgeForceGoal(edge, target=force_boundary) for edge in edges_boundary]
    goals_length = [EdgeLengthGoal(edge, target=length_interior) for edge in edges_interior]
    loss = Loss(
        MeanSquaredError(goals_force, name="BoundaryForce"),
        MeanSquaredError(goals_length, name="InteriorLength"),
    )
    res = run_constrained(mesh, LBFGSB(), loss, parameters=parameters, maxiter=maxiter, tol=tol)

    return finish_case(
        "jaxfdm_cablenet",
        res.datastructure,
        "examples/cablenet/cablenet.py",
        "A 10 x 10 square cable-net (tension only, no loads) with the four corners anchored "
        "and one diagonal pair of corners lifted by 5 m. The original example prestresses "
        "the net so that the boundary cables carry 20 kN and the interior cables are 1 m long.",
        original_run_record(
            res,
            "LBFGSB (maxiter=1000, tol=1e-8, bounds [0.1, inf))",
            "MeanSquaredError(EdgeForceGoal 20 on boundary edges) + "
            "MeanSquaredError(EdgeLengthGoal 1.0 on interior edges). Pure prestress: all nodal "
            "loads are zero.",
        ),
        bounds_from_example=(0.1, None),
    )


def _monkey_saddle_mesh_and_supports():
    mesh = FDMesh.from_json(os.path.join(UPSTREAM_DATA, "monkey_saddle.json"))
    mesh = mesh.subdivided(scheme="quad", k=2)

    corners = set(v for v in mesh.vertices() if mesh.vertex_degree(v) == 2)
    boundary = mesh.vertices_on_boundaries()[0]
    if boundary[0] == boundary[-1]:
        boundary = boundary[:-1]
    start = next(i for i, v in enumerate(boundary) if v in corners)
    boundary = boundary[start:] + boundary[:start]

    polyedges, polyedge = [], [boundary[0]]
    for v in boundary[1:] + [boundary[0]]:
        polyedge.append(v)
        if v in corners:
            polyedges.append(polyedge)
            polyedge = [v]
    polyedge2length = {
        tuple(pe): sum(mesh.edge_length((u, v)) for u, v in pairwise(pe)) for pe in polyedges
    }
    supports = []
    mean_length = sum(polyedge2length.values()) / len(polyedge2length)
    for pe, length in polyedge2length.items():
        if length < mean_length:
            supports += pe
    supports = set(supports)

    steps = {}
    adjacency = mesh.adjacency
    weight = {(u, v): 1.0 for u in adjacency for v in adjacency[u]}
    for v in supports:
        if v in corners:
            steps[v] = 0
        else:
            steps[v] = min(
                len(dijkstra_path(adjacency, weight, v, c)) - 1 for c in corners
            )
    max_step = max(steps.values())
    steps = {v: max_step - s for v, s in steps.items()}
    return mesh, supports, steps, max_step


def case_monkey_saddle() -> dict:
    """examples/monkey_saddle/monkey_saddle.py."""
    q0, pz, qmin, qmax = -2.0, -1.0, -20.0, -0.01
    rmin, rmax, r_exp = 4.0, 8.0, 0.5
    weight_length, weight_residual, alpha, alpha_lp = 1.0, 10.0, 0.1, 0.01
    maxiter, tol = 500, 1e-3

    mesh, supports, steps, max_step = _monkey_saddle_mesh_and_supports()
    for v in supports:
        mesh.vertex_support(v)
    mesh.vertices_loads([0.0, 0.0, pz], keys=list(mesh.vertices_free()))
    mesh.edges_forcedensities(q0)

    parameters = [EdgeForceDensityParameter(edge, qmin, qmax) for edge in mesh.edges()]
    goals_a = [
        EdgeLengthGoal(e, mesh.edge_length(e), weight=weight_length) for e in mesh.edges()
    ]
    goals_b = []
    for key in mesh.vertices_supports():
        step = steps[key]
        reaction = (1 - step / max_step) ** r_exp * (rmax - rmin) + rmin
        goals_b.append(VertexResidualForceGoal(key, reaction, weight=weight_residual))
    goals_c = [MeshLoadPathGoal()]

    loss = Loss(
        SquaredError(goals_a, alpha=1.0, name="EdgeLengthGoal"),
        SquaredError(goals_b, alpha=1.0, name="ReactionForceGoal"),
        PredictionError(goals_c, alpha=alpha_lp, name="LoadPathGoal"),
        L2Regularizer(alpha=alpha),
    )
    mesh = fdm(mesh)  # the script form-finds once before optimising
    res = run_constrained(mesh, LBFGSB(), loss, parameters=parameters, maxiter=maxiter, tol=tol)

    return finish_case(
        "jaxfdm_monkey_saddle",
        res.datastructure,
        "examples/monkey_saddle/monkey_saddle.py (data/json/monkey_saddle.json, quad-subdivided k=2)",
        "A monkey-saddle compression shell mesh (coarse mesh quad-subdivided twice) anchored "
        "along its three shorter boundary polyedges under uniform vertical loads. The "
        "original example keeps the flat-mesh edge lengths, grades the support reactions "
        "between 4 and 8 and penalises load path with an L2 regulariser.",
        original_run_record(
            res,
            "LBFGSB (maxiter=500, tol=1e-3, bounds [-20, -0.01])",
            "SquaredError(EdgeLengthGoal = initial lengths, w=1) + "
            "SquaredError(VertexResidualForceGoal graded 4..8 at supports, w=10) + "
            "0.01 * PredictionError(MeshLoadPathGoal) + L2Regularizer(alpha=0.1).",
        ),
        bounds_from_example=(-20.0, -0.01),
    )


def case_pillow() -> dict:
    """examples/pillow/pillow.py — SLSQP with hard constraints; random q0 seeded."""
    random.seed(0)
    l1, l2, divisions = 10.0, 10.0, 8
    q0, dq, pz = -2.0, 0.1, -100.0
    maxiter = 1000  # the upstream script defines tol=1e-3 but never passes it
    ratio_length_min, ratio_length_max = 0.5, 3.0
    force_min, force_max = -100.0, -1.0
    crv_min, crv_max = -100.0, -0.1

    mesh = FDMesh.from_meshgrid(dx=l1, nx=divisions, dy=l2, ny=divisions)
    for key in mesh.vertices():
        if mesh.is_vertex_on_boundary(key):
            mesh.vertex_support(key)
    for edge in mesh.edges():
        mesh.edge_forcedensity(edge, q0 + dq * (random.random() - 0.5))

    mesh_area = mesh.area()
    for key in mesh.vertices():
        mesh.vertex_load(key, load=[0.0, 0.0, pz * mesh.vertex_area(key) / mesh_area])

    goals = []
    for vertex in mesh.vertices_free():
        xyz = mesh.vertex_coordinates(vertex)
        line = [xyz, add_vectors(xyz, [0.0, 0.0, 1.0])]
        goals.append(VertexLineGoal(vertex, target=line, weight=1.0))
    loss = Loss(SquaredError(goals=goals))

    constraints = []
    average_length = np.mean([mesh.edge_length(e) for e in mesh.edges()])
    for edge in mesh.edges():
        constraints.append(
            EdgeLengthConstraint(
                edge,
                bound_low=ratio_length_min * average_length,
                bound_up=ratio_length_max * average_length,
            )
        )
    for edge in mesh.edges():
        constraints.append(EdgeForceConstraint(edge, bound_low=force_min, bound_up=force_max))
    stride = divisions + 1
    mid_column = divisions // 2
    for key in [mid_column * stride + row for row in range(1, divisions)]:
        polygon = mesh.vertex_neighbors(key, ordered=True)
        constraints.append(
            VertexCurvatureConstraint(key, polygon, bound_low=crv_min, bound_up=crv_max)
        )

    # the script runs an unconstrained SLSQP first; the constrained one is the design
    run_constrained(mesh, SLSQP(), loss, maxiter=maxiter)
    res = run_constrained(mesh, SLSQP(), loss, constraints=constraints, maxiter=maxiter)

    return finish_case(
        "jaxfdm_pillow",
        res.datastructure,
        "examples/pillow/pillow.py (random q0 perturbation seeded with random.seed(0))",
        "An 8 x 8 quad 'pillow' compression shell with the whole boundary supported and a "
        "total vertical load of 100 distributed by tributary area (frozen nodal loads). The "
        "original example keeps every free vertex on its vertical line while enforcing hard "
        "edge-length, edge-force and mid-column curvature constraints.",
        original_run_record(
            res,
            "SLSQP (maxiter=1000, tol=1e-6 default, unbounded q, hard constraints)",
            "SquaredError(VertexLineGoal vertical lines through the flat-mesh vertices) with "
            "EdgeLengthConstraint [0.5, 3.0] x mean length, EdgeForceConstraint [-100, -1] and "
            "VertexCurvatureConstraint [-100, -0.1] on the 7 mid-column vertices (dense solve). "
            "Loads on the supported boundary vertices are stored as zeros.",
        ),
        bounds_from_example=None,
    )


def case_gridshell_planarization() -> dict:
    """examples/gridshell/gridshell_planarization.py — includes sliding supports."""
    length, nx = 10.0, 8
    q0, q0_boundary, qmin, qmax = -1.0, -5.0, -50.0, -0.01
    planarity_weight, smooth_weight, shape_weight = 1.0, 0.03, 0.07

    mesh = FDMesh.from_meshgrid(length, nx=nx)
    mesh.transform(Translation.from_vector([-length / 2.0, -length / 2.0, 0.0]))
    corners = list(mesh.vertices_where(vertex_degree=2))
    side = list(mesh.vertices_where(x=length / 2.0))
    for vertex in corners + side:
        mesh.vertex_support(vertex)
    for vertex in mesh.vertices_free():
        mesh.vertex_load(vertex, [0.0, 0.0, -1.0])
    for edge in mesh.edges():
        if mesh.is_edge_on_boundary(edge) and not mesh.is_edge_fully_supported(edge):
            mesh.edge_forcedensity(edge, q0_boundary)
        else:
            mesh.edge_forcedensity(edge, q0)

    shell = fdm(mesh)
    shape = {v: shell.vertex_coordinates(v) for v in mesh.vertices_free()}

    parameters = [EdgeForceDensityParameter(edge, qmin, qmax) for edge in mesh.edges()]
    for vertex in side:
        x = shell.vertex_coordinates(vertex)[0]
        xtol = 0.1 * length
        parameters.append(VertexSupportXParameter(vertex, x - xtol, x + xtol))

    loss = Loss(
        PredictionError([MeshPlanarityGoal()], alpha=planarity_weight, name="Planarity"),
        MeanSquaredError(
            [VertexPointGoal(v, target=shape[v]) for v in mesh.vertices_free()],
            alpha=shape_weight,
            name="ShapeFidelity",
        ),
        PredictionError([MeshSmoothGoal()], alpha=smooth_weight, name="Smoothness"),
    )
    res = run_constrained(mesh, LBFGSB(), loss, parameters=parameters, maxiter=5000, tol=1e-8)

    return finish_case(
        "jaxfdm_gridshell_planarization",
        res.datastructure,
        "examples/gridshell/gridshell_planarization.py",
        "An 8 x 8 quad compression gridshell (10 m square) pinned at the four corners and "
        "along one full side, under unit vertical nodal loads. The original example "
        "planarises the quad faces while holding the funicular shape and smoothing, and "
        "also lets the pinned-side supports slide in x (the exported fixed coordinates are "
        "the optimised ones).",
        original_run_record(
            res,
            "LBFGSB (maxiter=5000, tol=1e-8, q bounds [-50, -0.01], support-x bounds +-1 m)",
            "PredictionError(MeshPlanarityGoal, a=1) + 0.07 * MeanSquaredError(VertexPointGoal "
            "= fdm() shape on free vertices) + 0.03 * PredictionError(MeshSmoothGoal); "
            "parameters: all edge force densities + VertexSupportXParameter on the pinned side.",
        ),
        bounds_from_example=(-50.0, -0.01),
    )


# ==========================================================================
# Summary table
# ==========================================================================

SKIPPED_EXAMPLES = [
    ("examples/animation/*", "visualisation-only scripts (viewer animations), no optimisation"),
    ("examples/arch/arch.py, arch_plotter.py", "forward fdm() only / plotting only; no inverse problem"),
    (
        "examples/monkey_saddle/monkey_saddle_constraints.py",
        "same mesh as monkey_saddle with TrustRegionConstrained + hard length constraints "
        "(tol 1e-2); dropped as a duplicate topology",
    ),
    (
        "examples/pringle/pringle_sequential.py, pringle_temporal*.py, dome/dome2.py, "
        "dome_constraints.py, dome_sequential.py, dome_temporal.py (historical)",
        "variants of the exported pringle/dome problems, removed upstream",
    ),
    (
        "examples/butt, examples/vault, examples/creased_vault (historical)",
        "the same 193-node / 324-edge structure as creased_shell under earlier names",
    ),
    (
        "tensegrity tower, rhon-klinikum cablenet, hexagonal/butt-hinged examples",
        "not present in any commit of arpastrana/jax_fdm; the only mixed-sign example "
        "upstream is truss/truss_equal_force.py",
    ),
]


def write_summary() -> str:
    paths = sorted(
        p for p in os.listdir(CASES_DIR) if p.startswith("jaxfdm_") and p.endswith(".json")
    )
    cases = [json.load(open(os.path.join(CASES_DIR, p))) for p in paths]

    def fmt(x: float) -> str:
        return f"{x:.3g}"

    lines = [
        "# jax_fdm baselines",
        "",
        f"Tool: `{cases[0]['external']['tool']}`, upstream commit `{JAX_FDM_COMMIT}` "
        f"({JAX_FDM_REPO}).",
        "",
        "Position fit: `jax_fdm.optimization.LBFGSB` (SciPy L-BFGS-B, `tol=1e-10` -> "
        "`ftol=gtol=1e-10`, `maxiter=1000`, SciPy defaults `maxfun=15000`, `maxls=20`, "
        "`maxcor=10`), `NodePointGoal` on every free node with weight 1 in a `SquaredError` "
        "loss, q0 = sign * median|q_ref| (seed sweep: 10x and 0.1x, clipped to the box). "
        "`final_err` = sqrt(sum of squared free-node distances) recomputed with an independent "
        "NumPy FDM solve at `q_final`. Times are the SciPy solve only (JIT compile excluded; "
        "the JSON files also carry the inclusive time).",
        "",
        "| case | n_nodes | n_edges | n_fixed | sign mix (T/C) | original goals | original opt | "
        "original time [s] / iters | fit time [s] / iters / final_err | seed 10x (iters / final_err) | "
        "seed 0.1x (iters / final_err) | caveats |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for c in cases:
        ext = c["external"]
        o = ext["original_run"]
        f = ext["position_fit_run"]
        s10, s01 = f["seed_sweep"]
        cav = list(c.get("caveats", []))
        for lab, s in (("median seed", f), ("seed 10x", s10), ("seed 0.1x", s01)):
            # runs with a large final_err are already flagged in the case caveats
            if s["final_err"] <= 1e-2 and "ITERATIONS REACHED LIMIT" in s["message"]:
                cav.append(
                    f"{lab}: hit the {FIT_MAXITER}-iteration cap (final_err {fmt(s['final_err'])})."
                )
        lines.append(
            "| {name} | {nn} | {ne} | {nf} | {t}/{cmp} | {goals} | {opt} | {ot} / {oi} | "
            "{ft} / {fi} / {fe} | {s10i} / {s10e} | {s01i} / {s01e} | {cav} |".format(
                name=c["name"],
                nn=len(c["nodes"]),
                ne=len(c["edges"]),
                nf=len(c["fixed"]),
                t=c["sign_mix"]["tension"],
                cmp=c["sign_mix"]["compression"],
                goals=o["notes"].replace("|", "/"),
                opt=o["optimizer"].replace("|", "/"),
                ot=fmt(o["time_s_excl_jit"]),
                oi=o["iterations"],
                ft=fmt(f["time_s_excl_jit"]),
                fi=f["iterations"],
                fe=fmt(f["final_err"]),
                s10i=s10["iterations"],
                s10e=fmt(s10["final_err"]),
                s01i=s01["iterations"],
                s01e=fmt(s01["final_err"]),
                cav=" ".join(cav).replace("|", "/") or "-",
            )
        )

    lines += ["", "## Mixed-sign cases", ""]
    mixed = [c["name"] for c in cases if c["sign_mix"]["tension"] and c["sign_mix"]["compression"]]
    lines += [f"- `{n}`" for n in mixed] or ["- none"]

    lines += ["", "## Examples not exported", ""]
    lines += [f"- {what}: {why}" for what, why in SKIPPED_EXAMPLES]
    lines += [
        "",
        "## Conventions",
        "",
        "- q > 0 tension, q < 0 compression (jax_fdm convention); 0-based node indices in "
        "`list(datastructure.nodes())` / `vertices()` order, edges in `edges()` order.",
        "- `nodes` is the reference equilibrium (the upstream example's own optimised result); "
        "`target == nodes` for every case; `target_original` is present only where the upstream "
        "example fits an explicit designer target.",
        "- `loads` are the frozen nodal loads actually used by the example; loads that the "
        "upstream script assigned to supported nodes are stored as zeros (they only affect "
        "reactions).",
        "- `bounds` are the upstream example's box where it had one (`null` = unbounded side), "
        "otherwise the sign-consistent default box; `bounds_source` in `position_fit_run` says "
        "which.",
        "- `external.reference_check_rel_err` is the relative error of an independent NumPy "
        "FDM solve reproducing `nodes` from (`q_ref`, fixed nodes, `loads`).",
        "",
    ]
    text = "\n".join(lines)
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    path = os.path.join(HERE, "results", "jax_fdm_summary.md")
    with open(path, "w") as fh:
        fh.write(text)
    print(f"\nwrote {os.path.relpath(path, HERE)}")
    return text


CASES = {
    "truss": case_truss_equal_force,
    "arch_reaction": case_arches,
    "arch_loadpath": case_arch_loadpath,
    "pringle": case_pringle,
    "dome": case_dome,
    "creased_shell": case_creased_shell,
    "cablenet": case_cablenet,
    "monkey_saddle": case_monkey_saddle,
    "pillow": case_pillow,
    "gridshell": case_gridshell_planarization,
}


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("cases", nargs="*", help=f"subset of: {', '.join(CASES)}")
    parser.add_argument(
        "--summary-only",
        action="store_true",
        help="only rebuild results/jax_fdm_summary.md from the existing case files",
    )
    args = parser.parse_args(argv)

    if args.summary_only:
        write_summary()
        return 0

    names = args.cases or list(CASES)
    unknown = [n for n in names if n not in CASES]
    if unknown:
        parser.error(f"unknown case(s): {unknown}")

    os.makedirs(CASES_DIR, exist_ok=True)
    print(f"jax_fdm {jax_fdm.__version__}, jax {jax.__version__}, x64={jax.config.jax_enable_x64}")

    failures = {}
    for name in names:
        print(f"\n=== {name} ===")
        t0 = perf_counter()
        try:
            CASES[name]()
        except Exception as exc:  # noqa: BLE001
            failures[name] = repr(exc)
            print(f"  FAILED: {exc!r}")
        print(f"  ({perf_counter() - t0:.1f}s)")

    write_summary()

    if failures:
        print("\nFailed cases:")
        for k, v in failures.items():
            print(f"  {k}: {v}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
