"""Export compas_cem (CEM) case studies as FDM benchmark cases.

Reproduces the published compas_cem examples (repo ``examples/`` and the
CAD 2023 paper's numerical-validation notebooks from
https://github.com/arpastrana/cem_ad_cad), runs the CEM form-finding and,
where the example defines one, the constrained optimisation, and writes one
JSON case per example to ``cases/cem_<name>.json`` plus a summary table to
``results/compas_cem_summary.md``.

Setup (run from ``bench/external/``)::

    uv venv .venv-cem --python 3.12
    uv pip install --python .venv-cem/bin/python --no-deps /tmp/compas_cem_src
    uv pip install --python .venv-cem/bin/python "compas>=2.15,<3" "numpy>=1.26" autograd "nlopt>2.7"

``/tmp/compas_cem_src`` is a clone of https://github.com/arpastrana/compas_cem
(the PyPI wheel 0.8.6 pins compas==1.17.10, which does not import on Python
3.12). The paper data is expected in a clone of ``cem_ad_cad`` (see
``CEM_AD_CAD`` below; override with the environment variable of the same
name). Viewer/plotter extras are not needed and are never imported.

Run (``--no-project`` keeps uv from picking up the jax_fdm ``pyproject.toml``
that lives in the same folder)::

    uv run --no-project --python .venv-cem/bin/python export_compas_cem.py            # all cases
    uv run --no-project --python .venv-cem/bin/python export_compas_cem.py tree_2d bridge_2d
    uv run --no-project --python .venv-cem/bin/python export_compas_cem.py --list

Conventions of the exported JSON (see ``build_case``):

* node indices are 0-based and follow the CEM node keys sorted ascending
  (the original keys are kept in ``external.cem.node_keys``);
* ``q_ref`` = axial force / geometric length at the reference equilibrium,
  tension positive, compression negative (CEM's own sign convention);
* CEM support nodes (including auto-generated auxiliary-trail supports) are
  the FDM ``fixed`` nodes; every other node is free;
* the reference equilibrium is re-solved with a tight CEM convergence
  threshold (``eta``) so that it is a genuine FDM equilibrium to ~1e-12.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import subprocess
import sys
import time
import warnings
from dataclasses import dataclass, field
from math import copysign, cos, floor, pi, sin, sqrt
from typing import Callable

import numpy as np

import compas_cem
from compas.geometry import Line, Polyline, Translation, offset_polyline
from compas.geometry import closest_point_on_line, closest_point_on_plane
from compas.itertools import pairwise
from compas_cem.diagrams import TopologyDiagram
from compas_cem.elements import DeviationEdge, Node, TrailEdge
from compas_cem.equilibrium import static_equilibrium
from compas_cem.loads import NodeLoad
from compas_cem.optimization import (
    DeviationEdgeParameter,
    LineGoal,
    Optimizer,
    OriginNodeYParameter,
    OriginNodeZParameter,
    PlaneGoal,
    PointGoal,
    TrailEdgeForceGoal,
    TrailEdgeParameter,
)
from compas_cem.supports import NodeSupport

HERE = os.path.dirname(os.path.abspath(__file__))
CASES_DIR = os.path.join(HERE, "cases")
RESULTS_DIR = os.path.join(HERE, "results")
COMPAS_CEM_SRC = os.environ.get("COMPAS_CEM_SRC", "/tmp/compas_cem_src")
CEM_AD_CAD = os.environ.get("CEM_AD_CAD", "/tmp/cem_ad_cad")

# Convergence thresholds tried (in order) for the exported reference equilibrium.
TIGHT_ETAS = (1e-13, 1e-12, 1e-11, 1e-10)
TIGHT_TMAX = 10000


# ------------------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------------------


def git_sha(path: str) -> str:
    try:
        out = subprocess.run(
            ["git", "-C", path, "rev-parse", "--short=12", "HEAD"],
            capture_output=True, text=True, check=True,
        )
        return out.stdout.strip()
    except Exception:  # noqa: BLE001
        return "unknown"


@contextlib.contextmanager
def quiet():
    """Silence compas_cem's print chatter and its internal 0/0 warnings.

    compas_cem normalises the trail vector with ``rvec / |rvec|`` and then
    checks for NaN, so a zero-force trail always emits a NumPy warning.
    """
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        yield


def tight_equilibrium(topology: TopologyDiagram):
    """Re-solve the CEM equilibrium with the tightest ``eta`` that converges."""
    last_exc = None
    for eta in TIGHT_ETAS:
        try:
            with quiet():
                form = static_equilibrium(topology, eta=eta, tmax=TIGHT_TMAX)
            return form, eta
        except ValueError as exc:  # "Over N iters. Residual ... > eta"
            last_exc = exc
    raise RuntimeError(f"CEM did not converge for any eta in {TIGHT_ETAS}: {last_exc}")


def timed_form_finding(topology: TopologyDiagram, **kwargs):
    """One CEM forward form-finding with the example's own settings, timed."""
    t0 = time.perf_counter()
    with quiet():
        form = static_equilibrium(topology, **kwargs)
    return form, time.perf_counter() - t0


def load_compas1_network(path: str):
    """Read a COMPAS 1 ``Network`` JSON (as shipped with the paper data)."""
    with open(path) as fh:
        data = json.load(fh)
    nodes = {int(k): (v["x"], v["y"], v["z"]) for k, v in data["node"].items()}
    edges = [(int(u), int(v)) for u, nbrs in data["edge"].items() for v in nbrs]
    return nodes, edges


# ------------------------------------------------------------------------------
# Optimisation wrapper
# ------------------------------------------------------------------------------


@dataclass
class OptRun:
    algorithm: str
    time_s: float
    iterations: int
    loss_start: float
    loss_end: float
    n_parameters: int
    n_goals: int
    status: str
    gradient_norm: float
    eps: float
    kappa: float
    iters_max: int
    grad: str
    notes: str = ""

    def to_json(self) -> dict:
        return {
            "algorithm": f"nlopt LD_{self.algorithm}",
            "time_s": self.time_s,
            "iterations": self.iterations,
            "loss_start": self.loss_start,
            "loss_end": self.loss_end,
            "n_parameters": self.n_parameters,
            "n_goals": self.n_goals,
            "status": self.status,
            "gradient_norm_end": None if self.gradient_norm != self.gradient_norm else self.gradient_norm,
            "stopval_eps": self.eps,
            "ftol_abs_kappa": self.kappa,
            "max_evals": self.iters_max,
            "gradient": "autograd (AD)" if self.grad == "AD" else "finite differences",
            "notes": self.notes,
        }


def run_optimizer(opt: Optimizer, topology: TopologyDiagram, algorithm: str, *,
                  iters: int, eps: float = 1e-6, kappa: float = 1e-8,
                  grad: str = "AD", tmax: int = 100, eta: float = 1e-6,
                  notes: str = "") -> OptRun:
    """Run ``Optimizer.solve`` and pin the topology to the returned optimum.

    ``Optimizer.solve`` mutates ``topology`` in place with whatever parameter
    vector nlopt evaluated *last*, which need not be the optimum. We therefore
    re-apply ``x_opt`` explicitly before the reference equilibrium is computed.
    """
    x0 = opt.optimization_parameters(topology)
    loss_start = float(opt._optimize_form(x0, topology.copy(), tmax, eta))

    with quiet():
        opt.solve(topology, algorithm=algorithm, grad=grad, iters=iters, eps=eps,
                  kappa=kappa, tmax=tmax, eta=eta, verbose=False)

    if opt.x_opt is None:
        raise RuntimeError("compas_cem optimisation failed (no x_opt)")
    opt._update_parameters(topology, opt.x_opt)

    return OptRun(
        algorithm=algorithm, time_s=float(opt.time_opt), iterations=int(opt.evals),
        loss_start=loss_start, loss_end=float(opt.penalty),
        n_parameters=opt.number_of_parameters(), n_goals=opt.number_of_goals(),
        status=str(opt.status), gradient_norm=float(opt.gradient_norm),
        eps=eps, kappa=kappa, iters_max=iters, grad=grad, notes=notes,
    )


def goal_targets(opt: Optimizer, form) -> dict[int, list[float]]:
    """Target point per node for the geometric goals (Point/Line/Plane)."""
    targets = {}
    for goal in opt.goals.values():
        node = goal.key()
        if isinstance(goal, PointGoal):
            targets[node] = [float(c) for c in goal.target()]
        elif isinstance(goal, LineGoal):
            targets[node] = [float(c) for c in closest_point_on_line(form.node_coordinates(node), goal._target)]
        elif isinstance(goal, PlaneGoal):
            targets[node] = [float(c) for c in closest_point_on_plane(form.node_coordinates(node), goal._target)]
    return targets


def describe_goals(opt: Optimizer | None) -> dict[str, int]:
    if opt is None:
        return {}
    counts: dict[str, int] = {}
    for goal in opt.goals.values():
        counts[type(goal).__name__] = counts.get(type(goal).__name__, 0) + 1
    return counts


def describe_parameters(opt: Optimizer | None) -> dict[str, int]:
    if opt is None:
        return {}
    counts: dict[str, int] = {}
    for prm in opt.parameters.values():
        counts[type(prm).__name__] = counts.get(type(prm).__name__, 0) + 1
    return counts


# ------------------------------------------------------------------------------
# Export + FDM validation
# ------------------------------------------------------------------------------


@dataclass
class CaseResult:
    name: str
    source_ref: str
    description: str
    topology: TopologyDiagram
    form: object                     # FormDiagram at the reference equilibrium
    eta_ref: float
    form_finding_time_s: float
    opt: Optimizer | None = None
    opt_run: OptRun | None = None
    paper: bool = False
    extra: dict = field(default_factory=dict)


def fdm_check(nodes, edges, fixed, loads, q):
    """Assemble C^T Q C x = p and check the CEM state is an FDM equilibrium."""
    X = np.asarray(nodes, float)
    P = np.asarray(loads, float)
    q = np.asarray(q, float)
    n, m = len(X), len(edges)
    C = np.zeros((m, n))
    for k, (i, j) in enumerate(edges):
        C[k, i] = 1.0
        C[k, j] = -1.0
    free = np.array(sorted(set(range(n)) - set(fixed)), int)
    fixed = np.array(sorted(fixed), int)
    D = C.T @ (q[:, None] * C)
    residual = D @ X - P
    res_free = float(np.abs(residual[free]).max()) if len(free) else 0.0
    res_free_rel = res_free / max(float(np.abs(q[:, None] * (C @ X)).max()), 1e-300)

    Dff = D[np.ix_(free, free)]
    rhs = P[free] - D[np.ix_(free, fixed)] @ X[fixed]
    x_free = np.linalg.solve(Dff, rhs)
    err = float(np.abs(x_free - X[free]).max())
    scale = float(np.abs(X).max())
    cond = float(np.linalg.cond(Dff)) if len(free) else 0.0
    return {
        "residual_free_max_abs": res_free,
        "residual_free_rel": res_free_rel,
        "solve_error_max_abs": err,
        "solve_error_rel": err / scale,
        "cond_Dff": cond,
        "n_free": int(len(free)),
    }


COND_MAX = 1e10


def free_block_cond(X: np.ndarray, edges, q, fixed) -> float:
    n = len(X)
    free = sorted(set(range(n)) - set(fixed))
    if not free:
        return 1.0
    C = np.zeros((len(edges), n))
    for k, (i, j) in enumerate(edges):
        C[k, i], C[k, j] = 1.0, -1.0
    Dff = (C.T @ (np.asarray(q)[:, None] * C))[np.ix_(free, free)]
    return float(np.linalg.cond(Dff))


def grow_fixed_nodes(X: np.ndarray, edges, q, fixed, planar: bool) -> list[int]:
    """Add well-spread fixed nodes until the FDM free block is well conditioned.

    Needed in two situations: (1) a self-stressed CEM result with no supports
    left (``C^T Q C x = 0`` has nullspace containing span{1, x, y(, z)}, so at
    least 3 (2D) / 4 (3D) nodes must be fixed) and (2) a mixed tension/
    compression state whose FDM matrix happens to be singular with the CEM
    supports alone (zero-stiffness mode). Farthest-point sampling picks the
    additional nodes.
    """
    n = len(X)
    chosen = list(fixed)
    minimum = 0 if chosen else (3 if planar else 4)
    while len(chosen) < n:
        if chosen and len(chosen) >= minimum and free_block_cond(X, edges, q, chosen) < COND_MAX:
            break
        if chosen:
            d = np.min([np.linalg.norm(X - X[c], axis=1) for c in chosen], axis=0)
        else:
            d = np.linalg.norm(X - X.mean(axis=0), axis=1)
        chosen.append(int(np.argmax(d)))
    return sorted(chosen)


def build_case(res: CaseResult, compas_cem_sha: str) -> dict:
    topo, form = res.topology, res.form

    # Auxiliary trails whose optimised force is machine zero are dropped: CEM
    # then places the auxiliary support on top of its parent node (zero-length
    # edge, undefined force density) and the edge carries nothing anyway.
    all_edges = list(form.edges())
    all_forces = {e: float(form.edge_force(e)) for e in all_edges}
    fmax = max(abs(f) for f in all_forces.values())
    dropped_edges, dropped_nodes = [], set()
    for origin, edge in topo.auxiliary_trails(keys=True):
        edge = tuple(edge)
        if edge not in all_forces:
            edge = (edge[1], edge[0])
        if abs(all_forces[edge]) <= 1e-9 * fmax:
            dropped_edges.append(edge)
            dropped_nodes.add(edge[1] if edge[0] == origin else edge[0])

    keys = sorted(k for k in form.nodes() if k not in dropped_nodes)
    index = {k: i for i, k in enumerate(keys)}
    X = np.array([form.node_coordinates(k) for k in keys], float)

    edge_keys = [e for e in all_edges if e not in dropped_edges]
    edges = [[index[u], index[v]] for u, v in edge_keys]
    forces = np.array([all_forces[e] for e in edge_keys], float)
    lengths = np.array([np.linalg.norm(X[i] - X[j]) for i, j in edges], float)
    if np.any(lengths < 1e-12):
        raise RuntimeError("zero-length edge in the reference equilibrium")
    q = forces / lengths

    signs = []
    for e, f in zip(edge_keys, forces):
        if f != 0.0:
            signs.append(int(copysign(1, f)))
        else:  # zero force: fall back to the prescribed combinatorial state
            state = topo.edge_attribute(e, "length") if topo.is_trail_edge(e) else topo.edge_attribute(e, "force")
            signs.append(int(copysign(1, state)) if state else 1)

    cem_fixed = sorted(index[k] for k in form.support_nodes() if k in index)
    fixed = cem_fixed
    fixed_note = "CEM support nodes"
    planar = bool(np.ptp(X[:, 2]) < 1e-9)
    if not fixed:
        fixed = grow_fixed_nodes(X, edges, q, [], planar)
        fixed_note = ("self-stressed CEM result with no supports left after dropping the zero-force "
                      f"auxiliary trails; {len(fixed)} well-spread nodes fixed for the FDM problem")
    elif free_block_cond(X, edges, q, fixed) >= COND_MAX:
        fixed = grow_fixed_nodes(X, edges, q, fixed, planar)
        extra = sorted(set(fixed) - set(cem_fixed))
        fixed_note = (f"CEM support nodes {cem_fixed} give a singular FDM matrix at the reference state "
                      f"(zero-stiffness mode); additionally fixed nodes {extra}")
    loads = [[float(c) for c in form.node_load(k)] for k in keys]
    check = fdm_check(X, edges, fixed, loads, q)

    aux_edges = sorted(index_edge for index_edge, e in enumerate(edge_keys) if topo.is_auxiliary_trail_edge(e))
    n_trail = sum(1 for e in edge_keys if topo.is_trail_edge(e))
    n_dev = len(edge_keys) - n_trail
    description = res.description
    if dropped_edges:
        description += (f" Export note: {len(dropped_edges)} auxiliary trail(s) with machine-zero optimised "
                        "force (and their auxiliary support nodes) were dropped.")
    if not fixed_note.startswith("CEM"):
        description += f" Export note: {fixed_note}."

    case = {
        "name": f"cem_{res.name}",
        "source": "compas_cem",
        "source_ref": f"{res.source_ref} (compas_cem commit {compas_cem_sha})",
        "description": description,
        "nodes": X.tolist(),
        "edges": edges,
        "fixed": fixed,
        "loads": loads,
        "q_ref": q.tolist(),
        "signs": signs,
        "target": X.tolist(),
        "bounds": None,
        "external": {
            "tool": f"compas_cem {compas_cem.__version__}",
            "form_finding_time_s": res.form_finding_time_s,
            "optimization_run": res.opt_run.to_json() if res.opt_run else None,
            "cem": {
                "node_keys": keys,
                "n_trail_edges": n_trail,
                "n_deviation_edges": n_dev,
                "n_indirect_deviation_edges": topo.number_of_indirect_deviation_edges(),
                "n_trails": topo.number_of_trails(),
                "n_auxiliary_trails": topo.number_of_auxiliary_trails(),
                "auxiliary_trail_edges": aux_edges,
                "dropped_zero_force_auxiliary_trails": len(dropped_edges),
                "fixed_nodes_note": fixed_note,
                "origin_nodes": sorted(index[k] for k in topo.origin_nodes() if k in index),
                "reference_eta": res.eta_ref,
                "reactions": {str(index[k]): [float(c) for c in form.reaction_force(k)]
                              for k in form.support_nodes() if k in index},
                "goals": describe_goals(res.opt),
                "parameters": describe_parameters(res.opt),
                "from_cad2023_paper": res.paper,
                **res.extra,
            },
            "fdm_check": check,
        },
    }

    if res.opt is not None:
        targets = goal_targets(res.opt, form)
        if targets:
            case["target_original"] = [targets.get(k) for k in keys]

    return case


# ------------------------------------------------------------------------------
# Case adapters (repo examples)
# ------------------------------------------------------------------------------


def case_braced_tower_2d() -> CaseResult:
    """examples/02_braced_tower_2d.py: form-finding only."""
    points = [(0, [0.0, 0.0, 0.0]), (1, [0.0, 1.0, 0.0]), (2, [0.0, 2.0, 0.0]),
              (3, [1.0, 0.0, 0.0]), (4, [1.0, 1.0, 0.0]), (5, [1.0, 2.0, 0.0])]
    topology = TopologyDiagram()
    for key, point in points:
        topology.add_node(Node(key, point))
    for u, v in [(0, 1), (1, 2), (3, 4), (4, 5)]:
        topology.add_edge(TrailEdge(u, v, length=-1.0))
    for u, v in [(1, 4), (2, 5)]:
        topology.add_edge(DeviationEdge(u, v, force=-1.0))
    for u, v in [(1, 5), (1, 3), (2, 4)]:  # indirect deviation edges (bracing)
        topology.add_edge(DeviationEdge(u, v, force=1.0))
    topology.add_support(NodeSupport(0))
    topology.add_support(NodeSupport(3))
    topology.add_load(NodeLoad(2, [0.0, -1.0, 0.0]))
    topology.add_load(NodeLoad(5, [0.0, -1.0, 0.0]))
    topology.build_trails()

    _, t_ff = timed_form_finding(topology, eta=1e-6, tmax=100)
    form, eta = tight_equilibrium(topology)
    return CaseResult(
        name="braced_tower_2d",
        source_ref="examples/02_braced_tower_2d.py",
        description=("Two-storey braced tower in 2D: two compression trails (columns), two "
                     "compression deviation edges (beams) and three tension indirect deviation "
                     "edges (bracing); unit downward loads at the top nodes. Pure CEM "
                     "form-finding, no optimisation."),
        topology=topology, form=form, eta_ref=eta, form_finding_time_s=t_ff,
    )


def case_bridge_2d() -> CaseResult:
    """examples/03_bridge_2d.py: PointGoals on the supports, SLSQP."""
    topology = TopologyDiagram.from_json(os.path.join(COMPAS_CEM_SRC, "examples", "03_bridge_2d.json"))
    topology.add_support(NodeSupport(1))
    topology.add_support(NodeSupport(5))
    topology.add_load(NodeLoad(2, [0.0, -1.0, 0.0]))
    topology.add_load(NodeLoad(6, [0.0, -1.0, 0.0]))
    topology.build_trails()
    _, t_ff = timed_form_finding(topology)

    opt = Optimizer()
    for node, target in zip([1, 5], [(-20.67, 42.7, 0.0), (15.7, 28.84, 0.0)]):
        opt.add_goal(PointGoal(node, target))
    for edge in topology.trail_edges():
        opt.add_parameter(TrailEdgeParameter(edge, bound_low=15.0, bound_up=5.0))
    for edge in topology.deviation_edges():
        opt.add_parameter(DeviationEdgeParameter(edge, bound_low=10.0, bound_up=10.0))
    run = run_optimizer(opt, topology, "SLSQP", iters=100, eps=1e-6)

    form, eta = tight_equilibrium(topology)
    return CaseResult(
        name="bridge_2d",
        source_ref="examples/03_bridge_2d.py + 03_bridge_2d.json",
        description=("Suspended bridge in 2D loaded from the repo JSON: two tension trails (cables) "
                     "rising from the loaded deck nodes 2 and 6 to the supports 1 and 5, cross-linked "
                     "by two tension ties and one compression strut (deviation edges); unit downward "
                     "loads at nodes 2 and 6. Optimisation: PointGoal on the two "
                     "support nodes (1, 5), parameters = all trail-edge lengths (bounds -15/+5) "
                     "and all deviation-edge forces (+-10), nlopt SLSQP, 100 evals, eps 1e-6."),
        topology=topology, form=form, eta_ref=eta, form_finding_time_s=t_ff, opt=opt, opt_run=run,
    )


def case_tree_2d() -> CaseResult:
    """examples/04_tree_2d.py: auxiliary trails, TrailEdgeForceGoal, SLSQP."""
    width = 4.0
    height = width / 2
    topology = TopologyDiagram()
    topology.add_node(Node(1, [-width / 2.0, height, 0.0]))
    topology.add_node(Node(2, [width / 2.0, height, 0.0]))
    topology.add_node(Node(3, [0.0, height / 2.5, 0.0]))
    topology.add_node(Node(4, [0.0, 0.0, 0.0]))
    topology.add_edge(TrailEdge(3, 4, length=-height / 2))
    topology.add_edge(DeviationEdge(1, 3, force=-sqrt(4.0)))
    topology.add_edge(DeviationEdge(2, 3, force=-sqrt(2.0)))
    topology.add_edge(DeviationEdge(1, 2, force=2.0))
    topology.add_support(NodeSupport(4))
    topology.add_load(NodeLoad(1, [0.0, -1.0, 0.0]))
    topology.add_load(NodeLoad(2, [0.0, -1.0, 0.0]))
    with quiet():
        topology.build_trails(auxiliary_trails=True)
    _, t_ff = timed_form_finding(topology, eta=1e-5, tmax=100)

    opt = Optimizer()
    for edge in topology.auxiliary_trail_edges():
        # The upstream example adds this goal twice per edge (effective weight 2); kept as is.
        opt.add_goal(TrailEdgeForceGoal(edge, force=0.0))
        opt.add_goal(TrailEdgeForceGoal(edge, force=0.0))
    opt.add_parameter(DeviationEdgeParameter((1, 2), 1.0, 10.0))
    opt.add_parameter(DeviationEdgeParameter((1, 3), 1.0, 10.0))
    opt.add_parameter(DeviationEdgeParameter((2, 3), 10.0, 1.0))
    run = run_optimizer(opt, topology, "SLSQP", iters=100, eps=1e-6,
                        notes="upstream example registers each TrailEdgeForceGoal twice (weight 2)")

    form, eta = tight_equilibrium(topology)
    return CaseResult(
        name="tree_2d",
        source_ref="examples/04_tree_2d.py",
        description=("Tree-like 2D structure: one compression trunk trail, two compression "
                     "branches and one tension tie as deviation edges, unit downward loads on "
                     "the two canopy nodes; the canopy nodes get automatic auxiliary trails. "
                     "Optimisation: TrailEdgeForceGoal(0) on the auxiliary trails (paper's "
                     "auxiliary-trail trick), parameters = the three deviation-edge forces, "
                     "nlopt SLSQP."),
        topology=topology, form=form, eta_ref=eta, form_finding_time_s=t_ff, opt=opt, opt_run=run,
    )


def tensegrity_wheel_topology(num_sides: int, tension_force: float, compression_force: float,
                              diameter: float = 1.0, appendix_length: float = 0.10) -> TopologyDiagram:
    assert num_sides % 2 == 0
    topology = TopologyDiagram()
    thetas = np.linspace(0.0, 2 * pi, num_sides + 1)[:-1]
    radius = diameter / 2.0
    for i, theta in enumerate(thetas):
        topology.add_node(Node(i, [radius * cos(theta), radius * sin(theta), 0.0]))
    for u, v in pairwise(list(range(num_sides)) + [0]):
        topology.add_edge(DeviationEdge(u, v, force=tension_force))
    half = num_sides // 2
    for u in range(half):
        topology.add_edge(DeviationEdge(u, u + half, force=compression_force))
    topology.auxiliary_trail_length = -appendix_length
    with quiet():
        topology.build_trails(auxiliary_trails=True)
    return topology


def case_tensegrity_wheel_2d() -> CaseResult:
    """examples/05_tensegrity_wheel_2d.py: 16 sides, LBFGS, AD gradients."""
    topology = tensegrity_wheel_topology(16, tension_force=1.0, compression_force=-0.5)
    _, t_ff = timed_form_finding(topology, eta=1e-6, tmax=100)

    opt = Optimizer()
    for edge in topology.auxiliary_trail_edges():
        opt.add_goal(TrailEdgeForceGoal(edge, force=0.0))
    for edge in topology.deviation_edges():
        opt.add_parameter(DeviationEdgeParameter(edge, 2.0, 2.0))
    run = run_optimizer(opt, topology, "LBFGS", iters=100, eps=1e-6, grad="AD")

    form, eta = tight_equilibrium(topology)
    return CaseResult(
        name="tensegrity_wheel_2d",
        source_ref="examples/05_tensegrity_wheel_2d.py",
        description=("Self-stressed planar tensegrity wheel with 16 sides: tension deviation edges "
                     "on the rim (+1), compression spokes through the centre (-0.5), no external "
                     "loads; every rim node gets an automatic auxiliary trail of length 0.1. "
                     "Optimisation: TrailEdgeForceGoal(0) on all 16 auxiliary trails, parameters "
                     "= all 24 deviation-edge forces (+-2), nlopt LBFGS with autograd gradients."),
        topology=topology, form=form, eta_ref=eta, form_finding_time_s=t_ff, opt=opt, opt_run=run,
    )


# ------------------------------------------------------------------------------
# Case adapters (CAD 2023 paper, cem_ad_cad/numerical_validation)
# ------------------------------------------------------------------------------


def case_tensegrity_wheel_paper(num_sides: int = 64) -> CaseResult:
    """cem_ad_cad wheel.ipynb: paper settings (+1/-1, bound 2, LBFGS, 1000 evals)."""
    topology = tensegrity_wheel_topology(num_sides, tension_force=1.0, compression_force=-1.0)
    _, t_ff = timed_form_finding(topology)

    opt = Optimizer()
    for edge in topology.auxiliary_trail_edges():
        opt.add_goal(TrailEdgeForceGoal(edge, force=0.0))
    for edge in topology.deviation_edges():
        opt.add_parameter(DeviationEdgeParameter(edge, bound_low=2.0, bound_up=2.0))
    run = run_optimizer(opt, topology, "LBFGS", iters=1000, eps=1e-6, grad="AD",
                        notes=f"paper sweeps num_sides = 2^2..2^8; this export uses {num_sides}")

    form, eta = tight_equilibrium(topology)
    return CaseResult(
        name=f"tensegrity_wheel_{num_sides}_paper",
        source_ref="cem_ad_cad/numerical_validation/wheel/wheel.ipynb (CAD 2023, Sec. 4.1)",
        description=(f"Self-stressed planar tensegrity wheel with {num_sides} sides as in the CAD 2023 "
                     "paper: rim tension +1, spoke compression -1, no loads, one auxiliary trail "
                     "(length 0.1) per rim node. Optimisation: TrailEdgeForceGoal(0) on all "
                     "auxiliary trails, parameters = all deviation-edge forces (+-2), nlopt LBFGS, "
                     "AD gradients, stopval 1e-6."),
        topology=topology, form=form, eta_ref=eta, form_finding_time_s=t_ff, opt=opt, opt_run=run,
        paper=True,
    )


def case_tree_canopy_3d(algorithm: str = "SLSQP") -> CaseResult:
    """cem_ad_cad tree.ipynb: 3D tree canopy, auxiliary trails everywhere."""
    force_compression, force_tension = -1.0, 0.5
    load = [0.0, 0.0, -0.5]
    base = os.path.join(CEM_AD_CAD, "numerical_validation", "tree")
    network_data = {"edges_compression.json": {"force": force_compression, "load_nodes": False},
                    "edges_tension.json": {"force": force_tension, "load_nodes": True}}

    topology = TopologyDiagram()
    for fname, data in network_data.items():
        nodes, edges = load_compas1_network(os.path.join(base, fname))
        for u, v in edges:
            topology.add_edge(DeviationEdge.from_line((nodes[u], nodes[v]), force=data["force"]))
        if data["load_nodes"]:
            for xyz in nodes.values():
                topology.add_load(NodeLoad.from_point_and_vector(xyz, load))
    with quiet():
        topology.build_trails(auxiliary_trails=True)
    _, t_ff = timed_form_finding(topology)

    opt = Optimizer()
    supports = list(topology.nodes_where({"z": (-0.01, 0.01)}))  # intended (ground) supports
    for edge in topology.auxiliary_trail_edges():
        u, v = edge
        if u in supports or v in supports:
            continue
        opt.add_goal(TrailEdgeForceGoal(edge, force=0.0))
    for edge in topology.deviation_edges():
        opt.add_parameter(DeviationEdgeParameter(edge, bound_up=5.0, bound_low=5.0))
    for node in topology.origin_nodes():
        if node in supports:
            continue
        opt.add_parameter(OriginNodeZParameter(node, bound_up=1.5, bound_low=1.5))
        opt.add_parameter(OriginNodeYParameter(node, bound_up=1.5, bound_low=1.5))
    run = run_optimizer(opt, topology, algorithm, iters=500, eps=1e-6, grad="AD",
                        notes="paper compares LBFGS, SLSQP and AUGLAG on this problem")

    form, eta = tight_equilibrium(topology)
    return CaseResult(
        name="tree_canopy_3d",
        source_ref="cem_ad_cad/numerical_validation/tree/tree.ipynb (CAD 2023, Sec. 4.2)",
        description=("3D tree canopy from the CAD 2023 paper: compression branch/trunk edges (-1) "
                     "and tension canopy edges (+0.5) all modelled as deviation edges; loads of "
                     "-0.5 z on the canopy nodes; every node receives an automatic auxiliary trail "
                     "(default length -1 along (1,1,1)). The auxiliary trails at the ground nodes "
                     "(z=0) act as the real supports. Optimisation: TrailEdgeForceGoal(0) on all "
                     "other auxiliary trails; parameters = all deviation forces (+-5) and the Y/Z "
                     "coordinates of the non-ground origin nodes (+-1.5); nlopt "
                     f"{algorithm}, 500 evals, stopval 1e-6."),
        topology=topology, form=form, eta_ref=eta, form_finding_time_s=t_ff, opt=opt, opt_run=run,
        paper=True, extra={"ground_support_nodes_cem_keys": sorted(supports)},
    )


def case_curved_bridge_3d(num_hangers: int = 10) -> CaseResult:
    """cem_ad_cad bridge.ipynb: curved bridge under torsion, SLSQP."""
    assert num_hangers % 2 == 0
    height = width = 1.0
    forces = {"force_chord_bottom": -5.0, "force_chord_top": 5.0, "force_hanger_vertical": -0.1,
              "force_hanger_deck": -1.5, "force_hanger_tie": 1.5}
    load_line, bound_force, bound_length, target_line_height = -1.0, 10.0, 0.5, 5.0

    with open(os.path.join(CEM_AD_CAD, "numerical_validation", "bridge", "curve.json")) as fh:
        guide_curve = Polyline(json.load(fh)["points"])
    guide_points = guide_curve.points
    target_lines = [Line((x, y, -target_line_height / 2.0), (x, y, target_line_height / 2.0))
                    for x, y, _ in (guide_points[0], guide_points[-1])]

    chord_bottom = Polyline(points=guide_curve.divide(num_hangers - 1))
    chord_top = chord_bottom.transformed(Translation.from_vector([0.0, 0.0, height]))
    chord_deck = Polyline(points=offset_polyline(chord_bottom.points, distance=-width))
    node_load = [0.0, 0.0, load_line * guide_curve.length / num_hangers]

    num_segments = num_hangers - 1
    topology = TopologyDiagram()
    for chord, force in zip((chord_bottom, chord_top), (forces["force_chord_bottom"], forces["force_chord_top"])):
        lines = list(chord.lines)
        line_middle = lines.pop(floor(num_segments / 2.0))  # middle edge is a deviation edge
        topology.add_edge(DeviationEdge.from_line(line_middle, force=force))
        for line in lines:
            topology.add_edge(TrailEdge.from_line(line, length=line.length * copysign(1.0, force)))
    hangers = [((chord_bottom, chord_top), forces["force_hanger_vertical"]),
               ((chord_top, chord_deck), forces["force_hanger_tie"]),
               ((chord_deck, chord_bottom), forces["force_hanger_deck"])]
    for chords, force in hangers:
        for a, b in zip(chords[0].points, chords[1].points):
            topology.add_edge(DeviationEdge.from_line((a, b), force=force))
    for chord in (chord_bottom, chord_top):
        for point in (chord.points[0], chord.points[-1]):
            topology.add_support(NodeSupport.from_point(point))
    for point in chord_deck.points:
        topology.add_load(NodeLoad.from_point_and_vector(point, vector=node_load))
    topology.auxiliary_trail_vector = [0.0, 0.0, 1.0]
    with quiet():
        topology.build_trails(auxiliary_trails=True)
    _, t_ff = timed_form_finding(topology)

    opt = Optimizer()
    for edge in topology.auxiliary_trail_edges():
        opt.add_goal(TrailEdgeForceGoal(edge, force=0.0))
    for chord in (chord_bottom, chord_top):
        for pt, line in zip((chord.points[0], chord.points[-1]), target_lines):
            opt.add_goal(LineGoal(topology.node_key(pt), line))
    for edge in topology.deviation_edges():
        opt.add_parameter(DeviationEdgeParameter(edge, bound_force, bound_force))
    for edge in topology.trail_edges():
        if topology.is_edge_supported(edge) and not topology.is_auxiliary_trail_edge(edge):
            opt.add_parameter(TrailEdgeParameter(edge, bound_length, bound_length))
    run = run_optimizer(opt, topology, "SLSQP", iters=100, eps=1e-6, grad="AD",
                        notes=f"paper sweeps num_hangers = 4..22; this export uses {num_hangers}")

    form, eta = tight_equilibrium(topology)
    return CaseResult(
        name=f"curved_bridge_{num_hangers}_3d",
        source_ref="cem_ad_cad/numerical_validation/bridge/bridge.ipynb + bridge/curve.json (CAD 2023, Sec. 4.3)",
        description=(f"Curved bridge under torsion from the CAD 2023 paper with {num_hangers} hangers: "
                     "a compression bottom chord (-5) and a tension top chord (+5) split by a "
                     "middle deviation edge into two trails each, vertical compression hangers "
                     "(-0.1), tension ties (+1.5) and compression deck struts (-1.5) as deviation "
                     "edges; the deck nodes carry a uniform downward line load and get auxiliary "
                     "trails along +z. Supports at both ends of both chords. Optimisation: "
                     "TrailEdgeForceGoal(0) on the auxiliary trails and LineGoal pulling the four "
                     "chord-end supports onto vertical target lines; parameters = all deviation "
                     "forces (+-10) and the lengths of the trail edges at the supports (+-0.5); "
                     "nlopt SLSQP, 100 evals, stopval 1e-6."),
        topology=topology, form=form, eta_ref=eta, form_finding_time_s=t_ff, opt=opt, opt_run=run,
        paper=True,
    )


# ------------------------------------------------------------------------------
# Registry / driver
# ------------------------------------------------------------------------------

CASES: dict[str, Callable[[], CaseResult]] = {
    "braced_tower_2d": case_braced_tower_2d,
    "bridge_2d": case_bridge_2d,
    "tree_2d": case_tree_2d,
    "tensegrity_wheel_2d": case_tensegrity_wheel_2d,
    "tensegrity_wheel_64_paper": lambda: case_tensegrity_wheel_paper(64),
    "tree_canopy_3d": case_tree_canopy_3d,
    "curved_bridge_10_3d": lambda: case_curved_bridge_3d(10),
    "curved_bridge_22_3d": lambda: case_curved_bridge_3d(22),
}

DROPPED = [
    ("spiral_staircase (CAD 2023 Sec. 5)", "only available as a binary Grasshopper definition "
     "(cem_ad_cad/case_study/staircase.gh, examples/ghpython/spiral_staircase.ghx); no Python source or geometry export"),
    ("examples/ghpython/{bridge_3d,dome,tensegrity_jessen}.ghx", "Grasshopper-only definitions, geometry internalised in Rhino data"),
    ("historical examples/07_arching_tower_3d.py, 06_surface_structure.py", "depend on compas_singular (CoarseQuadMesh densification), not on PyPI and COMPAS 1 only"),
    ("examples/01_quick_start.py", "compression-only three-bar chain (no tension members)"),
]


def summary_row(case: dict) -> dict:
    ext = case["external"]
    cem = ext["cem"]
    opt = ext["optimization_run"]
    q = np.array(case["q_ref"])
    aux = set(cem["auxiliary_trail_edges"])
    n_t = int(sum(1 for i, s in enumerate(case["signs"]) if s > 0))
    n_c = len(q) - n_t
    params = ", ".join(f"{v} {k}" for k, v in cem["parameters"].items()) or "-"
    goals = ", ".join(f"{v} {k}" for k, v in cem["goals"].items()) or "-"
    return {
        "case": case["name"],
        "paper": "yes" if cem["from_cad2023_paper"] else "no",
        "n_nodes": len(case["nodes"]),
        "n_edges": len(case["edges"]),
        "n_fixed": len(case["fixed"]),
        "tension/compression": f"{n_t}/{n_c}" + (f" ({len(aux)} aux)" if aux else ""),
        "cem_edges": f"{cem['n_trail_edges']} trail / {cem['n_deviation_edges']} dev ({cem['n_indirect_deviation_edges']} indirect)",
        "params": f"{opt['n_parameters']} ({params})" if opt else "-",
        "goals": goals,
        "ff_time_s": f"{ext['form_finding_time_s']:.4f}",
        "opt": (f"{opt['algorithm']}, {opt['time_s']:.1f}s, {opt['iterations']} evals, "
                f"loss {opt['loss_start']:.3g} -> {opt['loss_end']:.3g}, {opt['status']}") if opt else "-",
        "fdm_residual": f"{ext['fdm_check']['residual_free_max_abs']:.1e} abs / solve err {ext['fdm_check']['solve_error_rel']:.1e} rel",
    }


def write_summary(cases: list[dict], failures: dict[str, str], compas_cem_sha: str, cem_ad_cad_sha: str) -> str:
    cols = ["case", "paper", "n_nodes", "n_edges", "n_fixed", "tension/compression", "cem_edges",
            "params", "goals", "ff_time_s", "opt", "fdm_residual"]
    lines = [
        "# compas_cem export summary",
        "",
        f"- tool: compas_cem {compas_cem.__version__} (git {compas_cem_sha}), paper data cem_ad_cad {cem_ad_cad_sha}",
        "- `paper` = case reproduced from the CAD 2023 paper's numerical validation notebooks",
        "- `tension/compression` counts edges by the sign of the reference force; `aux` = auxiliary trail edges (force ~ 0 after optimisation)",
        "- `fdm_residual`: max |C^T Q C x - p| on free nodes, and relative error of the FDM re-solve of the free nodes",
        "- `opt`: nlopt algorithm, wall time, number of objective evaluations, loss at start -> end, nlopt status",
        "",
        "| " + " | ".join(cols) + " |",
        "|" + "|".join("---" for _ in cols) + "|",
    ]
    for case in cases:
        row = summary_row(case)
        lines.append("| " + " | ".join(str(row[c]) for c in cols) + " |")

    lines += [
        "",
        "## CEM parameterisation vs. FDM",
        "",
        "CEM parameters are signed trail-edge lengths (or plane offsets), deviation-edge force magnitudes, "
        "and optionally origin-node coordinates / node loads; FDM has one force density per edge.",
        "",
        "| case | n_edges (FDM q) | CEM design params in example | trail lengths | deviation forces | origin coords | converged | stop criterion hit | tolerance |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for case in cases:
        cem = case["external"]["cem"]
        opt = case["external"]["optimization_run"]
        prm = cem["parameters"]
        n_tl = cem["n_trail_edges"]
        n_dv = cem["n_deviation_edges"]
        if opt:
            conv = "yes" if opt["loss_end"] <= opt["stopval_eps"] else f"no (loss {opt['loss_end']:.2g} > eps)"
            tol = f"stopval {opt['stopval_eps']:g}, ftol_abs {opt['ftol_abs_kappa']:g}, max {opt['max_evals']} evals"
            lines.append(f"| {case['name']} | {len(case['edges'])} | {opt['n_parameters']} "
                         f"| {prm.get('TrailEdgeParameter', 0)} of {n_tl} | {prm.get('DeviationEdgeParameter', 0)} of {n_dv} "
                         f"| {sum(v for k, v in prm.items() if k.startswith('OriginNode'))} | {conv} | {opt['status']} | {tol} |")
        else:
            lines.append(f"| {case['name']} | {len(case['edges'])} | (form-finding only: {n_tl} trail lengths + {n_dv} deviation forces prescribed) "
                         f"| {n_tl} of {n_tl} | {n_dv} of {n_dv} | 0 | - | - | eta 1e-6 |")

    lines += ["", "## Dropped / not reproduced", ""]
    for name, why in DROPPED:
        lines.append(f"- {name}: {why}")
    for name, why in failures.items():
        lines.append(f"- {name}: FAILED - {why}")
    lines += [
        "",
        "## Notes",
        "",
        "- Sign convention: compas_cem stores forces positive in tension and negative in compression for both trail "
        "and deviation edges (`Diagram.edge_force`); `q_ref = force / length` needs no sign flip. At support nodes the FDM "
        "residual `C^T Q C x - p` equals compas_cem's `reaction_force` vector.",
        "- The reference geometry is re-solved with `static_equilibrium(eta<=1e-10, tmax=10000)`; with the examples' own "
        "`eta=1e-6` the indirect deviation edges leave nodal residuals of ~1e-7..1e-9.",
        "- `Optimizer.solve` leaves the topology at nlopt's *last evaluated* parameters; the export re-applies `x_opt` before "
        "computing the reference equilibrium.",
        "- Auxiliary trails (paper's extension) are exported as ordinary edges whose far node is fixed; after optimisation their "
        "forces are only as small as the loss threshold allows (`sum f^2 < 1e-6`).",
    ]
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("cases", nargs="*", help="case names (default: all)")
    parser.add_argument("--list", action="store_true", help="list available cases and exit")
    parser.add_argument("--no-summary", action="store_true", help="do not (re)write results/compas_cem_summary.md")
    args = parser.parse_args(argv)

    if args.list:
        print("\n".join(CASES))
        return 0

    names = args.cases or list(CASES)
    unknown = [n for n in names if n not in CASES]
    if unknown:
        parser.error(f"unknown case(s): {unknown}; available: {list(CASES)}")

    os.makedirs(CASES_DIR, exist_ok=True)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    compas_cem_sha = git_sha(COMPAS_CEM_SRC)
    cem_ad_cad_sha = git_sha(CEM_AD_CAD)

    exported: list[dict] = []
    failures: dict[str, str] = {}
    for name in names:
        print(f"=== {name}")
        t0 = time.perf_counter()
        try:
            res = CASES[name]()
            case = build_case(res, compas_cem_sha)
        except Exception as exc:  # noqa: BLE001
            failures[name] = f"{type(exc).__name__}: {exc}"
            print(f"    FAILED: {failures[name]}")
            continue
        path = os.path.join(CASES_DIR, f"{case['name']}.json")
        with open(path, "w") as fh:
            json.dump(case, fh, indent=1)
        chk = case["external"]["fdm_check"]
        opt = case["external"]["optimization_run"]
        print(f"    nodes {len(case['nodes'])}, edges {len(case['edges'])}, fixed {len(case['fixed'])}, "
              f"eta_ref {case['external']['cem']['reference_eta']:g}, ff {case['external']['form_finding_time_s']:.4f}s")
        if opt:
            print(f"    opt {opt['algorithm']}: {opt['time_s']:.1f}s, {opt['iterations']} evals, "
                  f"loss {opt['loss_start']:.4g} -> {opt['loss_end']:.4g}, {opt['status']}")
        print(f"    FDM residual free {chk['residual_free_max_abs']:.2e} (rel {chk['residual_free_rel']:.2e}), "
              f"solve err rel {chk['solve_error_rel']:.2e}, cond {chk['cond_Dff']:.2e}")
        print(f"    wrote {os.path.relpath(path, HERE)} in {time.perf_counter() - t0:.1f}s")
        exported.append(case)

    if not args.no_summary:
        # Summarise every case file present so partial runs do not drop rows.
        all_cases = []
        for fname in sorted(os.listdir(CASES_DIR)):
            if fname.startswith("cem_") and fname.endswith(".json"):
                with open(os.path.join(CASES_DIR, fname)) as fh:
                    all_cases.append(json.load(fh))
        summary = write_summary(all_cases, failures, compas_cem_sha, cem_ad_cad_sha)
        with open(os.path.join(RESULTS_DIR, "compas_cem_summary.md"), "w") as fh:
            fh.write(summary)
        print(summary)

    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
