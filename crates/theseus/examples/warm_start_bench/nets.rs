//! Synthetic cable / strut networks with a known funicular.
//!
//! Every [`Net`] carries a topology, a plan (nominal node positions, only the
//! fixed ones are used as anchors), vertical loads on the free nodes and a
//! "true" signed force-density vector `q_true`. The forward FDM solve at
//! `q_true` is the funicular; targets are built by perturbing it so that the
//! inverse problem has a known feasible answer inside the box.

use ndarray::Array2;
use std::panic::{catch_unwind, AssertUnwindSafe};
use theseus::fdm;
use theseus::sparse::SparseColMatOwned;
use theseus::types::{
    AnchorInfo, Bounds, FdmCache, NetworkTopology, ObjectiveTrait, Problem, SolverOptions,
};

/// Deterministic LCG in `[0, 1)`.
pub fn lcg(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 32) as u32) as f64 / u32::MAX as f64
}

/// Log-uniform sample in `[lo, hi]`.
fn logu(state: &mut u64, lo: f64, hi: f64) -> f64 {
    10f64.powf(lo.log10() + (hi.log10() - lo.log10()) * lcg(state))
}

/// Signed log-uniform force densities spanning `±spread` decades around 1.
fn assign_q(edges: usize, spread: f64, seed: u64, sign: f64) -> Vec<f64> {
    let mut s = seed;
    (0..edges)
        .map(|_| sign * logu(&mut s, 10f64.powf(-spread), 10f64.powf(spread)))
        .collect()
}

pub fn build_incidence(edges: &[(usize, usize)], num_nodes: usize) -> SparseColMatOwned {
    let mut rows = Vec::with_capacity(edges.len() * 2);
    let mut cols = Vec::with_capacity(edges.len() * 2);
    let mut values = Vec::with_capacity(edges.len() * 2);
    for (edge, &(start, end)) in edges.iter().enumerate() {
        rows.extend([edge, edge]);
        cols.extend([start, end]);
        values.extend([-1.0, 1.0]);
    }
    SparseColMatOwned::from_coo(edges.len(), num_nodes, &rows, &cols, &values).unwrap()
}

/// Forward FDM solve; returns free-node positions (n_free × 3). Panics on a
/// singular Laplacian — use [`try_forward`] where that is possible.
pub fn forward(problem: &Problem, q: &[f64]) -> Array2<f64> {
    try_forward(problem, q).expect("forward FDM solve failed")
}

/// Forward FDM solve returning `None` when the Laplacian is singular, the
/// result is non-finite, or the solver panics.
pub fn try_forward(problem: &Problem, q: &[f64]) -> Option<Array2<f64>> {
    if q.len() != problem.topology.num_edges || !q.iter().all(|v| v.is_finite()) {
        return None;
    }
    let out = catch_unwind(AssertUnwindSafe(|| {
        let mut cache = FdmCache::new(problem).ok()?;
        fdm::solve_fdm(&mut cache, q, problem, &Array2::zeros((0, 3)), 0.0).ok()?;
        let x = Array2::from_shape_fn((problem.topology.free_node_indices.len(), 3), |(i, d)| {
            cache.nf[[problem.topology.free_node_indices[i], d]]
        });
        Some(x)
    }))
    .ok()
    .flatten()?;
    if out.iter().all(|v| v.is_finite()) {
        Some(out)
    } else {
        None
    }
}

/// `‖x(q) − x*‖₂` over the free nodes. Panics on a failed forward solve.
pub fn geom_err(problem: &Problem, target: &Array2<f64>, q: &[f64]) -> f64 {
    let x = forward(problem, q);
    (x - target).iter().map(|v| v * v).sum::<f64>().sqrt()
}

/// Geometric error, NaN when the forward solve fails (singular D for
/// out-of-box q) or panics.
pub fn safe_err(problem: &Problem, target: &Array2<f64>, q: &[f64]) -> f64 {
    match try_forward(problem, q) {
        Some(x) => (x - target).iter().map(|v| v * v).sum::<f64>().sqrt(),
        None => f64::NAN,
    }
}

pub struct Net {
    pub name: String,
    pub edges: Vec<(usize, usize)>,
    pub free: Vec<usize>,
    pub fixed: Vec<usize>,
    /// Nominal position of every node; fixed ones are the anchors.
    pub plan: Vec<[f64; 3]>,
    /// Downward load magnitude per free node (negative = upward).
    pub loads: Vec<f64>,
    /// Signed force densities that define the funicular.
    pub q_true: Vec<f64>,
    /// Characteristic length used to normalise errors and depths.
    pub extent: f64,
    /// Edges whose force density the `reactions` subcommand reports
    /// (tie / bottom chord). Empty for most nets.
    pub tie_edges: Vec<usize>,
}

impl Net {
    pub fn n_nodes(&self) -> usize {
        self.plan.len()
    }

    pub fn named(mut self, name: &str) -> Self {
        self.name = name.to_string();
        self
    }

    pub fn fixed_positions(&self) -> Array2<f64> {
        Array2::from_shape_fn((self.fixed.len(), 3), |(i, d)| self.plan[self.fixed[i]][d])
    }

    pub fn problem(
        &self,
        fixed_positions: &Array2<f64>,
        objectives: Vec<Box<dyn ObjectiveTrait>>,
        lower: &[f64],
        upper: &[f64],
        solver: SolverOptions,
    ) -> Problem {
        let incidence = build_incidence(&self.edges, self.n_nodes());
        let free_incidence = incidence.extract_columns(&self.free);
        let fixed_incidence = incidence.extract_columns(&self.fixed);
        let mut loads = Array2::zeros((self.free.len(), 3));
        for i in 0..self.free.len() {
            loads[[i, 2]] = -self.loads[i];
        }
        Problem {
            topology: NetworkTopology {
                incidence,
                free_incidence,
                fixed_incidence,
                num_edges: self.edges.len(),
                num_nodes: self.n_nodes(),
                free_node_indices: self.free.clone(),
                fixed_node_indices: self.fixed.clone(),
            },
            free_node_loads: loads,
            fixed_node_positions: fixed_positions.clone(),
            anchors: AnchorInfo::all_fixed(fixed_positions.clone()),
            objectives,
            bounds: Bounds {
                lower: lower.to_vec(),
                upper: upper.to_vec(),
            },
            solver,
            self_weight: None,
            pressure: None,
        }
    }

    /// Problem with the loose box and no objectives (forward solves, inverse solves).
    pub fn plain_problem(&self) -> Problem {
        let (lo, hi) = self.box_from_true(100.0, 100.0);
        self.problem(
            &self.fixed_positions(),
            Vec::new(),
            &lo,
            &hi,
            SolverOptions::default(),
        )
    }

    /// Per-edge box: |q| in [min|q_group|/f_lo, max|q_group|*f_hi] with the
    /// sign of q_true, where the group is the tension or compression subset.
    /// f < 1 squeezes the box inside the true range (active bounds at the optimum).
    pub fn box_from_true(&self, f_lo: f64, f_hi: f64) -> (Vec<f64>, Vec<f64>) {
        let group_range = |sign: bool| {
            let v: Vec<f64> = self
                .q_true
                .iter()
                .filter(|&&q| (q > 0.0) == sign)
                .map(|q| q.abs())
                .collect();
            let (mn, mx) = (
                v.iter().cloned().fold(f64::INFINITY, f64::min),
                v.iter().cloned().fold(0.0, f64::max),
            );
            if f_lo < 1.0 {
                // "active" box: inset each end by (1 - f_lo) of the group's log-range
                let span = (mx / mn).max(1.0 + 1e-9).ln();
                let inset = (1.0 - f_lo) * span;
                (mn * inset.exp(), mx / inset.exp())
            } else {
                (mn, mx)
            }
        };
        let (pmin, pmax) = group_range(true);
        let (nmin, nmax) = group_range(false);
        let (f_lo, f_hi) = if f_lo < 1.0 { (1.0, 1.0) } else { (f_lo, f_hi) };
        let mut lo = Vec::new();
        let mut hi = Vec::new();
        for &q in &self.q_true {
            if q > 0.0 {
                lo.push(pmin / f_lo);
                hi.push(pmax * f_hi);
            } else {
                lo.push(-nmax * f_hi);
                hi.push(-nmin / f_lo);
            }
        }
        (lo, hi)
    }

    /// Free-node positions of the funicular (forward solve at `q_true`).
    pub fn funicular(&self) -> Array2<f64> {
        forward(&self.plain_problem(), &self.q_true)
    }

    /// Same as [`Net::funicular`] but `None` if the forward solve fails.
    pub fn try_funicular(&self) -> Option<Array2<f64>> {
        try_forward(&self.plain_problem(), &self.q_true)
    }

    /// Vertical extent of the funicular over the free nodes.
    pub fn depth(&self) -> f64 {
        depth_of(&self.funicular())
    }

    /// jitter: white z noise of `jitter`*depth; bump: smooth sin*sin push of `bump`*depth.
    pub fn target(&self, jitter: f64, bump: f64) -> Array2<f64> {
        let mut x = self.funicular();
        let depth = self.depth();
        let (xmin, xmax, ymin, ymax) = x.rows().into_iter().fold(
            (
                f64::INFINITY,
                f64::NEG_INFINITY,
                f64::INFINITY,
                f64::NEG_INFINITY,
            ),
            |a, r| (a.0.min(r[0]), a.1.max(r[0]), a.2.min(r[1]), a.3.max(r[1])),
        );
        let mut state = 0xc0ffee_u64;
        // mirror the perturbation for compression systems so tension/compression targets are exact mirrors
        let mirror = if self.q_true.iter().filter(|&&q| q < 0.0).count() * 2 > self.q_true.len() {
            -1.0
        } else {
            1.0
        };
        for i in 0..x.nrows() {
            let sx = (std::f64::consts::PI * (x[[i, 0]] - xmin) / (xmax - xmin).max(1e-9)).sin();
            let sy = (std::f64::consts::PI * (x[[i, 1]] - ymin) / (ymax - ymin).max(1e-9)).sin();
            x[[i, 2]] +=
                mirror * (jitter * depth * (2.0 * lcg(&mut state) - 1.0) - bump * depth * sx * sy);
        }
        x
    }

    /// Rescale loads so the funicular depth is `frac` of the extent (FDM is linear in p).
    pub fn normalized(mut self, frac: f64) -> Self {
        let depth = self.depth().max(1e-12);
        let s = frac * self.extent / depth;
        for l in &mut self.loads {
            *l *= s;
        }
        self
    }

    /// Exact tension/compression mirror: negate `q_true` and the loads. The
    /// box follows automatically because it is derived from `q_true`.
    pub fn comp(mut self) -> Self {
        for q in &mut self.q_true {
            *q = -*q;
        }
        for l in &mut self.loads {
            *l = -*l;
        }
        self.name = format!("{}_comp", self.name);
        self
    }

    /// Uniform magnitude with the true sign pattern.
    pub fn sign_seed(&self) -> Vec<f64> {
        self.q_true.iter().map(|q| q.signum()).collect()
    }

    /// Full positions (all nodes) from free-node positions and the anchors.
    pub fn full_positions(&self, x_free: &Array2<f64>) -> Array2<f64> {
        let mut xyz = Array2::zeros((self.n_nodes(), 3));
        for (i, &n) in self.fixed.iter().enumerate() {
            for d in 0..3 {
                xyz[[n, d]] = self.plan[self.fixed[i]][d];
            }
        }
        for (i, &n) in self.free.iter().enumerate() {
            for d in 0..3 {
                xyz[[n, d]] = x_free[[i, d]];
            }
        }
        xyz
    }

    /// Edge lengths for the given free-node positions.
    pub fn edge_lengths(&self, x_free: &Array2<f64>) -> Vec<f64> {
        let xyz = self.full_positions(x_free);
        self.edges
            .iter()
            .map(|&(a, b)| {
                (0..3)
                    .map(|d| (xyz[[b, d]] - xyz[[a, d]]).powi(2))
                    .sum::<f64>()
                    .sqrt()
            })
            .collect()
    }
}

pub fn depth_of(x: &Array2<f64>) -> f64 {
    let zmin = x.column(2).iter().cloned().fold(f64::INFINITY, f64::min);
    let zmax = x.column(2).iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    zmax - zmin
}

/// Member-force resultant at each fixed node, `Σ_e q_e (x_other − x_fixed)`
/// (n_fixed × 3). This is the force the members exert on the support, the
/// same convention as `FdmCache::reactions`; the support reaction that
/// balances it is its negative.
pub fn support_reactions(net: &Net, q: &[f64], x_free: &Array2<f64>) -> Array2<f64> {
    let xyz = net.full_positions(x_free);
    let mut slot = vec![usize::MAX; net.n_nodes()];
    for (i, &n) in net.fixed.iter().enumerate() {
        slot[n] = i;
    }
    let mut r = Array2::zeros((net.fixed.len(), 3));
    for (e, &(a, b)) in net.edges.iter().enumerate() {
        for (node, other) in [(a, b), (b, a)] {
            let s = slot[node];
            if s != usize::MAX {
                for d in 0..3 {
                    r[[s, d]] += q[e] * (xyz[[other, d]] - xyz[[node, d]]);
                }
            }
        }
    }
    r
}

/// Column sum of [`support_reactions`].
pub fn net_reaction(r: &Array2<f64>) -> [f64; 3] {
    let mut out = [0.0; 3];
    for row in r.rows() {
        for d in 0..3 {
            out[d] += row[d];
        }
    }
    out
}

// ───────────────────────── grid helper ─────────────────────────

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Anchoring {
    Corners,
    Boundary,
}

struct Grid {
    side: usize,
    edges: Vec<(usize, usize)>,
    free: Vec<usize>,
    fixed: Vec<usize>,
    plan: Vec<[f64; 3]>,
}

fn grid(side: usize, anchoring: Anchoring) -> Grid {
    let node = |r: usize, c: usize| r * side + c;
    let mut edges = Vec::new();
    for r in 0..side {
        for c in 0..side {
            if c + 1 < side {
                edges.push((node(r, c), node(r, c + 1)));
            }
            if r + 1 < side {
                edges.push((node(r, c), node(r + 1, c)));
            }
        }
    }
    let is_fixed = |r: usize, c: usize| match anchoring {
        Anchoring::Corners => (r == 0 || r == side - 1) && (c == 0 || c == side - 1),
        Anchoring::Boundary => r == 0 || r == side - 1 || c == 0 || c == side - 1,
    };
    let mut free = Vec::new();
    let mut fixed = Vec::new();
    let mut plan = Vec::new();
    for r in 0..side {
        for c in 0..side {
            plan.push([c as f64, r as f64, 0.0]);
            if is_fixed(r, c) {
                fixed.push(node(r, c));
            } else {
                free.push(node(r, c));
            }
        }
    }
    Grid {
        side,
        edges,
        free,
        fixed,
        plan,
    }
}

impl Grid {
    fn length(&self) -> f64 {
        (self.side - 1) as f64
    }
    fn into_net(self, name: String, q_true: Vec<f64>) -> Net {
        let extent = self.length();
        Net {
            name,
            loads: vec![1.0; self.free.len()],
            edges: self.edges,
            free: self.free,
            fixed: self.fixed,
            plan: self.plan,
            q_true,
            extent,
            tie_edges: Vec::new(),
        }
    }
}

// ───────────────────────── topologies ─────────────────────────

/// Square quad grid, anchored at the four corners or along the whole boundary;
/// `sign < 0` gives the all-compression version (loads still downward, so the
/// funicular is an arch).
pub fn quad(side: usize, corners_only: bool, sign: f64) -> Net {
    let g = grid(
        side,
        if corners_only {
            Anchoring::Corners
        } else {
            Anchoring::Boundary
        },
    );
    let q = assign_q(g.edges.len(), 0.75, 0x5eed, sign);
    let name = format!(
        "quad{side}{}{}",
        if corners_only { "c" } else { "b" },
        if sign < 0.0 { "_comp" } else { "" }
    );
    g.into_net(name, q)
}

/// Quad grid with a fraction of interior edges removed at random (irregular valence, still quads).
pub fn quad_holes(side: usize, corners_only: bool, frac: f64) -> Net {
    let mut net = quad(side, corners_only, 1.0);
    let mut state = 0xabcdef_u64;
    let mut degree = vec![0usize; net.n_nodes()];
    for &(a, b) in &net.edges {
        degree[a] += 1;
        degree[b] += 1;
    }
    let is_fixed: Vec<bool> = (0..net.n_nodes()).map(|n| net.fixed.contains(&n)).collect();
    let mut keep = Vec::new();
    for &(a, b) in &net.edges {
        let removable = !is_fixed[a] && !is_fixed[b] && degree[a] > 3 && degree[b] > 3;
        if removable && lcg(&mut state) < frac {
            degree[a] -= 1;
            degree[b] -= 1;
        } else {
            keep.push((a, b));
        }
    }
    net.edges = keep;
    net.q_true = assign_q(net.edges.len(), 0.75, 0x5eed, 1.0);
    net.name = format!("holes{side}{}", if corners_only { "c" } else { "b" });
    net
}

/// 45°-rotated lattice clipped to a square; boundary nodes anchored.
pub fn diamond(n: usize) -> Net {
    // nodes at (i, j) with i+j even, 0..=2n
    let m = 2 * n;
    let idx = |i: usize, j: usize| -> Option<usize> {
        if (i + j) % 2 == 0 && i <= m && j <= m {
            Some((i * (m + 1) + j) / 2)
        } else {
            None
        }
    };
    let mut plan = Vec::new();
    let mut fixed = Vec::new();
    let mut free = Vec::new();
    for i in 0..=m {
        for j in 0..=m {
            if let Some(k) = idx(i, j) {
                debug_assert_eq!(k, plan.len());
                plan.push([i as f64, j as f64, 0.0]);
                if i == 0 || j == 0 || i == m || j == m {
                    fixed.push(k);
                } else {
                    free.push(k);
                }
            }
        }
    }
    let mut edges = Vec::new();
    for i in 0..=m {
        for j in 0..=m {
            if let Some(k) = idx(i, j) {
                if let Some(l) = idx(i + 1, j + 1) {
                    edges.push((k, l));
                }
                if j >= 1 {
                    if let Some(l) = idx(i + 1, j - 1) {
                        edges.push((k, l));
                    }
                }
            }
        }
    }
    let ne = edges.len();
    Net {
        name: format!("diamond{n}"),
        edges,
        loads: vec![1.0; free.len()],
        free,
        fixed,
        plan,
        q_true: assign_q(ne, 0.75, 0x5eed, 1.0),
        extent: m as f64,
        tie_edges: Vec::new(),
    }
}

/// Spider web: hub + `rings` concentric polygons with `spokes` nodes each. Outer ring anchored
/// entirely (boundary) or at every 4th node.
pub fn radial(rings: usize, spokes: usize, all_outer_fixed: bool) -> Net {
    let mut plan = vec![[0.0, 0.0, 0.0]];
    let node = |r: usize, s: usize| 1 + (r - 1) * spokes + s;
    for r in 1..=rings {
        for s in 0..spokes {
            let a = 2.0 * std::f64::consts::PI * s as f64 / spokes as f64;
            plan.push([r as f64 * a.cos(), r as f64 * a.sin(), 0.0]);
        }
    }
    let mut edges = Vec::new();
    for s in 0..spokes {
        edges.push((0, node(1, s)));
    }
    for r in 1..=rings {
        for s in 0..spokes {
            edges.push((node(r, s), node(r, (s + 1) % spokes)));
            if r < rings {
                edges.push((node(r, s), node(r + 1, s)));
            }
        }
    }
    let mut fixed = Vec::new();
    let mut free = vec![0];
    for r in 1..=rings {
        for s in 0..spokes {
            if r == rings && (all_outer_fixed || s % 4 == 0) {
                fixed.push(node(r, s));
            } else {
                free.push(node(r, s));
            }
        }
    }
    let ne = edges.len();
    Net {
        name: format!(
            "radial{rings}x{spokes}{}",
            if all_outer_fixed { "b" } else { "c" }
        ),
        edges,
        loads: vec![1.0; free.len()],
        free,
        fixed,
        plan,
        q_true: assign_q(ne, 0.75, 0x5eed, 1.0),
        extent: 2.0 * rings as f64,
        tie_edges: Vec::new(),
    }
}

/// Boundary-anchored quad grid with a square oculus removed from the centre.
pub fn oculus(side: usize, hole: usize) -> Net {
    let g = grid(side, Anchoring::Boundary);
    let lo = (side - hole) / 2;
    let hi = lo + hole;
    let in_hole = |n: usize| {
        let (r, c) = (n / side, n % side);
        r >= lo && r < hi && c >= lo && c < hi
    };
    let edges: Vec<(usize, usize)> = g
        .edges
        .iter()
        .cloned()
        .filter(|&(a, b)| !in_hole(a) && !in_hole(b))
        .collect();
    let free: Vec<usize> = g.free.iter().cloned().filter(|&n| !in_hole(n)).collect();
    let ne = edges.len();
    Net {
        name: format!("oculus{side}h{hole}"),
        edges,
        loads: vec![1.0; free.len()],
        free,
        fixed: g.fixed,
        plan: g.plan,
        q_true: assign_q(ne, 0.75, 0x5eed, 1.0),
        extent: g.length(),
        tie_edges: Vec::new(),
    }
}

/// Boundary-anchored n×n quad net whose `q_true` is 8× larger along the
/// interior grid row `n/2` (a stiff "cable" that creates a ridge in the
/// funicular). Force densities elsewhere are log-uniform over ±0.5 decades.
pub fn crease(n: usize) -> Net {
    let g = grid(n, Anchoring::Boundary);
    let mid = n / 2;
    let mut q = assign_q(g.edges.len(), 0.5, 0x5eed, 1.0);
    for (e, &(a, b)) in g.edges.iter().enumerate() {
        if a / n == mid && b / n == mid {
            q[e] *= 8.0;
        }
    }
    g.into_net(format!("crease{n}"), q)
}

/// Like [`crease`] but the stiff cable follows the staircase path along the
/// main diagonal, (i,i) → (i,i+1) → (i+1,i+1), so it zig-zags across quads.
pub fn crease_diag(n: usize) -> Net {
    let g = grid(n, Anchoring::Boundary);
    let mut q = assign_q(g.edges.len(), 0.5, 0x5eed, 1.0);
    for (e, &(a, b)) in g.edges.iter().enumerate() {
        let (ra, ca) = (a / n, a % n);
        let (rb, cb) = (b / n, b % n);
        let on_stair = (ra == ca && rb == ra && cb == ca + 1) || (cb == rb && ca == cb && rb == ra + 1);
        if on_stair {
            q[e] *= 8.0;
        }
    }
    g.into_net(format!("creasediag{n}"), q)
}

/// Saddle (hypar) net anchored at the four corners only, with corner heights
/// ±0.15·extent. Interior edges are tension cables (q>0); the four boundary
/// chains are compression "edge arches" (q<0) that push outwards against the
/// inward pull of the net. Loads are vertical; the compression magnitude is
/// kept below the level at which D(q) becomes singular. Not normalised: the
/// corner heights dominate the depth.
pub fn hypar_mixed(n: usize) -> Net {
    let mut g = grid(n, Anchoring::Corners);
    let l = g.length();
    let h = 0.15 * l;
    for p in &mut g.plan {
        let (u, v) = (2.0 * p[0] / l - 1.0, 2.0 * p[1] / l - 1.0);
        p[2] = h * u * v;
    }
    let mut state = 0x4a7_u64;
    let on_boundary = |k: usize| {
        let (r, c) = (k / n, k % n);
        r == 0 || r == n - 1 || c == 0 || c == n - 1
    };
    let q: Vec<f64> = g
        .edges
        .iter()
        .map(|&(a, b)| {
            if on_boundary(a) && on_boundary(b) {
                -logu(&mut state, 1.5, 2.5)
            } else {
                logu(&mut state, 0.6, 1.6)
            }
        })
        .collect();
    let mut net = g.into_net(format!("hypar{n}m"), q);
    for l in &mut net.loads {
        *l = 0.02;
    }
    net
}

/// Geiger-style cable dome. The outer compression ring is modelled by fixed
/// nodes at radius `R`, z = 0. Ring `i` (1..=rings, inward) carries an upper
/// ridge node and a lower hoop node per spoke. Ridge cables (tension) run
/// radially between consecutive upper nodes, diagonals (tension) from each
/// lower node to the next-outer upper node, struts (compression, q<0) join the
/// upper and lower node of a ring, and tension hoops link the lower nodes; the
/// innermost upper nodes form a central tension ring. Loads act downward on
/// the upper nodes only; the struts make D(q) indefinite so the ridge rises
/// above the ring while the hoops hang below it. No triangulation between spokes.
pub fn cable_dome(rings: usize, spokes: usize) -> Net {
    let r_out = 10.0;
    let outer = |s: usize| s;
    let upper = |i: usize, s: usize| spokes + 2 * ((i - 1) * spokes + s);
    let lower = |i: usize, s: usize| upper(i, s) + 1;
    let mut plan = Vec::new();
    for s in 0..spokes {
        let a = 2.0 * std::f64::consts::PI * s as f64 / spokes as f64;
        plan.push([r_out * a.cos(), r_out * a.sin(), 0.0]);
    }
    for i in 1..=rings {
        let r = r_out * (1.0 - i as f64 / (rings as f64 + 1.0));
        for s in 0..spokes {
            let a = 2.0 * std::f64::consts::PI * s as f64 / spokes as f64;
            plan.push([r * a.cos(), r * a.sin(), 1.0]);
            plan.push([r * a.cos(), r * a.sin(), -1.0]);
        }
    }
    let mut edges = Vec::new();
    let mut q = Vec::new();
    let mut state = 0xd0e_u64;
    let jitter = |state: &mut u64| logu(state, 0.85, 1.15);
    for i in 1..=rings {
        for s in 0..spokes {
            let up_outer = if i == 1 { outer(s) } else { upper(i - 1, s) };
            edges.push((up_outer, upper(i, s))); // ridge cable
            q.push(1.0 * jitter(&mut state));
            edges.push((lower(i, s), up_outer)); // diagonal cable
            q.push(2.0 * jitter(&mut state));
            edges.push((upper(i, s), lower(i, s))); // strut
            q.push(-1.0 * jitter(&mut state));
            edges.push((lower(i, s), lower(i, (s + 1) % spokes))); // hoop
            q.push(3.0 * jitter(&mut state));
        }
    }
    for s in 0..spokes {
        edges.push((upper(rings, s), upper(rings, (s + 1) % spokes))); // central tension ring
        q.push(3.0 * jitter(&mut state));
    }
    let fixed: Vec<usize> = (0..spokes).collect();
    let mut free = Vec::new();
    let mut loads = Vec::new();
    for i in 1..=rings {
        for s in 0..spokes {
            free.push(upper(i, s));
            loads.push(1.0);
            free.push(lower(i, s));
            loads.push(0.0);
        }
    }
    Net {
        name: format!("cabledome{rings}x{spokes}"),
        edges,
        free,
        fixed,
        plan,
        loads,
        q_true: q,
        extent: 2.0 * r_out,
        tie_edges: Vec::new(),
    }
}

/// Lenticular (Jawerth-type) cable truss in the x–z plane spanning along x.
/// Top chord and bottom chord are tension cables that meet at the two fixed
/// abutments; vertical posts between the chords are in compression (q<0) and
/// the bays are un-triangulated quads. Loads act downward on the top chord.
/// `tie_edges` lists the bottom-chord edges.
pub fn cable_truss(n: usize) -> Net {
    let span = 20.0;
    let mut plan = Vec::new();
    // node 0 = left abutment, node n = right abutment (0..=n along the top),
    // bottom interior nodes n+1 .. 2n-1
    for i in 0..=n {
        let x = span * i as f64 / n as f64;
        plan.push([x, 0.0, 2.0 * (x / span) * (1.0 - x / span) * 4.0]);
    }
    for i in 1..n {
        let x = span * i as f64 / n as f64;
        plan.push([x, 0.0, -2.0 * (x / span) * (1.0 - x / span) * 4.0]);
    }
    let top = |i: usize| i;
    let bottom = |i: usize| if i == 0 || i == n { i } else { n + i };
    let mut edges = Vec::new();
    let mut q = Vec::new();
    let mut tie_edges = Vec::new();
    let mut state = 0x7ru55_u64;
    for i in 0..n {
        edges.push((top(i), top(i + 1))); // top chord: tension
        q.push(logu(&mut state, 3.0, 5.0));
    }
    for i in 0..n {
        tie_edges.push(edges.len());
        edges.push((bottom(i), bottom(i + 1))); // bottom chord: tension
        q.push(logu(&mut state, 6.0, 10.0));
    }
    for i in 1..n {
        edges.push((top(i), bottom(i))); // posts: compression
        q.push(-logu(&mut state, 0.8, 1.2));
    }
    let fixed = vec![0, n];
    let mut free: Vec<usize> = (1..n).collect();
    free.extend(n + 1..2 * n);
    let loads = free.iter().map(|&k| if k <= n { 1.0 } else { 0.0 }).collect();
    Net {
        name: format!("cabletruss{n}"),
        edges,
        free,
        fixed,
        plan,
        loads,
        q_true: q,
        extent: span,
        tie_edges,
    }
}

/// Compression grid shell (all q<0) over an `n × m` rectangular plan, arching
/// across x. The two springing lines x = 0 and x = n−1 are fixed; every third
/// bay a tension tie (q>0) joins the free nodes adjacent to opposite springing
/// lines, so the vault is partially self-tied. Loads act downward and the
/// vault rises. `tie_edges` lists the ties.
pub fn barrel_vault(n: usize, m: usize) -> Net {
    let node = |i: usize, j: usize| j * n + i;
    let mut plan = Vec::new();
    let mut fixed = Vec::new();
    let mut free = Vec::new();
    for j in 0..m {
        for i in 0..n {
            plan.push([i as f64, j as f64, 0.0]);
            if i == 0 || i == n - 1 {
                fixed.push(node(i, j));
            } else {
                free.push(node(i, j));
            }
        }
    }
    let mut edges = Vec::new();
    let mut q = Vec::new();
    let mut state = 0xba77e1_u64;
    for j in 0..m {
        for i in 0..n {
            if i + 1 < n {
                edges.push((node(i, j), node(i + 1, j)));
                q.push(-logu(&mut state, 0.6, 1.6));
            }
            if j + 1 < m {
                edges.push((node(i, j), node(i, j + 1)));
                q.push(-logu(&mut state, 0.6, 1.6));
            }
        }
    }
    let mut tie_edges = Vec::new();
    for j in (1..m).step_by(3) {
        if j + 1 >= m {
            break;
        }
        tie_edges.push(edges.len());
        edges.push((node(1, j), node(n - 2, j)));
        q.push(logu(&mut state, 0.08, 0.12));
    }
    Net {
        name: format!("barrel{n}x{m}"),
        edges,
        loads: vec![1.0; free.len()],
        free,
        fixed,
        plan,
        q_true: q,
        extent: (n - 1) as f64,
        tie_edges,
    }
}

/// Self-tied spoked wheel: outer compression ring (free except every `anchor_every`-th node),
/// inner tension ring, radial tension cables between them. Quads between radials (no diagonals).
pub fn wheel(spokes: usize, anchor_every: usize) -> Net {
    let r_out = 10.0;
    let r_in = 4.0;
    let mut plan = Vec::new();
    for s in 0..spokes {
        let a = 2.0 * std::f64::consts::PI * s as f64 / spokes as f64;
        plan.push([r_out * a.cos(), r_out * a.sin(), 0.0]);
    }
    for s in 0..spokes {
        let a = 2.0 * std::f64::consts::PI * s as f64 / spokes as f64;
        plan.push([r_in * a.cos(), r_in * a.sin(), -1.0]);
    }
    let mut edges = Vec::new();
    let mut q = Vec::new();
    let mut state = 0x77_u64;
    for s in 0..spokes {
        edges.push((s, (s + 1) % spokes)); // outer ring: compression
        q.push(-logu(&mut state, 15.0, 30.0));
    }
    for s in 0..spokes {
        edges.push((spokes + s, spokes + (s + 1) % spokes)); // inner ring: tension
        q.push(logu(&mut state, 3.0, 8.0));
    }
    for s in 0..spokes {
        edges.push((s, spokes + s)); // radial cables: tension
        q.push(logu(&mut state, 0.7, 2.0));
    }
    let mut fixed = Vec::new();
    let mut free = Vec::new();
    for s in 0..spokes {
        if s % anchor_every == 0 {
            fixed.push(s);
        } else {
            free.push(s);
        }
    }
    for s in 0..spokes {
        free.push(spokes + s);
    }
    let loads = free
        .iter()
        .map(|&n| if n < spokes { 0.2 } else { 1.0 })
        .collect();
    Net {
        name: format!("wheel{spokes}a{anchor_every}"),
        edges,
        free,
        fixed,
        plan,
        loads,
        q_true: q,
        extent: 2.0 * r_out,
        tie_edges: Vec::new(),
    }
}

/// Tied arch spanning along x: compression arch polyline over a tension tie,
/// vertical tension hangers; the two abutments are the only anchors and hanger
/// bays are un-triangulated quads. `tie_edges` lists the tie segments.
pub fn tied_arch(bays: usize) -> Net {
    let span = 20.0;
    let rise = 5.0;
    let mut plan = Vec::new();
    // arch nodes 0..=bays, tie nodes bays+1 .. 2*bays-1 (interior only; ends shared with the arch)
    for i in 0..=bays {
        let x = span * i as f64 / bays as f64;
        let z = 4.0 * rise * (x / span) * (1.0 - x / span);
        plan.push([x, 0.0, z]);
    }
    for i in 1..bays {
        plan.push([span * i as f64 / bays as f64, 0.0, 0.0]);
    }
    let tie = |i: usize| if i == 0 || i == bays { i } else { bays + i };
    let mut edges = Vec::new();
    let mut q = Vec::new();
    let mut tie_edges = Vec::new();
    let mut state = 0x99_u64;
    for i in 0..bays {
        edges.push((i, i + 1)); // arch: compression
        q.push(-logu(&mut state, 8.0, 16.0));
    }
    for i in 0..bays {
        tie_edges.push(edges.len());
        edges.push((tie(i), tie(i + 1))); // tie: tension
        q.push(logu(&mut state, 6.0, 12.0));
    }
    for i in 1..bays {
        edges.push((i, tie(i))); // hangers: tension
        q.push(logu(&mut state, 0.5, 1.5));
    }
    let fixed = vec![0, bays];
    let mut free: Vec<usize> = (1..bays).collect();
    free.extend(bays + 1..2 * bays);
    let loads = free
        .iter()
        .map(|&n| if n <= bays { 0.1 } else { 1.0 })
        .collect();
    Net {
        name: format!("tiedarch{bays}"),
        edges,
        free,
        fixed,
        plan,
        loads,
        q_true: q,
        extent: span,
        tie_edges,
    }
}

// ───────────────────────── registry ─────────────────────────

/// Nets used by the `suite` subcommand, in order.
pub fn suite_nets() -> Vec<Net> {
    vec![
        quad(21, true, 1.0).normalized(0.25),
        quad(21, true, 1.0).normalized(1.0).named("quad21c_d1"),
        quad(21, true, -1.0).normalized(0.25),
        quad_holes(21, true, 0.15).normalized(0.25),
        diamond(10).normalized(0.25),
        radial(8, 24, false).normalized(0.25),
        oculus(21, 7).normalized(0.25),
        crease(21).normalized(0.25),
        crease_diag(21).normalized(0.25),
        hypar_mixed(21),
        wheel(24, 4).normalized(0.25),
        tied_arch(16).normalized(0.25),
        cable_truss(16).normalized(0.25),
        cable_dome(4, 16).normalized(0.25),
        barrel_vault(16, 12).normalized(0.25),
    ]
}

/// Every net (the suite plus compression mirrors), for the `nets` sanity check.
pub fn all_nets() -> Vec<Net> {
    let mut nets = suite_nets();
    nets.push(radial(8, 24, false).normalized(0.25).comp());
    nets.push(oculus(21, 7).normalized(0.25).comp());
    nets
}
