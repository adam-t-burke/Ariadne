# Compliance-weighted warm starts for inverse force-density form finding

Presentation outline and technical write-up for the Stage-1/Stage-2 warm start
implemented in `crates/theseus/src/inverse.rs` and `inverse_extra.rs`. All
numbers below were produced by `cargo run --release -p theseus --example
warm_start_bench -- <subcommand>` on the branch that introduced this document;
the tables are reproducible with the commands quoted in each section.

---

## 1. Problem statement

Given a network topology, fixed node positions `x_f`, free-node loads `p`, and a
**target** geometry `x*` for the free nodes, find force densities `q` such that
the force-density forward solve lands on the target:

```
D(q) x(q) = p − D_f(q) x_f,      D(q) = Cₙᵀ diag(q) Cₙ,   D_f(q) = Cₙᵀ diag(q) C_f
x(q) ≈ x*   with   lo ≤ q ≤ hi
```

Downstream, a nonlinear optimiser (L-BFGS-B on the forward solve with the hard
box `[lo, hi]`) refines `q` against the full objective. The question addressed
here is **what to hand that optimiser as its starting point**, and at what cost.

Constraints that shaped the design:

* Networks up to ~150 k edges: everything must be sparse, `O(nnz)` memory, and
  should scale roughly like one sparse factorisation per step.
* Tension-only, compression-only, and **mixed-sign** systems (tied arches,
  cable trusses, cable domes, barrel vaults with ties).
* Reaction constraints (e.g. "no horizontal reaction at the supports of a tied
  arch") must be honoured and must be detected when they are inconsistent.

---

## 2. The identity everything rests on

At the target geometry the FDM residual in force-density coordinates is

```
r(q) = E(x*) q − p,      E(x*) = Cₙᵀ diag(Cₙ x* + C_f x_f)   (stacked per axis)
```

and because `E(x*) q = D(q) x* + D_f(q) x_f`, the forward solve satisfies

```
x(q) − x* = − D(q)⁻¹ r(q).
```

So the **geometric error is the force residual pre-conditioned by the
compliance `D(q)⁻¹`**. Two consequences:

1. Minimising `‖r(q)‖` (the classic "force residual" inverse, Schek-style
   least squares on the equilibrium matrix) is *not* minimising the geometric
   error. It weights every node equally in force, which weights stiff nodes too
   much and soft nodes too little. On shallow nets the vertical rows of `E` are
   small (`Δz ≪ Δx`), so the force residual barely sees the height error at all.
2. The correct linear model of the geometric error around `q_k` is

   ```
   x(q_k + Δ) − x* ≈ − D(q_k)⁻¹ ( r(q_k) + E(x(q_k)) Δ )
   ```

   which is a **compliance-weighted least squares (CWLS)** problem in `Δ`. With
   the Jacobian taken at the target, `E(x*)`, it is a *frozen* CWLS step; with
   the Jacobian re-assembled at the current geometry `x(q_k) = x* − D⁻¹r`, it is
   an exact Gauss--Newton step. Both share the same fixed point since
   `E(x(q)) = E(x*)` once the target is hit.

A useful invariance: `D(s·q) = s·D(q)`, so the CWLS weighting depends only on
the *pattern* of `q_k`, not on its scale. This is why a uniform seed with the
right sign pattern is a legitimate metric seed.

---

## 3. Pipeline

```
target x*, loads p, box [lo,hi]
        │
        ▼
 (0) non-dimensionalise ── positions / L, loads / P, q · L/P
        │
        ▼
 (1) Stage 1: force-residual particular  min ‖E(x*)q − p‖²  s.t. box
        │        (Clarabel QP in member-force coordinates; sparse)
        ▼
 (2) seed guard ── score Stage-1 seed vs scaled uniform sign seed on ‖x(q)−x*‖
        │           suspicious?  → race both through the frozen step
        ▼
 (3) frozen CWLS step(s)   min ‖D(q_k)⁻¹(E(x*)Δ + r_k)‖² + Σ Λ_j Δ_j²  s.t. box
        │                  (active-set BVLS on the sparse weighted saddle)
        ▼
 (4) damped Gauss--Newton step(s)   same, Jacobian at x(q_k), LM damping
        │                           accept only if the measured error drops
        ▼
 warm-start q  →  L-BFGS-B (DirectBoxBounds)
```

### 3.1 Non-dimensionalisation (`nondimensionalize`)

Positions are divided by the bounding-box diagonal `L`, loads by `P = max|p|`;
force densities transform as `q̃ = q·L/P`, the box likewise, and the Tikhonov
terms as `λ̃ = λ/L²` (Stage 1) and `λ̃_c = λ_c P²/L⁴` (Stage 2) so that every
objective is the original divided by a constant. The recovered `q` is identical
in exact arithmetic; what changes is the conditioning of the KKT systems handed
to LDL and Clarabel, which previously failed in millimetre units. It also makes
the relative weight of reaction rows against geometric rows unit-free.

### 3.2 Stage 1

Unchanged in spirit: the box-constrained least-squares particular of the
equilibrium matrix at the target, in member-force coordinates (`solve_for_q =
false`) so that the box is expressed on `t = q·ℓ` with the target lengths. It is
the only place in the pipeline where the self-stress pattern of a self-tied
system can be recovered, and it is cheap (one interior-point solve on a sparse
system with `3n_free (+ reaction rows) × n_edges` entries).

### 3.3 Seed guard (`seed_guard_margin`, default 3)

Stage 1 collapses on shallow targets with white-noise perturbations: it trades
height error against in-plane force balance, and returns a seed whose forward
geometry can be an order of magnitude further from the target than a uniform
`q` would be. The guard builds a uniform seed with the sign pattern of the box
(or of Stage 1 where the box allows both signs), fits its magnitude per sign
group by a golden-section search on the exact geometric error (20 Laplacian
factorisations per group), and compares. If Stage 1 is worse by more than the
margin it is *suspect*; both seeds then take the frozen step and the lower
measured error continues. Tested effect: on the loose-box shallow-quad jitter
case the warm start improves from `0.179·L` to `0.108·L` and L-BFGS-B reaches
5 % of its final loss in 110 evaluations instead of 427; on cases where Stage 1
is fine the race picks Stage 1 and nothing changes.

### 3.4 Frozen CWLS step

Solves the weighted least squares as a sparse saddle system

```
[ I      0      Sᵀ ] [ e ]   [ 0  ]
[ 0      Λ    −Eᵀ  ] [ Δ ] = [ 0  ]        S = blkdiag(D, D, D, I_reactions)
[ S     −E      0  ] [ y ]   [ r  ]
```

with a single sparse LDL factorisation (faer). Nothing dense is ever formed; the
saddle has `2m + n` rows with `m = 3n_free + n_reaction_rows` and `nnz ≈
nnz(D)·3 + 2·nnz(E) + m + n`.

### 3.5 Bounds inside Stage 2: active-set BVLS (`Stage2Method::ActiveSet`)

The box on `Δ` is handled by an add/release active set:

1. Solve the saddle with the currently bound columns held at their bound
   (their coupling entries are kept in the pattern but zeroed and their
   contribution moved to the right-hand side, so the **symbolic factorisation is
   reused** across passes and across outer steps).
2. Any free variable that left the box is fixed at the violated bound; repeat.
3. When nothing violates, evaluate the reduced gradient
   `∇f_j = Λ_j Δ_j − E_jᵀ y` on the bound variables and release those whose
   multiplier points into the box (at most three release rounds).
4. The bound status persists to the next outer step as a warm start. If the set
   does not settle in 24 passes, the step falls back to Clarabel.

For healthy seeds this is two or three numeric refactorisations instead of a
15–30-iteration interior-point solve, and it is exact (no barrier smoothing of
the bound). Clarabel remains selectable and is the fallback.

### 3.6 Damped Gauss--Newton (`lm_damping`, default 1e-4)

The exact Gauss--Newton step from a poor seed overshoots (a q that goes to zero
on part of the net makes `D` nearly singular and the step huge). Instead of
halving the step length, the direction is damped Levenberg--Marquardt style:
`Λ_j = λ_c + λ_LM · s_j` with `s_j ≈ diag(Jᵀ S⁻² J)_j` estimated from the local
resistance of edge `j` (`‖E_j‖²` divided by the Laplacian diagonal at its end
nodes), `λ_LM` multiplied by ten on a rejected step and divided by ten on an
accepted one. Every candidate is scored by an exact probe `‖D(q)⁻¹ r(q)‖`, so a
step is never accepted on faith.

### 3.7 Reaction constraints (`enforce_zero_r{x,y,z}`)

`enforce_zero_rx` appends one row per fixed node requiring its x-reaction to
vanish. For a pin–pin tied arch this is what forces the thrust into the tie
instead of into the supports. Because the rows pin *every* support, they are
consistent only if the applied load has no net component along that axis; the
solver now rejects `enforce_zero_rz` under vertical load with an explicit error
instead of silently returning a collapsed particular.

---

## 4. Benchmarks

*(Filled from the `warm_start_bench` runs; see the sections below.)*

---

## 5. Sparsity and scaling

*(Filled from `warm_start_bench scale` and `warm_start_bench dense`.)*

---

## 6. Related approaches and how they compare

*(Filled below.)*

---

## 7. Recommendations and open items

*(Filled below.)*
