# Frozen CWLS / Gauss–Newton versus short L-BFGS-B from `q*`

This note answers a concrete question about the Stage-2 warm start in
[`warm_start.md`](warm_start.md): after the Stage-1 force-residual particular
`q*` — which can sit at a *large* geometric residual — is a handful of
compliance-weighted linearisations worth the cost, or does handing `q*` to
L-BFGS-B on the exact target-geometry objective close the same gap in ten
accepted steps (enough to fill the compact BFGS history)?

The numerical evidence is produced by

```bash
cargo run --release -p theseus --example warm_start_bench -- tradeoff bench/figures/data
cd bench/figures && uv run render_tradeoff.py
```

and lives in [`bench/figures/data/tradeoff.json`](../bench/figures/data/tradeoff.json).
Figures are under [`figures/`](figures/). Seed guard is **off** for every
controlled linearisation so every method starts from the same `q*`; the
production pipeline with the guard is reported as `pipe_guard` only.

---

## 1. Three maps of the same residual

The forward force-density solve `D(q) x(q) = p − D_f(q) x_f` and the
equilibrium identity `E(x*) q = D(q) x* + D_f(q) x_f` give

```
r(q)  := E(x*) q − p                 (affine in q)
x(q) − x* = − D(q)⁻¹ r(q)            (exact for geometry-independent loads)
```

so the geometric residual *is* the force residual seen through the current
compliance. Differentiating the implicit residual `D(q) x + D_f(q) x_f − p = 0`
produces the true Jacobian of the forward map:

```
J(q) := ∂x/∂q = − D(q)⁻¹ E(x(q)).
```

`E` is assembled on the *current* edge directions `C_n x(q) + C_f x_f`. The
identity residual uses `E(x*)`. Those two matrices agree if and only if
`x(q) = x*`.

The geometric objective L-BFGS-B actually minimises is

```
f(q) = ½ ‖x(q) − x*‖² = ½ ‖D(q)⁻¹ r(q)‖².
```

Its gradient (the adjoint already implemented in Theseus) is exact:

```
∇f(q) = J(q)ᵀ (x(q) − x*)
      = E(x(q))ᵀ D(q)⁻ᵀ D(q)⁻¹ r(q)
      = E(x(q))ᵀ D(q)⁻² r(q).
```

Three first-order models of a step `q_k + Δ` now sit on the table.

### 1.1 Stage 1 — Euclidean force residual

```
q* = arg min  ‖E(x*) q − p‖²    s.t.  lo ≤ q ≤ hi.
```

This is the correct *pattern* of `q` (and the only cheap place the self-stress
of a self-tied net can be recovered), but it is the wrong *metric*. Nodes that
`D(q)` treats as soft are under-weighted; the shallow vertical rows of `E(x*)`
are almost invisible. The geometric error at `q*` is `‖D(q*)⁻¹ r(q*)‖`, which
can be orders of magnitude larger than `‖r(q*)‖` suggests.

### 1.2 Frozen CWLS — Gauss–Newton with `E` held at the target

Hold `D` at `q_k` and `E` at `x*`:

```
x(q_k + Δ) − x*  ≈  − D(q_k)⁻¹ ( r(q_k) + E(x*) Δ ).
```

The boxed least-squares problem in `Δ` is a sparse saddle. It is *exactly*
Stage 1 with the norm changed from `I` to `D(q_k)⁻²`. Because the identity
uses `E(x*)`, this linearisation is the derivative of the identity residual
with the metric frozen — not the derivative of `x(q)`. When `x(q*)` is
collapsed, `E(x*)` still has the target's well-shaped edge directions, which
is why one frozen step can rescue a seed that the true Jacobian at `x(q*)`
cannot.

### 1.3 Gauss–Newton — the true first-order model of `x(q)`

```
x(q_k + Δ) − x*  ≈  − D(q_k)⁻¹ ( r(q_k) + E(x(q_k)) Δ )
                 =  e_k + J(q_k) Δ.
```

The normal equations are the Gauss–Newton system for `f`:

```
( Jᵀ J ) Δ  =  − Jᵀ e_k
E(x)ᵀ D⁻² E(x) Δ  =  − E(x)ᵀ D⁻² r_k.
```

The neglected Hessian term `Σ_i e_i ∇² x_i` is large precisely when the
geometric residual is large, so a full Newton step on `f` is *not* the same
as GN, and GN itself can overshoot — hence the merit line search (step
halving) in Stage 2.

### 1.4 L-BFGS-B — limited-memory inverse Hessian of the *true* `f`

L-BFGS-B uses the exact `∇f` above (so `E` *is* assembled at `x(q)`, the same
Jacobian GN uses) and a compact history of `m = 10` pairs `(s_i, y_i)`.

* Iteration 0 has no history. The first search direction is scaled steepest
  descent `Δ = −α ∇f` in Euclidean `q`-space, not the GN Newton step
  `−(Jᵀ J)⁻¹ ∇f`. The compliance metric `D⁻²` is present in the *gradient*
  but is not inverted.
* After `k ≤ 10` accepted steps the inverse-Hessian approximation has rank at
  most `k` plus a scaled identity. After 10 accepted steps the history is
  saturated: further steps only replace the oldest pair.
* The pairs sample the *true* Hessian `∇²f = Jᵀ J + Σ e_i ∇² x_i` along the
  nonlinear trajectory, including the residual term GN drops. That can help
  once the residual is small and hurt while `x(q)` is still collapsed.
* Each evaluation is one Laplacian factor of `D` (size `n_free`) plus an
  adjoint reuse of the same factor, plus line-search extras. A Stage-2 step
  is 2–5 numeric LDLs of a saddle of size `∼ 2·(3 n_free) + n_edges`.

The first L-BFGS step and the frozen CWLS step are therefore different
*directions*, not just different lengths. Frozen / GN invert a structured
`n_edges × n_edges` Gramian; the first L-BFGS step cannot.

---

## 2. What the comparison is allowed to claim

Write `e(q) = ‖x(q) − x*‖ / L`. Hypotheses, to be read off the suite:

**H1.** Ten L-BFGS-B steps from `q*` do **not** match one frozen CWLS step on
mixed-sign or collapsed Stage-1 seeds. The first L-BFGS direction is
`−∇f(q*)`, which uses `E(x(q*))`; if that geometry is folded the gradient is
a linearisation of the wrong mesh. Frozen uses `E(x*)`.

**H2.** On a well-behaved single-sign hanging net, where `q*` is already in a
mildly nonlinear basin and the missing scale of `q` is low-dimensional, ten
L-BFGS-B steps *can* catch a frozen step and sometimes a GN step, because a
rank-10 Hessian is enough for a handful of soft modes.

**H3.** Pure GN from `q*` (no frozen first) is not uniformly better than
frozen. When `x(q*)` is collapsed, `E(x(q*))` is a bad Jacobian
(`degenerate_linearizations` in the diagnostics). Frozen is the safer first
linearisation; GN earns its keep *after* the geometry has moved toward `x*`.

**H4.** Frozen-then-10-L-BFGS is the interesting cost tradeoff against
frozen-then-2-GN. After one CWLS step the residual is often small enough that
L-BFGS is in the quadratic basin where a rank-10 history is a legitimate
quasi-Newton method. That is a different claim from “skip the linearisation
entirely”.

**H5.** Force residual `‖r‖` and geometric residual `‖D⁻¹ r‖` can move in
opposite directions. Stage 1 is near-optimal for `‖r‖`; every geometric
method is allowed to *increase* `‖r‖` while decreasing `e(q)`. Reporting
only one of them hides the metric change.

---

## 3. Experiment

For each case below, Stage 1 is run once. From that `q*`:

| label | what |
|---|---|
| `s1` | Stage-1 particular, clipped |
| `frozen1`, `frozen2` | 1–2 frozen CWLS steps, guard off |
| `gn1`, `gn2`, `gn3` | 1–3 Gauss–Newton steps, **no** frozen phase, guard off |
| `pipe` | 1 frozen + 2 GN, guard off |
| `pipe_guard` | production pipeline (guard 3) |
| `s1+lb{k}` | L-BFGS-B from `q*` for `k ∈ {1,3,10,20}` accepted iterations |
| `frozen1+lb{k}` | the same from the frozen-1 seed |
| `gn1+lb{k}` | the same from the first GN seed |

Every row reports:

* `geom/L`, `max/L`, `rms/L` — Euclidean / nodal-max / RMS geometric error
  over the characteristic length
* `‖r‖`, `‖r‖/‖p‖` — force residual at the *target* equilibrium matrix
* `‖∇f‖∞` — projected gradient of the SSE objective (the L-BFGS merit)
* `active` — edges on a bound
* wall time, Stage-2 factorisations, accepted L-BFGS iterations and
  evaluations (the loss trace is one entry per `value_and_gradient` call,
  including line search)

The compact history size is the solver default `m = 10`
(`crates/theseus/src/optimizer.rs`). Ten accepted steps therefore saturate
it; twenty only recycle the oldest pairs.

Cases: the presentation showcase plus the easy hanging quad `quad21c`
(jittered target, loose box).

---

## 4. Results

*Filled in after `warm_start_bench tradeoff` on this branch. Numbers are
`geom/L` unless noted.*

### 4.1 Cross-case geometric error

<!-- TRADEOFF_SUMMARY -->

*(table inserted after the run)*

### 4.2 Residual anatomy on a well-behaved net and a collapsed one

<!-- TRADEOFF_DETAIL -->

*(per-case residual tables inserted after the run)*

### 4.3 Direction cosine of the first step

The cosine between `q_frozen − q*` and `q_{L-BFGS,1} − q*` tests §1.4: if
the first L-BFGS step were a damped CWLS step the cosine would be near 1.

<!-- TRADEOFF_COSINE -->

---

## 5. Reading the tradeoff

Cost, qualitatively, on these 400–800-edge nets:

* one frozen / GN step: a few numeric LDLs of the weighted saddle, typically
  a handful of milliseconds, but the constant is the *larger* KKT and the
  active-set passes
* one L-BFGS-B iteration: one `D`-factor + adjoint + 0–few line-search
  forwards, typically cheaper per iteration, 10–20 of them to fill history

The question is not “is L-BFGS cheaper per iteration” — it is — but whether
ten of those iterations substitute for the change of metric that one CWLS
step performs in closed form.

Predicted outcome, to be confirmed or refuted by §4:

1. Skipping the linearisation and taking 10 L-BFGS-B steps from `q*` is
   **not** a substitute for one frozen CWLS step on mixed-sign / collapsed
   seeds. The gradient at a bad `x(q*)` does not invert `D⁻²` and does not
   see `E(x*)`.
2. On a shallow hanging quad, 10 L-BFGS-B steps from `q*` can land close to
   a frozen or GN warm start, because the missing correction is
   low-dimensional.
3. Two GN steps after a frozen step still win on the hardest seeds (cable
   dome, hypar, truss). Ten L-BFGS-B steps *after* the frozen step are the
   fair “trade the remaining GN budget for quasi-Newton” comparison, and
   that comparison can go either way once the geometry is no longer
   collapsed.
4. Force residual at `q*` is the wrong stopping test. Geometric methods
   routinely raise `‖r‖` while cutting `e(q)` by 10–100×.

---

## 6. Regenerating

```bash
cargo build --release -p theseus --example warm_start_bench
B=target/release/examples/warm_start_bench
$B tradeoff bench/figures/data          # writes tradeoff.json, prints residual tables
cd bench/figures && uv run render_tradeoff.py
```

`BENCH_NETS=cabledome4x16,hypar21m` restricts the run.
