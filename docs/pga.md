# Projective geometric algebra branch

This is an experimental branch, not a Theseus rewrite. PGA is a geometry
language for Euclidean flats. The force-density solve is not.

## Verdict

Yes: keep a long-lived experimental branch. No: do not merge it into `dev`
until a Theseus consumer exists and finite-difference gradients still match.

Plane-based 3D PGA (`G(3,0,1)`, Klein / Gunn `P(R*_{3,0,1})`) is the right
algebra if the goal is points, lines, planes, incidence, and rigid motions.
It is the wrong algebra if the goal is a faster FDM kernel, better `q`
parameterization, or NURBS supports.

## What Ariadne already is

Theseus is a sparse Euclidean FDM solver with hand-coded adjoints:

- coordinates are `R^3` arrays
- variable supports are an enum of *charts*: sphere, axis-aligned roller box,
  rail segment, NURBS curve, NURBS surface
- `PlanarConstraintAlongDirection` is a pull along a vector `d`, not
  orthogonal incidence
- faces already have a Euclidean normal (Newell)
- motors / rigid subassemblies are absent

Those support kinds are mostly *flats plus bounds*. PGA talks about the flats.
The bounds stay optimizer coordinates (`t`, `uv`, box limits, L-BFGS-B).

## What PGA is good for here

1. **Incidence of nodes on lines and planes.** A rail is a point on a line. A
   2-axis roller is a point on a plane. Meet and join replace one-off cross
   products and plane equations, including infinity (ideal line / ideal plane).
2. **Composition of constraints.** Two planes meet in a rail. Three planes meet
   in a point. That composition is already the algebra, rather than another
   `VariableSupportKind` variant.
3. **Oriented geometry for faces.** Three nodes join to a plane; four nodes
   fail to be coplanar by a trivector residual. Useful for planarity
   objectives and for hydrostatic face geometry.
4. **Motors later.** The even subalgebra is dual quaternions. If Ariadne grows
   rigid frames attached to tensile nets, sandwich products are the native
   rigid-motion API. Not needed for today's FDM.
5. **Graphic statics as a research horizon.** Reciprocal diagrams and
   Maxwell–Cremona constructions are projective. PGA is a plausible language
   for form/force polarity. That is a new solver family, not a refactor of
   `solve_fdm`.

## What PGA is a poor fit for

| Theseus piece | Why PGA does not eat it |
| --- | --- |
| Sparse `A(q) x = b` and Cholesky | Euclidean linear algebra on coordinates. |
| L-BFGS-B box charts | Interval bounds on `q`, `t`, `uv`. Not flats. |
| NURBS curve / surface supports | Rational splines. Keep `nurbsbook`. PGA can hold the tangent line or osculating plane, not `C(t)`. |
| Sphere support | A metric ball in a latent chart, not a geometric sphere. Spheres are CGA, which is 32-D and not worth it. |
| `PlanarConstraintAlongDirection` | Pull *along* `d`. PGA incidence is orthogonal distance. Same plane, different residual. |
| Hydrostatic follower loads | Need current face area/normal. Newell already does this. |
| Force densities `q` | Scalars on edges. Weighted lines are a stretch, not a design. |

CGA (`G(4,1)`) is the usual next suggestion. Skip it unless geometric spheres
and circles become first-class supports. VGA (`G(3,0)`) cannot represent
points as algebra elements.

## Mapping from current supports

| Theseus kind | PGA object | Leftover chart |
| --- | --- | --- |
| Fixed node | Point | none |
| Rail | Line through `start`, `end` | `t ∈ (0, 1)` |
| Roller, 1 free axis | Line | axis interval |
| Roller, 2 free axes | Plane | in-plane box |
| Roller, 0 free axes | Point | none |
| Roller, 3 free axes | not a flat | box / free XYZ |
| Sphere | not a flat | radial box |
| NURBS curve / surface | not a flat | `t` / `uv` |

The leftover charts are why this crate must not swallow `variable_supports.rs`.
PGA can compute incidence residuals and construct the supporting flat. The
optimizer still needs a 1-D or 2-D parameterization with bounds.

## Branch policy

Keep the crate standalone (`publish = false`), same pattern as `nurbsbook`
before it was consumed by Theseus.

Do in this branch:

- G(3,0,1) kernel with Euclidean oracles
- typed Point / Line / Plane / Motor
- explicit mapping tests against Ariadne primitives
- no FFI, no Grasshopper components, no `theseus` dependency

Do next, still on this branch, only if the kernel stays stable:

- optional Theseus objective: orthogonal point-on-plane / point-on-line
  incidence, *alongside* `PlanarConstraintAlongDirection`, with FD gradient
  checks
- rail construction as `join(start, end)` inside a pure helper, without
  changing the latent chart

Do not do:

- rewrite `fdm.rs` in GA
- replace NURBS
- add GH components before a solver consumer
- generic n-D codegen or a full Klein SIMD port

## Risks

- **Dual vs metric dual.** Join must use the Poincaré complement. Tests that
  compare against Euclidean incidence are the contract.
- **False unification.** Axis-aligned roller *boxes* are not planes. Calling
  them PGA planes would drop the bounds.
- **Gradient surface area.** Any Theseus wiring needs hand-coded adjoints, same
  as every other objective. PGA residuals are not free in the optimizer.
- **Graphic statics scope creep.** Reciprocal diagrams deserve their own
  crate/branch after incidence is real.

## Recommendation

Treat PGA as the geometry layer for flats, the same way `nurbsbook` is the
geometry layer for splines. Leave Theseus as the solver. Merge only after an
objective or support helper uses this crate and the existing FD tests still
pass.
