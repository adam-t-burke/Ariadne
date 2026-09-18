# Ariadne

[![API reference](https://img.shields.io/badge/API-reference-8ca0b4?style=flat-square)](https://adam-t-burke.github.io/Ariadne/api/)

Ariadne is a Grasshopper plugin for **inverse design of form-found structures** in Rhino 8. It uses the force density method (FDM) and the Theseus solver for directed form-finding of tensile networks.

## Features

- **FDM networks** — Build and solve cable networks (forward solve and inverse optimization).
- **Optimization** — Minimize objectives (target lengths, forces, reactions) over force densities with configurable bounds and tolerances.
- **Streaming preview** — View intermediate results during optimization or output only the final result.
- **Design workflow** — Graph and network construction, Q-by-layer assignment, and visualization utilities.

## Requirements

- **Rhino 8** with Grasshopper
- **.NET 8** and **Rust 1.87 or newer** (only when building from source)

## Installation

**Recommended:** Install via the Package Manager in Rhino (search for *Ariadne*).

**From source on Windows:** Clone the repo, run `.\build.ps1`, then open
`Ariadne.sln` in Visual Studio and build in Release. Copy `Ariadne.gha`,
`theseus.dll`, `Ariadne-LICENSE.txt`, `NATIVE_ARTIFACTS.md`, and every
`basin-*` and `ariadne-lbfgsb-*` notice from `bin\Release\net8.0\` into your Grasshopper
Libraries folder. On macOS, run `./build.sh release`; it replaces the tracked
bootstrap dylib with a universal binary built from the checked-out source.
Release packages should use CI-built native artifacts, not untouched bootstrap
binaries from the source tree.

## Quick start

Define an FDM network (nodes, branches, fixed nodes), set initial force densities and loads, and connect **Theseus Solve**. For equilibrium only, leave the Opt Config input disconnected. For inverse form-finding, connect **Optimization Config** with objectives and bounds.

## Solver usage notes

- **Thread safety:** Each `TheseusSolver` instance is independent. Do not share a handle across threads; create one per solve.
- **Memory:** `TheseusSolver` implements `IDisposable`. Always use `using` or call `Dispose()` to free the native handle.
- **Layout:** All flattened arrays (xyz, loads, targets) use row-major layout: `array[index * 3 + dim]`.
- **Native library:** `theseus.dll` / `libtheseus.dylib` is built from the Rust workspace under `crates/` via `build.ps1` or `build.sh`.
- **q bounds:** Optimization defaults to `DirectBoxBounds` (Basin L-BFGS-B; strict finite lower/upper bounds on every edge). Switch to `DirectSoftBounds` from the OptConfig right-click menu for one-sided or infinite bounds. Legacy saved configs using the removed implicit bounded mode are migrated to box bounds when all q bounds are finite, otherwise to soft bounds.
- **Convergence tolerances:** In `DirectBoxBounds`, **AbsTol** bounds the projected-gradient infinity norm and **RelTol** bounds the relative accepted-iterate reduction `(f[k]-f[k+1])/max(|f[k]|,|f[k+1]|,1)`. Initial parameters outside finite bounds are projected before the first objective evaluation. `DirectSoftBounds` uses **AbsTol** for the Euclidean gradient norm and **RelTol** for the absolute change in objective value, preserving the existing settings behavior. Zero disables either tolerance.
- **Iteration limit:** Positive **MaxIter** values are applied directly as the maximum number of accepted optimization iterations; no smaller hidden default cap is applied.
- **Optimizer comparison:** See the [Basin benchmark report](crates/theseus/BENCHMARKS.md) for reproducible timings, numerical results, and stopping criteria against the previous solvers.
- **Variable supports:** In `DirectBoxBounds`, roller, rail, and NURBS parameters use normalized hard boxes without sigmoid saturation. `DirectSoftBounds` keeps the smooth unbounded latent maps.
- **q cache:** With **Cache Q** enabled, every completed solve with finite force densities updates the cache, including useful best-so-far results that stop at the iteration limit. Use **Reset Cache** to return to the input q values.
- **Hydrostatic follower loads:** `Pressure Load (Hydrostatic)` applies `p = ρgh` along each current face normal. For a forward solve with fixed positive `q`, Theseus uses Newton iterations with load continuation and returns the loaded equilibrium with member force `F = qL`. Face winding sets the pressure direction.

### Recovering unloaded formwork from a loaded hydrostatic state

The hydrostatic FDM result is a loaded funicular geometry, not a fabric cutting
pattern. To estimate an unloaded formwork state:

1. Solve the hydrostatic equilibrium with the chosen positive force densities.
2. Read each loaded member length `L` and force `F = qL`.
3. Choose a material law. For a linear axial member,
   `F = EA(L - L₀)/L₀`, so `L₀ = L/(1 + F/EA)`.
4. Use the recovered rest lengths `L₀` as constraints for constructing the
   unloaded 3D formwork. Flattening panels to 2D is a separate cutting-pattern
   problem and requires a membrane/material model.

If the loaded shape must match a prescribed target, first optimize `q` against
target-position objectives with the hydrostatic load connected, then recover
rest lengths from that matched loaded state.

## Linear solver backends

Every forward solve and every gradient evaluation solves the sparse FDM system
`A x = b` (and its adjoint) for three right-hand sides. Theseus offers three
backends for that step, chosen explicitly per definition; nothing is selected
automatically and the solver never falls back from one backend to another.

| Kind | Menu label | Message | Status |
| --- | --- | --- | --- |
| `Direct` (default) | Direct | `lin: Direct` | Available. Sparse Cholesky factorization; exact, robust, the right choice up to a few tens of thousands of nodes. |
| `IterativeCpu` | Iterative (CPU) | `lin: Iter-CPU` | **Not yet available in this build.** Matrix-free flexible CG with an aggregation-AMG preconditioner, for large meshes where factorization time or memory dominates. |
| `IterativeGpu` | Iterative (GPU) | `lin: Iter-GPU` | **Not yet available in this build.** The same method on a wgpu compute device (Vulkan, DirectX 12, Metal) with a host or device f64 outer loop. |

**Enabling a backend.** Right-click **Optimization Config** → **Linear solver**
and pick a kind. The choice is saved with the definition (key
`LinearSolverKind`; definitions saved before the toggle existed load as
`Direct`) and shows up in the component message as `lin: Direct`,
`lin: Iter-CPU` or `lin: Iter-GPU`. Selecting an iterative kind adds an
optional **Iterative Options** input; feed it from the **Iterative Solver
Options** component to change the tolerance schedule, iteration budget,
multigrid cycle, smoother degree, aggregation passes, coarsest level size,
preconditioner precision or GPU adapter preference. Leave it empty for the
defaults (adaptive tolerance `clamp(1e-2 × projected-gradient ratio, 1e-10,
1e-6)`, 200 iterations, K-cycle, degree-2 Chebyshev smoother, two aggregation
passes, coarsest level of 2000 nodes, backend-default precision). The same
knobs are exposed to C# through `TheseusSolver.SetLinearSolver` and
`TheseusSolver.SetIterativeOptions`, and to the config record through
`OptimizationConfig.LinearSolver` / `.IterativeOptions` (default `Direct`,
`null`).

**Reading the diagnostics.** After every solve **Theseus Solve** adds a remark
of the form `Linear solver: IterativeCpu, 12 solve(s), 340 iteration(s),
15.0 ms`. The same totals are available in C# as `SolverResult.LinearSolver`
(`LinearSolverStats`: backend, solve count, total and maximum CG iterations,
solve and setup milliseconds, whether every solve converged) and to native
callers through `theseus_get_linear_solver_stats`. For `Direct` the iteration
counts are always zero; in this build the direct factorization does not yet
report per-solve statistics either, so a `Direct` run shows
`Linear solver: Direct, 0 solve(s), 0 iteration(s)`.

**Failures are loud.** Selecting `Iterative (CPU)` or `Iterative (GPU)` in this
build makes the next solve fail with the **red** component error
`Linear solver: iterative linear solvers are not yet available in this build
(requested 'IterativeCpu') (native code -4)` (`TheseusException.NativeCode ==
-4`). Once the backends land, `Iterative (GPU)` first runs a GPU probe
(`GpuProbe.Query()` / `theseus_gpu_probe`) and, when no usable adapter is
found, fails the same way with the probe's reason and the list of adapters it
saw (native code -5). Iterative solves that do not converge within their
budget report native code -3 and out-of-memory on the device -6. In every case
the component shows an error instead of silently switching to `Direct`; switch
back to **Direct** in the context menu to continue.

## Objective reference

| C# method | Description |
|-----------|-------------|
| `AddTargetXyz` | Minimize 3D distance to target positions (SSE / MSE / RMSE) |
| `AddTargetXy` | Minimize 2D (XY) distance to target positions (SSE / MSE / RMSE) |
| `AddTargetPlane` | Minimize (u,v) distance to target positions on a plane (SSE / MSE / RMSE) |
| `AddPlanarConstraintAlongDirection` | Pull nodes onto a plane along a direction |
| `AddTargetLength` | Minimize difference from target edge lengths |
| `AddLengthVariation` | Minimize range of edge lengths |
| `AddForceVariation` | Minimize range of member forces |
| `AddSumForceLength` | Minimize sum of force × length products |
| `AddMinLength` / `AddMaxLength` | Barrier penalties for edge length thresholds |
| `AddMinForce` / `AddMaxForce` | Barrier penalties for force thresholds |
| `AddRigidSetCompare` | Compare pairwise distances between node sets |
| `AddReactionDirection` | Align anchor reaction directions |
| `AddReactionDirectionMagnitude` | Align reaction directions and magnitudes |
| `AddReactionMagnitude` | Match reaction magnitudes with configurable behavior |

## License

The Ariadne project is MIT licensed; see [LICENSE.txt](LICENSE.txt). Theseus uses
Basin, available under MIT or Apache-2.0, with preserved
[third-party notices](third-party/basin/README.md). The
standalone translated solver under `crates/lbfgsb` is BSD-3-Clause and carries
its own [license](crates/lbfgsb/LICENSE), preserved
[upstream license](crates/lbfgsb/UPSTREAM_LICENSE.txt), and
[third-party notices](crates/lbfgsb/THIRD_PARTY_NOTICES.md). Native
redistributions must include the generated license/notice bundle beside the
DLL or dylib.

**Author:** Adam Burke — [aburke3@mit.edu](mailto:aburke3@mit.edu)
