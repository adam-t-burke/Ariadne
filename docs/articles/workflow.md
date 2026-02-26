# Workflow

A typical Ariadne workflow: define the network, then solve (forward or optimization).

## High-level flow

1. **Topology and geometry** — Build the graph (nodes, edges, fixed nodes) and assemble the FDM network.
2. **Initial q and loads** — Set force densities (e.g. uniform or by layer) and loads on free nodes.
3. **Solve** — Use **Theseus Solve** for a single equilibrium (forward) or connect **Opt Config** for inverse form-finding.
4. **Use results** — Export geometry, forces, lengths; optionally feed optimized **Q** back or into visualization.

## Component roles (visual)

The main design components under the **Ariadne → Design** tab include:

| Step | Component | Role |
|------|------------|------|
| Graph | Construct Graph / Deconstruct Graph | Build connectivity and node roles (free/fixed). |
| Network | Construct Network | Build FDM network from graph and geometry. |
| Q | Q by Layer (optional) | Assign initial force densities by layer or grouping. |
| Solve | **Theseus Solve** | Forward solve or optimization; outputs network, nodes, edges, Q, forces, lengths. |
| Config | **Optimization Config** | Objectives, q bounds, tolerances, Run, Stream Preview. |

Objective components (node and edge objectives) feed into **Optimization Config**; connect their outputs to the **Objectives** input of Opt Config.

For a **workflow diagram with component icons**, the docs site shows the same sequence with small icons next to each component name. Icons are in the repository under `Resources/` and are copied to `docs/assets/icons/` when building the site (see [Component guide](component-guide.md) for per-component icons and descriptions).

## Forward vs optimization

- **Forward solve** — Opt Config input left unconnected. One linear solve; fast. Use when you only need equilibrium for a given *q* and loads.
- **Optimization** — Opt Config connected. Solver minimizes objectives over *q* subject to bounds; multiple equilibrium solves. Use for inverse form-finding (target lengths, forces, reactions, etc.).

For more detail on what is being optimized, see [Theory](theory.html).
