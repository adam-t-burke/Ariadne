# Ariadne

Ariadne is a Grasshopper plugin for the **inverse design of form-found structures in Rhino 8**. It uses the **force density method** (FDM) via the Theseus solver engine for directed form-finding of tensile network structures.

**Reference:** Schek, H.-J. "The force density method for form finding and computation of general networks." *Computer Methods in Applied Mechanics and Engineering*, Vol. 3, No. 1, 1974, pp. 115–134. For a mathematical introduction and an overview of optimization over force densities, see the [Theory](https://adam-t-burke.github.io/Ariadne/articles/theory.html) section in the documentation.

## Features

- **FDM networks** — Build and solve cable/networks with the force density method (forward solve and inverse optimization).
- **Optimization** — Minimize objectives (target lengths, forces, reactions, performance) over force densities with configurable bounds and convergence tolerances.
- **Streaming preview** — See intermediate results during optimization; optional single-shot output.
- **Design workflow** — Graph construction, network assembly, Q by layer, and visualization utilities.

## Requirements

- **Rhino 8** (with Grasshopper).
- **.NET 8** (only if building from source).

## Installation

1. **From Rhino (recommended)** — Open the Package Manager in Rhino and search for **Ariadne**; install from there.
2. **From source** — Clone this repository, open `Ariadne.sln` in Visual Studio, and build in Release. Copy the built `Ariadne.gha` and `theseus.dll` from `bin\Release\net8.0\` into your Grasshopper Libraries folder (or use the output path as a library path in Grasshopper).

## Quick start

Define an FDM network (nodes, branches, fixed nodes), supply initial force densities and loads, then connect **Theseus Solve**. For equilibrium-only (forward solve), leave the Opt Config input disconnected. For inverse form-finding, connect an **Optimization Config** with objectives and bounds; the solver will optimize force densities to meet your criteria. See the [documentation](https://adam-t-burke.github.io/Ariadne/) and [examples](examples/) for detailed workflows.

## Documentation

Full documentation (theory, getting started, workflow, component guide, API reference) is published at:

**https://adam-t-burke.github.io/Ariadne/**

## Examples

Example Grasshopper files (`.gh`) are in the [examples](examples/) folder. Each file is a how-to you can open in Rhino; the [documentation](https://adam-t-burke.github.io/Ariadne/) lists them with short descriptions and download links.

## Project structure

- **Solver** — Theseus Solve component, OptConfig, objectives (node and edge).
- **FDM** — Network construction/deconstruction, Q by layer.
- **Graphs** — Graph construction and experimental controls.
- **Utilities** — Visualization, serialization, tagging, and helpers.
- **Theseus** — FFI and wrapper around the Rust force-density solver.

The [API reference](https://adam-t-burke.github.io/Ariadne/api/) is generated from the C# source.

## License and author

**License:** MIT — see [LICENSE.txt](LICENSE.txt).

**Author:** Adam Burke — [aburke3@mit.edu](mailto:aburke3@mit.edu)
