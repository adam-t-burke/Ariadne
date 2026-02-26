# Ariadne

[![docs: latest](https://img.shields.io/badge/docs-latest-8ca0b4?style=flat-square)](https://adam-t-burke.github.io/Ariadne/)

Ariadne is a Grasshopper plugin for **inverse design of form-found structures** in Rhino 8. It uses the force density method (FDM) and the Theseus solver for directed form-finding of tensile networks. The method is based on Schek (1974); see the [documentation](https://adam-t-burke.github.io/Ariadne/) for theory and optimization details.

## Features

- **FDM networks** — Build and solve cable networks (forward solve and inverse optimization).
- **Optimization** — Minimize objectives (target lengths, forces, reactions) over force densities with configurable bounds and tolerances.
- **Streaming preview** — View intermediate results during optimization or output only the final result.
- **Design workflow** — Graph and network construction, Q-by-layer assignment, and visualization utilities.

## Requirements

- **Rhino 8** with Grasshopper
- **.NET 8** (only when building from source)

## Installation

**Recommended:** Install via the Package Manager in Rhino (search for *Ariadne*).

**From source:** Clone the repo, open `Ariadne.sln` in Visual Studio, and build in Release. Copy `Ariadne.gha` and `theseus.dll` from `bin\Release\net8.0\` into your Grasshopper Libraries folder.

## Quick start

Define an FDM network (nodes, branches, fixed nodes), set initial force densities and loads, and connect **Theseus Solve**. For equilibrium only, leave the Opt Config input disconnected. For inverse form-finding, connect **Optimization Config** with objectives and bounds. See the [documentation](https://adam-t-burke.github.io/Ariadne/) and [examples](examples/) for step-by-step workflows.

## Examples

Example Grasshopper files (`.gh`) are in the [examples](examples/) folder. The [documentation](https://adam-t-burke.github.io/Ariadne/) lists each example with a short description and download link.

## Project structure

| Area      | Description |
|-----------|-------------|
| **Solver** | Theseus Solve, OptConfig, node and edge objectives |
| **FDM**    | Network construction/deconstruction, Q by layer |
| **Graphs** | Graph construction and experimental controls |
| **Utilities** | Visualization, serialization, tagging, helpers |
| **Theseus** | FFI and wrapper for the Rust force-density solver |

The [API reference](https://adam-t-burke.github.io/Ariadne/api/) is generated from the C# source.

## License

MIT — see [LICENSE.txt](LICENSE.txt).  
**Author:** Adam Burke — [aburke3@mit.edu](mailto:aburke3@mit.edu)
