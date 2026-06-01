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
- **.NET 8** (only when building from source)

## Installation

**Recommended:** Install via the Package Manager in Rhino (search for *Ariadne*).

**From source:** Clone the repo, open `Ariadne.sln` in Visual Studio, and build in Release. Copy `Ariadne.gha` and `theseus.dll` from `bin\Release\net8.0\` into your Grasshopper Libraries folder.

## Quick start

Define an FDM network (nodes, branches, fixed nodes), set initial force densities and loads, and connect **Theseus Solve**. For equilibrium only, leave the Opt Config input disconnected. For inverse form-finding, connect **Optimization Config** with objectives and bounds.

## Solver usage notes

- **Thread safety:** Each `TheseusSolver` instance is independent. Do not share a handle across threads; create one per solve.
- **Memory:** `TheseusSolver` implements `IDisposable`. Always use `using` or call `Dispose()` to free the native handle.
- **Layout:** All flattened arrays (xyz, loads, targets) use row-major layout: `array[index * 3 + dim]`.
- **Native library:** `theseus.dll` / `libtheseus.dylib` is built from the Rust workspace under `crates/` via `build.ps1` or `build.sh`.

## Objective reference

| C# method | Description |
|-----------|-------------|
| `AddTargetXyz` | Minimize 3D distance to target positions |
| `AddTargetXy` | Minimize 2D (XY) distance to target positions |
| `AddTargetPlane` | Minimize (u,v) distance to target positions on a plane |
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

MIT — see [LICENSE.txt](LICENSE.txt).  
**Author:** Adam Burke — [aburke3@mit.edu](mailto:aburke3@mit.edu)
