# Getting started

## Prerequisites

- **Rhino 8** with Grasshopper.
- **.NET 8** only if you are building Ariadne from source.

## Installation

1. **From Rhino (recommended)** — In Rhino, open the Package Manager and search for **Ariadne**. Install from there.
2. **From source** — Clone the [repository](https://github.com/adam-t-burke/Ariadne), open `Ariadne.sln` in Visual Studio, and build in **Release**. Copy `Ariadne.gha` and `theseus.dll` from `bin\Release\net8.0\` into your Grasshopper Libraries folder.

## First steps

### 1. Define a network

Use Ariadne’s graph and FDM components to define:

- **Nodes** — Free and fixed (support) nodes.
- **Branches** — Connectivity (which nodes each edge connects).
- Optionally, **initial force densities** (e.g. from **Q by Layer** or a list).

The **Construct Network** component (and related graph components) produce an **FDM Network** data structure that the solver expects.

### 2. Forward solve (equilibrium only)

- Connect the **FDM Network** to **Theseus Solve**.
- Supply **Force Densities (q)** and **Loads** (e.g. gravity).
- Leave **Opt Config** disconnected.
- The component solves the equilibrium once and outputs the updated network, node positions, edge curves, force densities, forces, and lengths.

### 3. Inverse optimization (form finding)

- Connect an **Optimization Config** to the **Opt Config** input of **Theseus Solve**.
- In Opt Config, add **objectives** (e.g. target length, target force, reaction targets) and set **bounds** on force densities (*q_min*, *q_max*).
- Set **Run** (or use a button) to start optimization. With **Stream Preview** on, you’ll see intermediate results; the final result appears when the solver converges or hits max iterations.

For the underlying math, see [Theory: Force density method and optimization](theory.md). For a visual workflow and component roles, see [Workflow](workflow.md). For ready-made files, see [Examples](examples.md).
