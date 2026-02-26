# Component guide

Overview of Ariadne Grasshopper components. Icons are 24×24 and live in the docs under `assets/icons/` (sourced from the plugin’s `Resources/` when the site is built).

## Design tab (Ariadne → Design)

### Theseus Solve

- **Icon:** `Create.png`
- **Description:** Solves the FDM network using the Theseus engine. Performs a **forward solve** (single equilibrium) when Opt Config is disconnected, or **optimization** over force densities when Opt Config is connected.
- **Inputs:** Network, Force Densities (q), Loads, Opt Config (optional).
- **Outputs:** Network, Nodes, Edges, Force Densities (Q), Member Forces, Member Lengths, Iterations, Convex.

Use a **Run** trigger (e.g. from Opt Config) to start optimization; with Stream Preview on, intermediate results are streamed to the outputs.

### Optimization Config (Opt Config)

- **Icon:** `parameters.png`
- **Description:** Bundles optimization settings for Theseus Solve: objectives, q bounds, tolerances, max iterations, barrier parameters, Report Frequency, Run toggle, Stream Preview.
- **Inputs:** Objectives, Lower Bounds (qMin), Upper Bounds (qMax), Max Iterations, AbsTol, RelTol, Barrier Weight, Barrier Sharpness, Report Frequency, Run, Stream Preview.
- **Output:** Config (connect to Theseus Solve’s Opt Config input).

### Network construction

- **Construct Network** — Builds the FDM network from graph and geometry. Icon: `Construct_Network.png`.
- **Deconstruct Network** — Breaks the network back into parts. Icon: `Deconstruct_Network.png`.

### Graph

- **Construct Graph** — Icon: `Construct_Graph.png`.
- **Deconstruct Graph** — Icon: `Deconstruct_Graph.png`.

### Q by Layer

- **Icon:** `Qbylayer.png`  
- Assigns initial force densities by layer or grouping.

## Objectives (node and edge)

Connect these to the **Objectives** input of **Optimization Config**. They define what the solver minimizes (e.g. deviation from target length, target force, or reaction targets).

### Edge objectives (examples)

- Length variance, force variance, target length, min/max length, min/max force, performance-style terms. Icons: `lengthvar.png`, `Forcevar.png`, `Target_Length.png`, `minlength.png`, `maxlength.png`, `minforce.png`, `maxforce.png`, `performance.png`.

### Node objectives (examples)

- Target (point), Target XY, Target UV, rigid point set, reactions. Icons: `Target.png`, `Target_XY.png`, `Target_UV.png`, `PointSet.png`, `Reactions.png`.

## Utilities

- **Visualize** — `visualize.png`
- **Network Info** — `Information.png`
- **Dummy Network** — `Dummy_Network.png`
- **Serialize / export** — `Serialize.png`, `export.png`
- **Update Refs** — `update.png`
- **Tagger** — `tagger.png`
- **Baker** — `baker.png`
- **Curve to Curve** — `curve2curve.png`

For exact input/output names and types, use the tooltips in Grasshopper or the [API reference](../../api/index.md).
