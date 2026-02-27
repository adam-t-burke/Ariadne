# Grasshopper / Ariadne Integration Guide

This document describes how to replace the Julia WebSocket communication path
in Ariadne with direct P/Invoke calls to `theseus.dll`.

## Architecture Change

### Before (Julia WebSocket)

```
Grasshopper Component
  → serialize problem to JSON
  → send over WebSocket to Julia server (port 2000)
  → Julia parses JSON, runs FDM optimization
  → Julia sends back JSON result {X, Y, Z, Q, Loss, Losstrace, ...}
  → Grasshopper parses JSON response
```

### After (Rust DLL via P/Invoke)

```
Grasshopper Component
  → build incidence COO arrays from network topology
  → call TheseusSolver.Create(...)
  → add objectives (AddTargetXyz, AddLengthVariation, etc.)
  → call solver.Optimize() or solver.OptimizeConstrained()
  → read results directly from SolverResult struct
```

No network, no JSON serialization, no Julia runtime.

## Step-by-step Migration

### 1. Add `Theseus.Interop` to your Grasshopper project

Reference the `Theseus.Interop` class library project, or copy
`TheseusInterop.cs` and `TheseusSolver.cs` into your existing project.

### 2. Deploy theseus.dll

Place `theseus.dll` (from `rust/target/release/`) alongside your `.gha` file.
The .NET runtime will search the assembly's directory for native DLLs.

Alternatively, use `NativeLibrary.SetDllImportResolver` for custom load paths:

```csharp
using System.Reflection;
using System.Runtime.InteropServices;

NativeLibrary.SetDllImportResolver(
    typeof(TheseusSolver).Assembly,
    (name, assembly, path) =>
    {
        if (name == "theseus")
        {
            var dir = Path.GetDirectoryName(assembly.Location)!;
            return NativeLibrary.Load(Path.Combine(dir, "theseus.dll"));
        }
        return IntPtr.Zero;
    });
```

### 3. Replace problem serialization with direct construction

#### Julia / WebSocket path (current)

Your Grasshopper component currently builds a JSON object like:

```json
{
  "Edges": [[0,1],[1,2],...],
  "FreeNodes": [1,2,3,...],
  "FixedNodes": [0,6],
  "Loads": [[0,0,-1], ...],
  "FixedPositions": [[0,0,0], [6,0,0]],
  "Q": [1.0, 1.0, ...],
  "LowerBounds": [0.1, ...],
  "UpperBounds": [100, ...],
  "Objectives": [...]
}
```

#### P/Invoke path (new)

Build the same data as flat arrays and pass directly:

```csharp
// Build COO incidence from your edge list
var cooRows = new List<int>();
var cooCols = new List<int>();
var cooVals = new List<double>();

for (int e = 0; e < edges.Count; e++)
{
    var (start, end) = edges[e];
    cooRows.Add(e); cooCols.Add(start); cooVals.Add(-1.0);
    cooRows.Add(e); cooCols.Add(end);   cooVals.Add(+1.0);
}

// Flatten loads: freeNodes.Count × 3
double[] loads = new double[freeNodes.Count * 3];
for (int i = 0; i < freeNodes.Count; i++)
{
    loads[i * 3 + 0] = nodeLoads[freeNodes[i]].X;
    loads[i * 3 + 1] = nodeLoads[freeNodes[i]].Y;
    loads[i * 3 + 2] = nodeLoads[freeNodes[i]].Z;
}

// Similarly flatten fixed positions
double[] fixedPos = new double[fixedNodes.Count * 3];
for (int i = 0; i < fixedNodes.Count; i++)
{
    fixedPos[i * 3 + 0] = anchorPositions[fixedNodes[i]].X;
    fixedPos[i * 3 + 1] = anchorPositions[fixedNodes[i]].Y;
    fixedPos[i * 3 + 2] = anchorPositions[fixedNodes[i]].Z;
}

using var solver = TheseusSolver.Create(
    numEdges, numNodes, freeNodes.Count,
    cooRows.ToArray(), cooCols.ToArray(), cooVals.ToArray(),
    freeNodes.ToArray(), fixedNodes.ToArray(),
    loads, fixedPos,
    qInit, lowerBounds, upperBounds);
```

### 4. Replace objective setup

#### Julia path

```json
{ "type": "TargetXYZ", "weight": 1.0, "indices": [1,2,3], "target": [...] }
```

#### P/Invoke path

```csharp
solver.AddTargetXyz(1.0, nodeIndices, targetXyzFlat);
solver.AddLengthVariation(0.5, edgeIndices, sharpness: 20.0);
solver.AddMinLength(1.0, edgeIndices, thresholds, sharpness: 10.0);
// ... etc
```

### 5. Replace solve + result handling

#### Julia path

The WebSocket response contains:

```json
{
  "Finished": true,
  "Iter": 42,
  "Loss": 0.001,
  "Q": [...],
  "X": [...], "Y": [...], "Z": [...],
  "Losstrace": [...]
}
```

#### P/Invoke path

```csharp
solver.SetSolverOptions(maxIterations: 500);
var result = solver.Optimize();

// result.Xyz is a flat row-major array: xyz[node * 3 + dim]
// Convert to Rhino Point3d:
var points = new Point3d[numNodes];
for (int i = 0; i < numNodes; i++)
{
    points[i] = new Point3d(
        result.Xyz[i * 3 + 0],
        result.Xyz[i * 3 + 1],
        result.Xyz[i * 3 + 2]);
}

// Other result fields:
// result.MemberLengths   — double[numEdges]
// result.MemberForces    — double[numEdges]
// result.ForceDensities  — double[numEdges] (the optimised q)
// result.Reactions       — double[numNodes * 3] (row-major)
// result.Iterations      — int
// result.Converged       — bool
```

### 6. Constrained optimization

If using cable inextensibility constraints (MaxLength):

```csharp
solver.AddConstraintMaxLength(edgeIndices, maxLengths);
var result = solver.OptimizeConstrained(
    muInit: 10.0,
    muFactor: 5.0,
    muMax: 1e6,
    maxOuterIters: 15,
    constraintTol: 1e-3);

// result.ConstraintMaxViolation gives final max g+(k)
```

## Objective Reference

| Method | Rust type | Description |
|--------|-----------|-------------|
| `AddTargetXyz` | TargetXYZ | Minimize 3D distance to target positions |
| `AddTargetXy` | TargetXY | Minimize 2D (XY) distance to target positions |
| `AddTargetPlane` | TargetPlane | Minimize (u,v) distance to target positions on a plane |
| `AddPlanarConstraintAlongDirection` | PlanarConstraintAlongDirection | Planar constraint: pull nodes onto a plane along a direction; no target positions |
| `AddTargetLength` | TargetLength | Minimize difference from target edge lengths |
| `AddLengthVariation` | LengthVariation | Minimize range of edge lengths |
| `AddForceVariation` | ForceVariation | Minimize range of member forces |
| `AddSumForceLength` | SumForceLength | Minimize sum of force x length products |
| `AddMinLength` | MinLength | Barrier penalty for edges below threshold |
| `AddMaxLength` | MaxLength | Barrier penalty for edges above threshold |
| `AddMinForce` | MinForce | Barrier penalty for forces below threshold |
| `AddMaxForce` | MaxForce | Barrier penalty for forces above threshold |
| `AddRigidSetCompare` | RigidSetCompare | Compare pairwise distances between node sets |
| `AddReactionDirection` | ReactionDirection | Align anchor reaction directions |
| `AddReactionDirectionMagnitude` | ReactionDirectionMagnitude | Align reaction directions and magnitudes |

## Important Notes

- **Thread safety**: Each `TheseusSolver` instance is independent. Do not share
  a handle across threads. Create one per solve.
- **Memory**: `TheseusSolver` implements `IDisposable`. Always use `using` or
  call `Dispose()` explicitly to free the native handle.
- **usize/nuint**: Rust `usize` is 8 bytes on 64-bit Windows. The wrapper
  handles the `int` → `nuint` conversion internally.
- **Row-major layout**: All 2D arrays (xyz, loads, targets) are flattened
  row-major: `array[row * 3 + col]`.
