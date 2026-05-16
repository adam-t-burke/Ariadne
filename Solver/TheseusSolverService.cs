namespace Ariadne.Solver;

using System;
using System.Collections.Generic;
using System.Linq;
using Ariadne.FDM;
using Ariadne.Graphs;
using Rhino.Geometry;
using Theseus.Interop;

/// <summary>
/// Encapsulates all Theseus solver operations. Provides a clean interface
/// for solving FDM networks without Grasshopper dependencies.
/// </summary>
public static class TheseusSolverService
{
    /// <summary>
    /// Solve an FDM network with optimization.
    /// </summary>
    /// <param name="network">The FDM network to solve.</param>
    /// <param name="inputs">Solver inputs including loads, bounds, and objectives.</param>
    /// <param name="options">Optional solver options (iterations, tolerances, etc.).</param>
    /// <param name="progressCallback">
    /// Optional callback invoked every ReportFrequency accepted L-BFGS iterations with
    /// (majorIteration, loss, xyz[numNodes*3], q[numEdges]).
    /// Return <c>true</c> to continue, <c>false</c> to cancel.
    /// </param>
    public static SolveResult Solve(
        FDM_Network network,
        SolverInputs inputs,
        SolverOptions? options = null,
        Func<int, double, double[], double[], bool>? progressCallback = null)
    {
        ValidateCommon(network, inputs);
        ValidateOptimizationBounds(inputs);
        options ??= new SolverOptions();

        var context = BuildContext(network);
        var data = BuildSolverData(network, inputs, context);

        using var solver = TheseusSolver.Create(
            data.NumEdges, data.NumNodes, data.NumFree,
            data.CooRows, data.CooCols, data.CooVals,
            data.FreeIndices, data.FixedIndices,
            data.Loads, data.FixedPositions,
            data.QInit, data.LowerBounds, data.UpperBounds,
            data.VariableNodeIndices, data.VariableSupportKinds, data.SphereRadii,
            data.RollerEnabled, data.RollerLower, data.RollerUpper, data.RailStart, data.RailEnd);

        foreach (var objective in inputs.Objectives)
        {
            objective.ApplyTo(solver, context);
        }

        ApplyLoadConfig(solver, inputs, context);

        solver.SetSolverOptions(
            maxIterations: options.MaxIterations,
            absTol: options.AbsTol,
            relTol: options.RelTol,
            barrierWeight: options.BarrierWeight,
            barrierSharpness: options.BarrierSharpness);

        if (progressCallback != null)
            solver.SetProgressCallback(progressCallback, options.ReportFrequency);

        var result = solver.Optimize();

        return BuildResult(network, result, context);
    }

    /// <summary>
    /// Solve an FDM network without optimization (forward solve only).
    /// Supplies unconstrained bounds so the Rust solver selects LDL
    /// factorization, which handles mixed-sign q values correctly.
    /// </summary>
    public static SolveResult SolveForward(FDM_Network network, SolverInputs inputs)
    {
        ValidateCommon(network, inputs);

        var context = BuildContext(network);
        var data = BuildSolverData(network, inputs, context);

        using var solver = TheseusSolver.Create(
            data.NumEdges, data.NumNodes, data.NumFree,
            data.CooRows, data.CooCols, data.CooVals,
            data.FreeIndices, data.FixedIndices,
            data.Loads, data.FixedPositions,
            data.QInit, data.LowerBounds, data.UpperBounds,
            data.VariableNodeIndices, data.VariableSupportKinds, data.SphereRadii,
            data.RollerEnabled, data.RollerLower, data.RollerUpper, data.RailStart, data.RailEnd);

        ApplyLoadConfig(solver, inputs, context);

        var result = solver.SolveForward();

        return BuildResult(network, result, context);
    }

    /// <summary>
    /// Solve for force densities via pseudoinverse of the equilibrium system.
    /// Given target free-node positions, finds q that best satisfies equilibrium,
    /// then performs a forward FDM solve with the resulting q.
    /// </summary>
    public static SolveResult SolvePseudoinverse(
        FDM_Network network,
        SolverInputs inputs,
        double[] targetFreeXyz,
        double regularization = 1e-6,
        bool useL2 = true,
        int maxL1Iter = 20,
        bool useAugmented = false,
        bool enforceZeroRx = false,
        bool enforceZeroRy = false,
        bool enforceZeroRz = false,
        bool solveForQ = true)
    {
        ValidateCommon(network, inputs);
        var context = BuildContext(network);
        var data = BuildSolverData(network, inputs, context);

        using var solver = TheseusSolver.Create(
            data.NumEdges, data.NumNodes, data.NumFree,
            data.CooRows, data.CooCols, data.CooVals,
            data.FreeIndices, data.FixedIndices,
            data.Loads, data.FixedPositions,
            data.QInit, data.LowerBounds, data.UpperBounds,
            data.VariableNodeIndices, data.VariableSupportKinds, data.SphereRadii,
            data.RollerEnabled, data.RollerLower, data.RollerUpper, data.RailStart, data.RailEnd);

        ApplyLoadConfig(solver, inputs, context);

        var result = solver.SolvePseudoinverse(targetFreeXyz, regularization, useL2, maxL1Iter, useAugmented,
            enforceZeroRx, enforceZeroRy, enforceZeroRz, solveForQ);
        return BuildResult(network, result, context);
    }

    /// <summary>
    /// Solve for non-negative force densities via NNLS (spectral projected gradient).
    /// Given target free-node positions, finds q >= 0 minimising the equilibrium residual,
    /// then performs a forward FDM solve with the resulting q.
    /// </summary>
    public static SolveResult SolveNnls(
        FDM_Network network,
        SolverInputs inputs,
        double[] targetFreeXyz,
        int maxIter = 500,
        double tol = 1e-6)
    {
        ValidateCommon(network, inputs);
        var context = BuildContext(network);
        var data = BuildSolverData(network, inputs, context);

        using var solver = TheseusSolver.Create(
            data.NumEdges, data.NumNodes, data.NumFree,
            data.CooRows, data.CooCols, data.CooVals,
            data.FreeIndices, data.FixedIndices,
            data.Loads, data.FixedPositions,
            data.QInit, data.LowerBounds, data.UpperBounds,
            data.VariableNodeIndices, data.VariableSupportKinds, data.SphereRadii,
            data.RollerEnabled, data.RollerLower, data.RollerUpper, data.RailStart, data.RailEnd);

        ApplyLoadConfig(solver, inputs, context);

        var result = solver.SolveNnls(targetFreeXyz, maxIter, tol);
        return BuildResult(network, result, context);
    }

    /// <summary>
    /// Resolves a list of load-node positions to their indices in the network's
    /// free-node list (0-based). Matches by proximity using the network's ATol.
    /// </summary>
    /// <returns>List of free-node indices corresponding to each input point.</returns>
    /// <exception cref="ArgumentException">Thrown when a point doesn't match any free node.</exception>
    public static List<int> ResolveLoadNodeIndices(FDM_Network network, IReadOnlyList<Point3d> loadNodePoints)
    {
        var indices = new List<int>(loadNodePoints.Count);
        double tol = network.ATol;

        for (int p = 0; p < loadNodePoints.Count; p++)
        {
            var pt = loadNodePoints[p];
            int bestIdx = -1;
            double bestDist = double.MaxValue;

            for (int i = 0; i < network.Free.Count; i++)
            {
                double dist = pt.DistanceTo(network.Free[i].Value);
                if (dist < bestDist)
                {
                    bestDist = dist;
                    bestIdx = i;
                }
            }

            if (bestIdx < 0 || bestDist > tol)
                throw new ArgumentException(
                    $"Load node at ({pt.X:G4}, {pt.Y:G4}, {pt.Z:G4}) does not match any free node within tolerance {tol}.");

            indices.Add(bestIdx);
        }

        return indices;
    }

    #region Validation

    private static void ValidateCommon(FDM_Network network, SolverInputs inputs)
    {
        if (network == null)
            throw new ArgumentNullException(nameof(network));
        if (!network.Valid)
            throw new ArgumentException("Network is not valid. Check anchor definitions.", nameof(network));
        if (inputs.Loads == null || inputs.Loads.Count == 0)
            throw new ArgumentException("Loads list cannot be empty.", nameof(inputs));
        if (inputs.QInit == null || inputs.QInit.Count == 0)
            throw new ArgumentException("Initial force densities cannot be empty.", nameof(inputs));
    }

    private static void ValidateOptimizationBounds(SolverInputs inputs)
    {
        if (inputs.LowerBounds == null || inputs.LowerBounds.Count == 0)
            throw new ArgumentException("Lower bounds cannot be empty for optimization.", nameof(inputs));
        if (inputs.UpperBounds == null || inputs.UpperBounds.Count == 0)
            throw new ArgumentException("Upper bounds cannot be empty for optimization.", nameof(inputs));
    }

    #endregion

    #region Private Methods

    private static SolverContext BuildContext(FDM_Network network)
    {
        var nodeMap = new Dictionary<Node, int>(network.Graph.Nn);
        for (int i = 0; i < network.Graph.Nn; i++)
            nodeMap[network.Graph.Nodes[i]] = i;

        var edgeMap = new Dictionary<Edge, int>(network.Graph.Ne);
        for (int i = 0; i < network.Graph.Ne; i++)
            edgeMap[network.Graph.Edges[i]] = i;

        return new SolverContext
        {
            Network = network,
            NodeIndexMap = nodeMap,
            EdgeIndexMap = edgeMap
        };
    }

    private static SolverData BuildSolverData(FDM_Network network, SolverInputs inputs, SolverContext context)
    {
        var graph = network.Graph;
        int numEdges = graph.Ne;
        int numNodes = graph.Nn;
        int numFree = network.FreeNodes.Count;
        int numFixed = network.FixedNodes.Count;

        var cooRows = new int[numEdges * 2];
        var cooCols = new int[numEdges * 2];
        var cooVals = new double[numEdges * 2];

        for (int e = 0; e < numEdges; e++)
        {
            var edge = graph.Edges[e];
            int startIdx = context.NodeIndexMap[edge.Start];
            int endIdx = context.NodeIndexMap[edge.End];

            cooRows[e * 2] = e;
            cooCols[e * 2] = startIdx;
            cooVals[e * 2] = -1.0;
            cooRows[e * 2 + 1] = e;
            cooCols[e * 2 + 1] = endIdx;
            cooVals[e * 2 + 1] = 1.0;
        }

        double[] loads = new double[numFree * 3];
        if (inputs.LoadNodeIndices is { Count: > 0 } targetIndices)
        {
            if (inputs.Loads.Count != 1 && inputs.Loads.Count != targetIndices.Count)
                throw new ArgumentException(
                    $"When load nodes are specified ({targetIndices.Count}), " +
                    $"provide either 1 load or exactly {targetIndices.Count} loads (got {inputs.Loads.Count}).");

            for (int i = 0; i < targetIndices.Count; i++)
            {
                int freeIdx = targetIndices[i];
                var load = inputs.Loads.Count == 1 ? inputs.Loads[0] : inputs.Loads[i];
                loads[freeIdx * 3 + 0] = load.X;
                loads[freeIdx * 3 + 1] = load.Y;
                loads[freeIdx * 3 + 2] = load.Z;
            }
        }
        else
        {
            for (int i = 0; i < numFree; i++)
            {
                var load = i < inputs.Loads.Count ? inputs.Loads[i] : inputs.Loads[^1];
                loads[i * 3 + 0] = load.X;
                loads[i * 3 + 1] = load.Y;
                loads[i * 3 + 2] = load.Z;
            }
        }

        double[] fixedPos = new double[numFixed * 3];
        for (int i = 0; i < numFixed; i++)
        {
            var pos = network.Fixed[i].Value;
            fixedPos[i * 3 + 0] = pos.X;
            fixedPos[i * 3 + 1] = pos.Y;
            fixedPos[i * 3 + 2] = pos.Z;
        }

        double[] lowerBounds = inputs.LowerBounds != null
            ? Expand(inputs.LowerBounds, numEdges)
            : ExpandConstant(double.NegativeInfinity, numEdges);

        double[] upperBounds = inputs.UpperBounds != null
            ? Expand(inputs.UpperBounds, numEdges)
            : ExpandConstant(double.PositiveInfinity, numEdges);

        var supportData = BuildVariableSupportData(inputs.VariableSupports, network, context);

        return new SolverData(
            numEdges, numNodes, numFree,
            cooRows, cooCols, cooVals,
            [.. network.FreeNodes], [.. network.FixedNodes],
            loads, fixedPos,
            Expand(inputs.QInit, numEdges),
            lowerBounds,
            upperBounds,
            supportData.VariableNodeIndices,
            supportData.VariableSupportKinds,
            supportData.SphereRadii,
            supportData.RollerEnabled,
            supportData.RollerLower,
            supportData.RollerUpper,
            supportData.RailStart,
            supportData.RailEnd);
    }

    private readonly record struct VariableSupportData(
        int[] VariableNodeIndices,
        int[] VariableSupportKinds,
        double[] SphereRadii,
        byte[] RollerEnabled,
        double[] RollerLower,
        double[] RollerUpper,
        double[] RailStart,
        double[] RailEnd);

    private static VariableSupportData BuildVariableSupportData(
        List<VariableSupportConfig>? supports,
        FDM_Network network,
        SolverContext context)
    {
        if (supports == null || supports.Count == 0)
        {
            return new VariableSupportData(
                [],
                [],
                [],
                [],
                [],
                [],
                [],
                []);
        }

        var nodeIndices = new List<int>();
        var kinds = new List<int>();
        var sphereRadii = new List<double>();
        var rollerEnabled = new List<byte>();
        var rollerLower = new List<double>();
        var rollerUpper = new List<double>();
        var railStart = new List<double>();
        var railEnd = new List<double>();

        var seen = new HashSet<int>();
        foreach (var support in supports)
        {
            if (support.Nodes == null || support.Nodes.Count == 0)
                continue;

            foreach (var node in support.Nodes)
            {
                if (!context.NodeIndexMap.TryGetValue(node, out int nodeIdx))
                    throw new ArgumentException("Variable support node is not part of the current network.");
                if (!network.FixedNodes.Contains(nodeIdx))
                    throw new ArgumentException("Variable supports can only target fixed nodes.");
                if (!seen.Add(nodeIdx))
                    throw new ArgumentException($"Node index {nodeIdx} has multiple variable support definitions.");

                nodeIndices.Add(nodeIdx);

                switch (support)
                {
                    case SphereVariableSupport sphere:
                        if (!double.IsFinite(sphere.Radius) || sphere.Radius <= 0.0)
                            throw new ArgumentException("Sphere variable support radius must be positive and finite.");
                        kinds.Add(0);
                        sphereRadii.Add(sphere.Radius);
                        rollerEnabled.AddRange([0, 0, 0]);
                        rollerLower.AddRange([0.0, 0.0, 0.0]);
                        rollerUpper.AddRange([0.0, 0.0, 0.0]);
                        railStart.AddRange([0.0, 0.0, 0.0]);
                        railEnd.AddRange([0.0, 0.0, 0.0]);
                        break;

                    case RollerVariableSupport roller:
                        if (!(roller.FreeX || roller.FreeY || roller.FreeZ))
                            throw new ArgumentException("Roller variable support must free at least one axis.");
                        if ((roller.FreeX && !roller.DomainX.IsIncreasing) ||
                            (roller.FreeY && !roller.DomainY.IsIncreasing) ||
                            (roller.FreeZ && !roller.DomainZ.IsIncreasing))
                            throw new ArgumentException("Roller enabled-axis domains must be increasing.");

                        kinds.Add(1);
                        sphereRadii.Add(0.0);
                        rollerEnabled.Add((byte)(roller.FreeX ? 1 : 0));
                        rollerEnabled.Add((byte)(roller.FreeY ? 1 : 0));
                        rollerEnabled.Add((byte)(roller.FreeZ ? 1 : 0));
                        rollerLower.Add(roller.DomainX.T0);
                        rollerLower.Add(roller.DomainY.T0);
                        rollerLower.Add(roller.DomainZ.T0);
                        rollerUpper.Add(roller.DomainX.T1);
                        rollerUpper.Add(roller.DomainY.T1);
                        rollerUpper.Add(roller.DomainZ.T1);
                        railStart.AddRange([0.0, 0.0, 0.0]);
                        railEnd.AddRange([0.0, 0.0, 0.0]);
                        break;

                    case RailVariableSupport rail:
                        if (!rail.Rail.IsValid || rail.Rail.Length <= 1e-12)
                            throw new ArgumentException("Rail variable support requires a non-degenerate line segment.");
                        kinds.Add(2);
                        sphereRadii.Add(0.0);
                        rollerEnabled.AddRange([0, 0, 0]);
                        rollerLower.AddRange([0.0, 0.0, 0.0]);
                        rollerUpper.AddRange([0.0, 0.0, 0.0]);
                        railStart.Add(rail.Rail.FromX);
                        railStart.Add(rail.Rail.FromY);
                        railStart.Add(rail.Rail.FromZ);
                        railEnd.Add(rail.Rail.ToX);
                        railEnd.Add(rail.Rail.ToY);
                        railEnd.Add(rail.Rail.ToZ);
                        break;

                    default:
                        throw new ArgumentException($"Unsupported variable support type {support.GetType().Name}.");
                }
            }
        }

        return new VariableSupportData(
            [.. nodeIndices],
            [.. kinds],
            [.. sphereRadii],
            [.. rollerEnabled],
            [.. rollerLower],
            [.. rollerUpper],
            [.. railStart],
            [.. railEnd]);
    }

    private static void ApplyLoadConfig(TheseusSolver solver, SolverInputs inputs, SolverContext context)
    {
        if (inputs.SelfWeight is SelfWeightConfig.Prescribed prescribed)
        {
            int numEdges = context.EdgeIndexMap.Count;
            var mu = new double[numEdges];
            for (int i = 0; i < numEdges; i++)
                mu[i] = i < prescribed.LinearDensities.Count
                    ? prescribed.LinearDensities[i]
                    : prescribed.LinearDensities[^1];

            var g = new[] { prescribed.Gravity.X, prescribed.Gravity.Y, prescribed.Gravity.Z };
            solver.SetSelfWeightPrescribed(mu, g,
                prescribed.MaxIters, prescribed.Tolerance, prescribed.Relaxation);
        }
        else if (inputs.SelfWeight is SelfWeightConfig.Sizing sizing)
        {
            var g = new[] { sizing.Gravity.X, sizing.Gravity.Y, sizing.Gravity.Z };
            solver.SetSelfWeightSizing(sizing.Rho, sizing.Sigma, g,
                sizing.MaxIters, sizing.Tolerance, sizing.Relaxation);
        }

        if (inputs.Pressure is { } pressure)
        {
            int numFaces = pressure.Faces.Count;
            var faces = BuildFaceArray(pressure.Faces, context);

            switch (pressure)
            {
                case PressureConfig.Normal normal:
                {
                    var press = ExpandPerFace(normal.Pressures, numFaces);
                    solver.SetPressure(faces, press,
                        normal.MaxIters, normal.Tolerance, normal.Relaxation);
                    break;
                }
                case PressureConfig.Hydrostatic hydro:
                {
                    var up = new[] { hydro.UpDirection.X, hydro.UpDirection.Y, hydro.UpDirection.Z };
                    solver.SetPressureHydrostatic(faces,
                        hydro.RhoFluid, hydro.GMagnitude, hydro.ZDatum, up,
                        hydro.MaxIters, hydro.Tolerance, hydro.Relaxation);
                    break;
                }
                case PressureConfig.Directional dir:
                {
                    var press = ExpandPerFace(dir.Pressures, numFaces);
                    var d = new[] { dir.Direction.X, dir.Direction.Y, dir.Direction.Z };
                    solver.SetPressureDirectional(faces, press, d,
                        dir.MaxIters, dir.Tolerance, dir.Relaxation);
                    break;
                }
            }
        }
    }

    private static int[][] BuildFaceArray(List<List<int>> faces, SolverContext context)
    {
        int numFaces = faces.Count;
        var result = new int[numFaces][];
        for (int f = 0; f < numFaces; f++)
        {
            var faceVerts = faces[f];
            result[f] = new int[faceVerts.Count];
            for (int v = 0; v < faceVerts.Count; v++)
            {
                var node = context.Network.Graph.Nodes[faceVerts[v]];
                result[f][v] = context.NodeIndexMap[node];
            }
        }
        return result;
    }

    private static double[] ExpandPerFace(List<double> values, int numFaces)
    {
        var result = new double[numFaces];
        for (int f = 0; f < numFaces; f++)
            result[f] = f < values.Count ? values[f] : values[^1];
        return result;
    }

    private static SolveResult BuildResult(FDM_Network oldNetwork, SolverResult result, SolverContext context)
    {
        var newNetwork = BuildSolvedNetwork(oldNetwork, result, context);

        return new SolveResult
        {
            Network = newNetwork,
            ForceDensities = result.ForceDensities,
            MemberForces = result.MemberForces,
            MemberLengths = result.MemberLengths,
            Reactions = result.Reactions,
            Iterations = result.Iterations,
            Converged = result.Converged,
            TerminationReason = result.TerminationReason
        };
    }

    private static FDM_Network BuildSolvedNetwork(FDM_Network oldNetwork, SolverResult result, SolverContext context)
    {
        var oldGraph = oldNetwork.Graph;
        int numNodes = oldGraph.Nn;
        int numEdges = oldGraph.Ne;

        var newNodes = new List<Node>(numNodes);
        for (int i = 0; i < numNodes; i++)
        {
            newNodes.Add(new Node
            {
                Index = i,
                Anchor = oldGraph.Nodes[i].Anchor,
                Value = new Point3d(
                    result.Xyz[i * 3 + 0],
                    result.Xyz[i * 3 + 1],
                    result.Xyz[i * 3 + 2])
            });
        }

        for (int i = 0; i < numNodes; i++)
        {
            newNodes[i].Neighbors = oldGraph.Nodes[i].Neighbors
                .Select(n => newNodes[context.NodeIndexMap[n]])
                .ToList();
        }

        var newEdges = new List<Edge>(numEdges);
        for (int i = 0; i < numEdges; i++)
        {
            int startIdx = context.NodeIndexMap[oldGraph.Edges[i].Start];
            int endIdx = context.NodeIndexMap[oldGraph.Edges[i].End];

            var edge = new Edge
            {
                Start = newNodes[startIdx],
                End = newNodes[endIdx],
                Q = result.ForceDensities[i],
                ReferenceID = oldGraph.Edges[i].ReferenceID
            };
            edge.Value = new LineCurve(edge.Start.Value, edge.End.Value);
            newEdges.Add(edge);
        }

        var newGraph = new Graph
        {
            Tolerance = oldGraph.Tolerance,
            Nodes = newNodes,
            Edges = newEdges,
            EdgeInputMap = oldGraph.EdgeInputMap
        };

        newGraph.EdgeIndicesToTree();
        newGraph.AdjacencyListToTree();
        newGraph.BuildOutputEdgeTree();

        return new FDM_Network
        {
            Graph = newGraph,
            Valid = true,
            FreeNodes = oldNetwork.FreeNodes,
            FixedNodes = oldNetwork.FixedNodes,
            ATol = oldNetwork.ATol,
            ETol = oldNetwork.ETol,
            Free = oldNetwork.FreeNodes.Select(i => newNodes[i]).ToList(),
            Fixed = oldNetwork.FixedNodes.Select(i => newNodes[i]).ToList(),
            Anchors = oldNetwork.FixedNodes.Select(i => newNodes[i].Value).ToList()
        };
    }

    private static double[] Expand(IReadOnlyList<double> values, int length)
    {
        if (values.Count == 0)
            throw new ArgumentException("Values list cannot be empty", nameof(values));

        var result = new double[length];
        for (int i = 0; i < length; i++)
            result[i] = i < values.Count ? values[i] : values[^1];
        return result;
    }

    private static double[] ExpandConstant(double value, int length)
    {
        var result = new double[length];
        Array.Fill(result, value);
        return result;
    }

    #endregion
}
