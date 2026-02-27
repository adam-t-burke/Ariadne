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
    /// <param name="network">FDM network to solve.</param>
    /// <param name="inputs">Solver inputs (bounds, objectives, etc.).</param>
    /// <param name="options">Solver options (optional).</param>
    /// <param name="progressCallback">
    /// Optional callback invoked every ReportFrequency evaluations with
    /// (iteration, loss, xyz[numNodes*3], q[numEdges]).
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
            data.QInit, data.LowerBounds, data.UpperBounds);

        foreach (var objective in inputs.Objectives)
        {
            objective.ApplyTo(solver, context);
        }

        solver.SetSolverOptions(
            maxIterations: options.MaxIterations,
            absTol: options.AbsTol,
            relTol: options.RelTol);

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
            data.QInit, data.LowerBounds, data.UpperBounds);

        var result = solver.SolveForward();

        return BuildResult(network, result, context);
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
        for (int i = 0; i < numFree; i++)
        {
            var load = i < inputs.Loads.Count ? inputs.Loads[i] : inputs.Loads[^1];
            loads[i * 3 + 0] = load.X;
            loads[i * 3 + 1] = load.Y;
            loads[i * 3 + 2] = load.Z;
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

        return new SolverData(
            numEdges, numNodes, numFree,
            cooRows, cooCols, cooVals,
            [.. network.FreeNodes], [.. network.FixedNodes],
            loads, fixedPos,
            Expand(inputs.QInit, numEdges),
            lowerBounds,
            upperBounds);
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
