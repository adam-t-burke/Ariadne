namespace Ariadne.FEA;

using System;
using System.Collections.Generic;
using System.Linq;
using Ariadne.Graphs;
using Rhino.Geometry;
using Theseus.Interop;

/// <summary>
/// High-level API for FEA solving. Marshals FEA_Network + loads into flat arrays
/// for the Rust FFI, creates the solver, applies objectives, and returns results.
/// </summary>
public static class FeaSolverService
{
    internal static BarSolveData SolveForward(
        FEA_Network network,
        List<Vector3d> loads,
        List<int>? loadNodeIndices,
        bool includeSelfWeight)
    {
        var (solver, context) = BuildSolver(network, loads, loadNodeIndices, includeSelfWeight);
        using (solver)
        {
            var result = solver.SolveForward();
            return BuildResult(network, result);
        }
    }

    internal static BarSolveData Solve(
        FEA_Network network,
        List<Vector3d> loads,
        List<int>? loadNodeIndices,
        bool includeSelfWeight,
        FeaOptimizationConfig config)
    {
        var (solver, context) = BuildSolver(network, loads, loadNodeIndices, includeSelfWeight);
        using (solver)
        {
            foreach (var obj in config.Objectives)
            {
                if (obj.IsValid)
                    obj.ApplyTo(solver, context);
            }

            solver.SetSolverOptions(
                config.MaxIterations, config.AbsTol, config.RelTol,
                config.BarrierWeight, config.BarrierSharpness);

            // Set variable mask
            var nodeIdx = config.NodeIndices?.ToArray()
                ?? Enumerable.Range(0, network.Graph.Nn)
                    .Where(n => !network.SupportNodeIndices.Contains(n))
                    .ToArray();
            var areaIdx = config.AreaIndices?.ToArray()
                ?? Enumerable.Range(0, network.Sections.Count).ToArray();
            var supIdx = config.SupportNodeIndices?.ToArray()
                ?? network.SupportNodeIndices.ToArray();

            solver.SetVariableMask(
                config.OptimizeNodePositions, nodeIdx,
                config.OptimizeAreas, areaIdx,
                config.OptimizeSupportPositions, supIdx);

            var result = solver.Optimize();
            return BuildResult(network, result);
        }
    }

    private static (FeaSolver solver, FeaSolverContext context) BuildSolver(
        FEA_Network network,
        List<Vector3d> loads,
        List<int>? loadNodeIndices,
        bool includeSelfWeight)
    {
        var graph = network.Graph;
        int nn = graph.Nn;
        int ne = graph.Ne;

        // Node positions
        var nodePos = new double[nn * 3];
        for (int i = 0; i < nn; i++)
        {
            var pt = graph.Nodes[i].Value;
            nodePos[i * 3] = pt.X;
            nodePos[i * 3 + 1] = pt.Y;
            nodePos[i * 3 + 2] = pt.Z;
        }

        // Edge connectivity
        var edgeNodes = new int[ne * 2];
        var nodeMap = new Dictionary<Node, int>();
        for (int i = 0; i < nn; i++)
            nodeMap[graph.Nodes[i]] = i;
        var edgeMap = new Dictionary<Edge, int>();
        for (int e = 0; e < ne; e++)
        {
            var edge = graph.Edges[e];
            edgeMap[edge] = e;
            edgeNodes[e * 2] = edge.Start.Index;
            edgeNodes[e * 2 + 1] = edge.End.Index;
        }

        // Materials
        int numMat = network.Materials.Count;
        var matData = new double[numMat * 3];
        for (int i = 0; i < numMat; i++)
        {
            matData[i * 3] = network.Materials[i].E;
            matData[i * 3 + 1] = network.Materials[i].Density;
            matData[i * 3 + 2] = network.Materials[i].YieldStress;
        }

        // Sections
        int numSec = network.Sections.Count;
        var secData = new double[numSec];
        for (int i = 0; i < numSec; i++)
            secData[i] = network.Sections[i].Area;

        // Element properties
        var elemProps = new int[ne * 2];
        for (int e = 0; e < ne; e++)
        {
            elemProps[e * 2] = network.MaterialAssignment[e];
            elemProps[e * 2 + 1] = network.SectionAssignment[e];
        }

        // Supports
        int numSup = network.Supports.Count;
        var supData = new int[numSup * 4];
        for (int s = 0; s < numSup; s++)
        {
            supData[s * 4] = network.SupportToNodeMap[s];
            supData[s * 4 + 1] = network.Supports[s].FixX ? 1 : 0;
            supData[s * 4 + 2] = network.Supports[s].FixY ? 1 : 0;
            supData[s * 4 + 3] = network.Supports[s].FixZ ? 1 : 0;
        }

        // Loads
        int numLoads = loads.Count;
        var loadForces = new double[numLoads * 3];
        var loadNodeIdx = new int[numLoads];

        if (loadNodeIndices != null && loadNodeIndices.Count == numLoads)
        {
            for (int i = 0; i < numLoads; i++)
            {
                loadNodeIdx[i] = loadNodeIndices[i];
                loadForces[i * 3] = loads[i].X;
                loadForces[i * 3 + 1] = loads[i].Y;
                loadForces[i * 3 + 2] = loads[i].Z;
            }
        }
        else
        {
            // Apply loads to first N free nodes
            var freeNodes = Enumerable.Range(0, nn)
                .Where(n => !network.SupportNodeIndices.Contains(n))
                .ToList();
            for (int i = 0; i < numLoads; i++)
            {
                loadNodeIdx[i] = i < freeNodes.Count ? freeNodes[i] : freeNodes[^1];
                loadForces[i * 3] = loads[i].X;
                loadForces[i * 3 + 1] = loads[i].Y;
                loadForces[i * 3 + 2] = loads[i].Z;
            }
        }

        var gravity = new double[] { 0.0, 0.0, -9.81 };

        var solver = FeaSolver.Create(
            nn, ne, nodePos, edgeNodes,
            numMat, matData, numSec, secData, elemProps,
            numSup, supData, numLoads, loadForces, loadNodeIdx,
            includeSelfWeight, gravity);

        var context = new FeaSolverContext
        {
            Network = network,
            NodeIndexMap = nodeMap,
            EdgeIndexMap = edgeMap,
        };

        return (solver, context);
    }

    private static BarSolveData BuildResult(FEA_Network network, FeaSolverResult result)
    {
        var graph = network.Graph;
        int nn = graph.Nn;
        int ne = graph.Ne;

        // Convert flat displacement array to Vector3d per node
        var displacements = new Vector3d[nn];
        for (int i = 0; i < nn; i++)
        {
            displacements[i] = new Vector3d(
                result.Displacements[i * 3],
                result.Displacements[i * 3 + 1],
                result.Displacements[i * 3 + 2]);
        }

        // Extract reactions at support nodes only
        int numSupports = network.Supports.Count;
        var reactions = new Vector3d[numSupports];
        var reactionNodeIndices = new int[numSupports];
        for (int s = 0; s < numSupports; s++)
        {
            int nodeIdx = network.SupportToNodeMap[s];
            reactionNodeIndices[s] = nodeIdx;
            reactions[s] = new Vector3d(
                result.Reactions[nodeIdx * 3],
                result.Reactions[nodeIdx * 3 + 1],
                result.Reactions[nodeIdx * 3 + 2]);
        }

        // Build deformed node positions
        var deformedNodes = new Point3d[nn];
        for (int i = 0; i < nn; i++)
        {
            deformedNodes[i] = new Point3d(
                result.DeformedXyz[i * 3],
                result.DeformedXyz[i * 3 + 1],
                result.DeformedXyz[i * 3 + 2]);
        }

        var deformedEdges = new LineCurve[ne];
        for (int e = 0; e < ne; e++)
        {
            var edge = graph.Edges[e];
            deformedEdges[e] = new LineCurve(deformedNodes[edge.Start.Index], deformedNodes[edge.End.Index]);
        }

        return new BarSolveData
        {
            Displacements = displacements,
            Reactions = reactions,
            ReactionNodeIndices = reactionNodeIndices,
            AxialForces = result.AxialForces,
            Stresses = result.Stresses,
            Strains = result.Strains,
            Utilization = result.Utilization,
            DeformedNodes = deformedNodes,
            DeformedEdges = deformedEdges,
            Iterations = result.Iterations,
            Converged = result.Converged,
        };
    }
}
