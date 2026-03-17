namespace Ariadne.FEA;

using System;
using System.Collections.Generic;
using System.Linq;
using Rhino;
using Rhino.Geometry;
using Theseus.Interop;

/// <summary>
/// High-level API for solid FEA solving. Marshals FEA_SolidNetwork + loads into flat arrays
/// for the Rust FFI, creates the solver, and returns results.
/// </summary>
public static class SolidSolverService
{
    internal static SolidSolveData SolveForward(
        FEA_SolidNetwork network,
        List<Vector3d>? loads = null,
        List<int>? loadNodeIndices = null,
        bool includeSelfWeight = false)
    {
        using var solver = BuildSolver(network, loads, loadNodeIndices, includeSelfWeight);
        var result = solver.SolveForward();
        return BuildResult(network, result);
    }

    private static SolidSolver BuildSolver(
        FEA_SolidNetwork network,
        List<Vector3d>? loads,
        List<int>? loadNodeIndices,
        bool includeSelfWeight)
    {
        var mesh = network.Mesh;
        int nn = mesh.Nn;
        int ne = mesh.Ne;
        double toMeters = GetLengthScaleToMeters();

        var nodePos = new double[nn * 3];
        for (int i = 0; i < nn; i++)
        {
            var pt = mesh.TetNodes[i];
            nodePos[i * 3] = pt.X * toMeters;
            nodePos[i * 3 + 1] = pt.Y * toMeters;
            nodePos[i * 3 + 2] = pt.Z * toMeters;
        }

        var elements = new int[ne * 4];
        for (int e = 0; e < ne; e++)
        {
            var elem = mesh.Elements[e];
            elements[e * 4] = elem[0];
            elements[e * 4 + 1] = elem[1];
            elements[e * 4 + 2] = elem[2];
            elements[e * 4 + 3] = elem[3];
        }

        int numMat = network.Materials.Count;
        var matData = new double[numMat * 4];
        for (int i = 0; i < numMat; i++)
        {
            matData[i * 4] = network.Materials[i].E;
            matData[i * 4 + 1] = network.Materials[i].Nu;
            matData[i * 4 + 2] = network.Materials[i].Density;
            matData[i * 4 + 3] = network.Materials[i].YieldStress;
        }

        var elemProps = new int[ne];
        for (int e = 0; e < ne; e++)
            elemProps[e] = network.MaterialAssignment[e];

        // Supports: use FeaSupport with matched node indices
        int numSupports = network.Supports.Count;
        var supData = new int[numSupports * 4];
        for (int s = 0; s < numSupports; s++)
        {
            supData[s * 4] = network.SupportToNodeMap[s];
            supData[s * 4 + 1] = network.Supports[s].FixX ? 1 : 0;
            supData[s * 4 + 2] = network.Supports[s].FixY ? 1 : 0;
            supData[s * 4 + 3] = network.Supports[s].FixZ ? 1 : 0;
        }

        // Loads: passed as separate vectors + node indices
        int numLoads = loads?.Count ?? 0;
        var loadForces = new double[numLoads * 3];
        var loadNodes = new int[numLoads];

        if (loads != null && numLoads > 0)
        {
            var freeNodes = Enumerable.Range(0, nn)
                .Where(n => !network.SupportNodeIndices.Contains(n))
                .ToList();

            for (int i = 0; i < numLoads; i++)
            {
                loadForces[i * 3] = loads[i].X;
                loadForces[i * 3 + 1] = loads[i].Y;
                loadForces[i * 3 + 2] = loads[i].Z;

                if (loadNodeIndices != null && i < loadNodeIndices.Count)
                    loadNodes[i] = loadNodeIndices[i];
                else if (i < freeNodes.Count)
                    loadNodes[i] = freeNodes[i];
                else
                    loadNodes[i] = freeNodes.Count > 0 ? freeNodes[^1] : 0;
            }
        }

        var gravity = new double[] { 0.0, 0.0, -9.81 };

        return SolidSolver.Create(
            nn, ne, nodePos, elements,
            numMat, matData, elemProps,
            numSupports, supData,
            numLoads, loadForces, loadNodes,
            includeSelfWeight, gravity);
    }

    private static SolidSolveData BuildResult(FEA_SolidNetwork network, SolidSolverResult result)
    {
        var mesh = network.Mesh;
        int nn = mesh.Nn;
        int ne = mesh.Ne;
        int ngp = SolidConstants.NumGP;
        double fromMeters = GetLengthScaleFromMeters();

        var displacements = new Vector3d[nn];
        for (int i = 0; i < nn; i++)
        {
            displacements[i] = new Vector3d(
                result.Displacements[i * 3] * fromMeters,
                result.Displacements[i * 3 + 1] * fromMeters,
                result.Displacements[i * 3 + 2] * fromMeters);
        }

        int numSup = network.Supports.Count;
        var reactions = new Vector3d[numSup];
        var reactionNodeIndices = new int[numSup];
        for (int s = 0; s < numSup; s++)
        {
            int ni = network.SupportToNodeMap[s];
            reactionNodeIndices[s] = ni;
            reactions[s] = new Vector3d(
                result.Reactions[ni * 3],
                result.Reactions[ni * 3 + 1],
                result.Reactions[ni * 3 + 2]);
        }

        // Average GP stresses to one per element for downstream consumers
        var stresses = new double[ne, 6];
        double inv = 1.0 / ngp;
        for (int e = 0; e < ne; e++)
            for (int c = 0; c < 6; c++)
            {
                double sum = 0;
                for (int gp = 0; gp < ngp; gp++)
                    sum += result.Stresses[e * ngp * 6 + gp * 6 + c];
                stresses[e, c] = sum * inv;
            }

        // Average GP von Mises to one per element for visualization
        var vonMises = new double[ne];
        for (int e = 0; e < ne; e++)
        {
            double sum = 0;
            for (int gp = 0; gp < ngp; gp++)
                sum += result.VonMises[e * ngp + gp];
            vonMises[e] = sum * inv;
        }

        var deformedNodes = new Point3d[nn];
        for (int i = 0; i < nn; i++)
        {
            deformedNodes[i] = new Point3d(
                result.DeformedXyz[i * 3] * fromMeters,
                result.DeformedXyz[i * 3 + 1] * fromMeters,
                result.DeformedXyz[i * 3 + 2] * fromMeters);
        }

        return new SolidSolveData
        {
            Displacements = displacements,
            Reactions = reactions,
            ReactionNodeIndices = reactionNodeIndices,
            Stresses = stresses,
            VonMises = vonMises,
            DeformedNodes = deformedNodes,
        };
    }

    private static double GetLengthScaleToMeters()
    {
        var doc = RhinoDoc.ActiveDoc;
        return doc == null ? 1.0 : RhinoMath.UnitScale(doc.ModelUnitSystem, UnitSystem.Meters);
    }

    private static double GetLengthScaleFromMeters()
    {
        var doc = RhinoDoc.ActiveDoc;
        return doc == null ? 1.0 : RhinoMath.UnitScale(UnitSystem.Meters, doc.ModelUnitSystem);
    }
}

internal sealed record SolidSolveData
{
    public required Vector3d[] Displacements { get; init; }
    public required Vector3d[] Reactions { get; init; }
    public required int[] ReactionNodeIndices { get; init; }
    public required double[,] Stresses { get; init; }
    public required double[] VonMises { get; init; }
    public required Point3d[] DeformedNodes { get; init; }
}
