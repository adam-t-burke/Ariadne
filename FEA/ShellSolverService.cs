namespace Ariadne.FEA;

using System;
using System.Collections.Generic;
using System.Linq;
using Rhino;
using Rhino.Geometry;
using Theseus.Interop;

/// <summary>
/// High-level API for Shell FEA solving. Marshals shell data into flat arrays
/// for the Rust FFI and returns results in managed records.
/// </summary>
public static class ShellSolverService
{
    internal static ShellSolveData SolveForward(ShellSolverData data, bool includeSelfWeight = false)
    {
        double toMeters = GetLengthScaleToMeters();

        var scaledPositions = new double[data.NodePositions.Length];
        for (int i = 0; i < scaledPositions.Length; i++)
            scaledPositions[i] = data.NodePositions[i] * toMeters;

        var gravity = new double[] { 0.0, 0.0, -9.81 };

        using var solver = ShellSolver.Create(
            data.NumNodes, data.NumElements,
            scaledPositions, data.NodeThicknesses, data.Elements, data.NumNodesPerElement,
            data.NumMaterials, data.Materials,
            data.NumSections, data.Sections,
            data.ElementProps,
            data.NumSupports, data.Supports,
            data.NumLoads, data.LoadForces, data.LoadMoments, data.LoadNodes,
            includeSelfWeight, gravity);

        var result = solver.SolveForward();
        return BuildResult(data.NumNodes, data.NumElements, data.NodePositions, result);
    }

    private static ShellSolveData BuildResult(
        int numNodes, int numElements, double[] originalPositions, ShellSolverResult result)
    {
        double fromMeters = GetLengthScaleFromMeters();

        var displacements = new Vector3d[numNodes];
        var rotations = new Vector3d[numNodes];
        var reactionForces = new Vector3d[numNodes];
        var reactionMoments = new Vector3d[numNodes];
        var deformedNodes = new Point3d[numNodes];
        var s1 = new double[numElements];
        var s2 = new double[numElements];

        for (int i = 0; i < numNodes; i++)
        {
            displacements[i] = new Vector3d(
                result.Displacements[i * 6] * fromMeters,
                result.Displacements[i * 6 + 1] * fromMeters,
                result.Displacements[i * 6 + 2] * fromMeters);

            rotations[i] = new Vector3d(
                result.Displacements[i * 6 + 3],
                result.Displacements[i * 6 + 4],
                result.Displacements[i * 6 + 5]);

            reactionForces[i] = new Vector3d(
                result.Reactions[i * 6],
                result.Reactions[i * 6 + 1],
                result.Reactions[i * 6 + 2]);

            reactionMoments[i] = new Vector3d(
                result.Reactions[i * 6 + 3],
                result.Reactions[i * 6 + 4],
                result.Reactions[i * 6 + 5]);

            deformedNodes[i] = new Point3d(
                originalPositions[i * 3] + displacements[i].X,
                originalPositions[i * 3 + 1] + displacements[i].Y,
                originalPositions[i * 3 + 2] + displacements[i].Z);
        }

        for (int e = 0; e < numElements; e++)
        {
            s1[e] = result.PrincipalStresses[e * 2];
            s2[e] = result.PrincipalStresses[e * 2 + 1];
        }

        return new ShellSolveData
        {
            Displacements = displacements,
            Rotations = rotations,
            ReactionForces = reactionForces,
            ReactionMoments = reactionMoments,
            DeformedNodes = deformedNodes,
            S1 = s1,
            S2 = s2,
            TopStresses = result.TopStresses,
            BottomStresses = result.BottomStresses,
            VonMises = result.VonMises,
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
