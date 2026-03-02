using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;
using Rhino.Geometry;
using Ariadne.FDM;

namespace Ariadne.Solver.Components.Experimental;

/// <summary>
/// Experimental: finds force densities via pseudoinverse of the FDM equilibrium
/// system, then performs a forward solve to produce equilibrium geometry closest
/// to the target shape.  Single-pass — no iterative optimisation.
/// </summary>
public class PseudoinverseComponent : GH_Component
{
    public PseudoinverseComponent()
        : base("Pseudoinverse Solve", "Pinv",
            "Find force densities via pseudoinverse of the equilibrium system, then forward-solve. " +
            "Single-pass alternative to gradient-based optimisation.",
            "Ariadne", "Experimental")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Network", "Network", "FDM Network (topology + anchors)", GH_ParamAccess.item);
        pManager.AddPointParameter("Target Points", "Target", "Desired free-node positions (one per free node, matching order)", GH_ParamAccess.list);
        pManager.AddNumberParameter("Force Densities", "q", "Initial force densities (used for forward-solve bounds only)", GH_ParamAccess.list, 10.0);
        pManager.AddVectorParameter("Loads", "Loads", "Loads on free nodes", GH_ParamAccess.list, new Vector3d(0, 0, -1));
        pManager.AddNumberParameter("Regularization", "λ", "Tikhonov regularization parameter (larger = more stable, smaller = closer fit)", GH_ParamAccess.item, 1e-6);
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Network", "Network", "Solved network with updated geometry", GH_ParamAccess.item);
        pManager.AddPointParameter("Nodes", "Nodes", "Solved node positions", GH_ParamAccess.list);
        pManager.AddCurveParameter("Edges", "Edges", "Solved edge curves", GH_ParamAccess.list);
        pManager.AddNumberParameter("Force Densities", "Q", "Computed force densities", GH_ParamAccess.list);
        pManager.AddNumberParameter("Member Forces", "Forces", "Member forces (Q x Length)", GH_ParamAccess.list);
        pManager.AddNumberParameter("Member Lengths", "Lengths", "Member lengths", GH_ParamAccess.list);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        FDM_Network? network = null;
        List<Point3d> targetPoints = [];
        List<double> q = [];
        List<Vector3d> loads = [];
        double regularization = 1e-6;

        if (!DA.GetData(0, ref network)) return;
        if (!DA.GetDataList(1, targetPoints)) return;
        DA.GetDataList(2, q);
        DA.GetDataList(3, loads);
        DA.GetData(4, ref regularization);

        if (network == null || !network.Valid)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid or null network.");
            return;
        }

        int numFree = network.Free.Count;
        if (targetPoints.Count != numFree)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                $"Target points count ({targetPoints.Count}) must match free node count ({numFree}).");
            return;
        }

        // Pack target free-node positions (row-major, nn_free × 3)
        double[] targetFreeXyz = new double[numFree * 3];
        for (int i = 0; i < numFree; i++)
        {
            targetFreeXyz[i * 3 + 0] = targetPoints[i].X;
            targetFreeXyz[i * 3 + 1] = targetPoints[i].Y;
            targetFreeXyz[i * 3 + 2] = targetPoints[i].Z;
        }

        var inputs = new SolverInputs
        {
            QInit = q,
            Loads = loads,
        };

        try
        {
            var result = TheseusSolverService.SolvePseudoinverse(
                network, inputs, targetFreeXyz, regularization);

            DA.SetData(0, result.Network);
            DA.SetDataList(1, result.NodePositions);
            DA.SetDataList(2, result.EdgeCurves);
            DA.SetDataList(3, result.ForceDensities);
            DA.SetDataList(4, result.MemberForces);
            DA.SetDataList(5, result.MemberLengths);
        }
        catch (Exception ex)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, ex.Message);
        }
    }

    protected override Bitmap Icon => Properties.Resources.parameters;

    public override Guid ComponentGuid => new("E1F2A3B4-C5D6-7890-E1F2-A3B4C5D60001");
}
