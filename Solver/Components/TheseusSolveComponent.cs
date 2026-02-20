using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;
using Rhino.Geometry;
using Ariadne.FDM;
using Theseus.Interop;

namespace Ariadne.Solver;

/// <summary>
/// Grasshopper component that solves FDM networks using the Theseus optimization engine.
/// Replaces the WebSocket-based Julia solver with direct P/Invoke calls.
/// </summary>
public class TheseusSolveComponent : GH_Component
{
    private readonly TheseusSolverService _solver = new();
    private SolveResult? _cachedResult;

    public TheseusSolveComponent()
        : base("Theseus Solve", "Theseus",
            "Solve FDM network using Theseus optimization engine.",
            "Ariadne", "Design")
    { }

    protected override void RegisterInputParams(GH_InputParamManager pManager)
    {
        pManager.AddGenericParameter("Network", "Network", "FDM Network to solve", GH_ParamAccess.item);
        pManager.AddNumberParameter("Force Densities", "q", "Initial force densities", GH_ParamAccess.list, 10.0);
        pManager.AddVectorParameter("Loads", "Loads", "Loads on free nodes", GH_ParamAccess.list, new Vector3d(0, 0, -1));
        pManager.AddGenericParameter("Objectives", "OBJ", "Objective functions to minimize", GH_ParamAccess.list);
        pManager.AddNumberParameter("Lower Bounds", "qMin", "Lower bounds on force densities", GH_ParamAccess.list, 0.1);
        pManager.AddNumberParameter("Upper Bounds", "qMax", "Upper bounds on force densities", GH_ParamAccess.list, 100.0);
        pManager.AddIntegerParameter("Max Iterations", "MaxIter", "Maximum solver iterations", GH_ParamAccess.item, 500);
        pManager.AddBooleanParameter("Optimize", "Opt", "Run optimization (when objectives are present)", GH_ParamAccess.item, false);

        pManager[3].Optional = true;
    }

    protected override void RegisterOutputParams(GH_OutputParamManager pManager)
    {
        pManager.AddGenericParameter("Network", "Network", "Solved network with updated geometry", GH_ParamAccess.item);
        pManager.AddPointParameter("Nodes", "Nodes", "Solved node positions", GH_ParamAccess.list);
        pManager.AddCurveParameter("Edges", "Edges", "Solved edge curves", GH_ParamAccess.list);
        pManager.AddNumberParameter("Force Densities", "Q", "Optimized force densities", GH_ParamAccess.list);
        pManager.AddNumberParameter("Member Forces", "Forces", "Member forces (Q x Length)", GH_ParamAccess.list);
        pManager.AddNumberParameter("Member Lengths", "Lengths", "Member lengths", GH_ParamAccess.list);
        pManager.AddIntegerParameter("Iterations", "Iter", "Number of solver iterations", GH_ParamAccess.item);
        pManager.AddBooleanParameter("Converged", "Conv", "Did the solver converge?", GH_ParamAccess.item);
    }

    protected override void SolveInstance(IGH_DataAccess DA)
    {
        FDM_Network? network = null;
        List<double> q = [];
        List<Vector3d> loads = [];
        List<Objective> objectives = [];
        List<double> lb = [];
        List<double> ub = [];
        int maxIter = 500;
        bool optimize = false;

        if (!DA.GetData(0, ref network)) return;
        DA.GetDataList(1, q);
        DA.GetDataList(2, loads);
        DA.GetDataList(3, objectives);
        DA.GetDataList(4, lb);
        DA.GetDataList(5, ub);
        DA.GetData(6, ref maxIter);
        DA.GetData(7, ref optimize);

        if (network == null || !network.Valid)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid or null network");
            return;
        }

        try
        {
            var inputs = new SolverInputs
            {
                QInit = q,
                Loads = loads,
                LowerBounds = lb,
                UpperBounds = ub,
                Objectives = objectives
            };

            if (optimize && objectives.Count > 0)
            {
                var options = new SolverOptions { MaxIterations = maxIter };
                _cachedResult = _solver.Solve(network, inputs, options);

                OutputResult(DA, _cachedResult);

                if (_cachedResult.Converged)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Remark,
                        $"Converged in {_cachedResult.Iterations} iterations");
                }
                else
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                        $"Did not converge after {_cachedResult.Iterations} iterations");
                }
            }
            else
            {
                _cachedResult = _solver.SolveForward(network, inputs);
                OutputResult(DA, _cachedResult);
            }
        }
        catch (TheseusException ex)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Theseus error: {ex.Message}");
        }
        catch (Exception ex)
        {
            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Error: {ex.Message}");
        }
    }

    private static void OutputResult(IGH_DataAccess DA, SolveResult result)
    {
        DA.SetData(0, result.Network);
        DA.SetDataList(1, result.NodePositions);
        DA.SetDataList(2, result.EdgeCurves);
        DA.SetDataList(3, result.ForceDensities);
        DA.SetDataList(4, result.MemberForces);
        DA.SetDataList(5, result.MemberLengths);
        DA.SetData(6, result.Iterations);
        DA.SetData(7, result.Converged);
    }

    protected override Bitmap Icon => Properties.Resources.Create;

    public override Guid ComponentGuid => new("F8A7B2C1-3D4E-5F60-A1B2-C3D4E5F60718");
}
