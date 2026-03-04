using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    public class FeaSolveComponent : GH_Component
    {
        public FeaSolveComponent()
            : base("FEA Solve", "FEA Solve",
                "Solve a linear elastic FEA problem (forward or with optimization)",
                "Theseus-FEA", "Design")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Network", "N", "FEA network", GH_ParamAccess.item);
            pManager.AddVectorParameter("Loads", "L", "Point loads (force vectors)", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Load Nodes", "LN", "Node indices for loads (optional)", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Self Weight", "SW", "Include self-weight", GH_ParamAccess.item, false);
            pManager.AddGenericParameter("Opt Config", "Opt", "Optimization config (optional)", GH_ParamAccess.item);

            pManager[2].Optional = true;
            pManager[4].Optional = true;
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Network", "N", "Solved FEA network", GH_ParamAccess.item);
            pManager.AddNumberParameter("Displacements", "U", "Nodal displacements (flat: x0,y0,z0,...)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Reactions", "R", "Reaction forces (flat: x0,y0,z0,...)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Axial Forces", "F", "Axial force per element", GH_ParamAccess.list);
            pManager.AddNumberParameter("Stresses", "σ", "Axial stress per element", GH_ParamAccess.list);
            pManager.AddNumberParameter("Strains", "ε", "Axial strain per element", GH_ParamAccess.list);
            pManager.AddPointParameter("Deformed Nodes", "DN", "Deformed node positions", GH_ParamAccess.list);
            pManager.AddCurveParameter("Deformed Edges", "DE", "Deformed edge curves", GH_ParamAccess.list);
            pManager.AddNumberParameter("Utilization", "Util", "Stress utilization per element", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Iterations", "Iter", "Optimization iterations", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Converged", "Conv", "Whether optimization converged", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            FEA_Network network = null;
            var loads = new List<Vector3d>();
            var loadNodes = new List<int>();
            bool selfWeight = false;
            FeaOptimizationConfig optConfig = null;

            if (!DA.GetData(0, ref network)) return;
            if (!DA.GetDataList(1, loads)) return;
            DA.GetDataList(2, loadNodes);
            DA.GetData(3, ref selfWeight);
            DA.GetData(4, ref optConfig);

            if (network == null || !network.Valid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid FEA network.");
                return;
            }

            if (loads.Count == 0 && !selfWeight)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No loads applied.");
                return;
            }

            try
            {
                FeaSolveResult result;
                List<int>? loadNodeIdx = loadNodes.Count > 0 ? loadNodes : null;

                if (optConfig != null && optConfig.Run)
                {
                    result = FeaSolverService.Solve(network, loads, loadNodeIdx, selfWeight, optConfig);
                }
                else
                {
                    result = FeaSolverService.SolveForward(network, loads, loadNodeIdx, selfWeight);
                }

                DA.SetData(0, result.Network);
                DA.SetDataList(1, result.Displacements);
                DA.SetDataList(2, result.Reactions);
                DA.SetDataList(3, result.AxialForces);
                DA.SetDataList(4, result.Stresses);
                DA.SetDataList(5, result.Strains);
                DA.SetDataList(6, result.DeformedNodes);
                DA.SetDataList(7, result.DeformedEdges);
                DA.SetDataList(8, result.Utilization);
                DA.SetData(9, result.Iterations);
                DA.SetData(10, result.Converged);
            }
            catch (Exception ex)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"FEA solve failed: {ex.Message}");
            }
        }

        protected override Bitmap Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-000000000005");
    }
}
