using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Ariadne.Graphs;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    public class FeaSolveComponent : GH_Component
    {
        public FeaSolveComponent()
            : base("FEA Solve", "FEASolve",
                "Solve a linear elastic FEA problem for bar, solid, or coupled models",
                "Theseus-FEA", "Solver")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Model", "M", "FEA model (bar, solid, or coupled)", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Self Weight", "SW", "Include self-weight", GH_ParamAccess.item, false);
            pManager.AddGenericParameter("Opt Config", "Opt", "Optimization config (optional, bar only)", GH_ParamAccess.item);

            pManager[2].Optional = true;
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Results", "Res", "Full unified FEA solve result (contains DeformedModel)", GH_ParamAccess.item);
            pManager.AddVectorParameter("Displacements", "U", "Displacement vector per node", GH_ParamAccess.list);
            pManager.AddVectorParameter("Reactions", "R", "Reaction force vectors at supports", GH_ParamAccess.list);
            pManager.AddGenericParameter("Reaction Nodes", "RN", "Node objects corresponding to each reaction", GH_ParamAccess.list);
            pManager.AddNumberParameter("Stresses", "σ", "Axial stress (bar) or von Mises stress (solid) per element", GH_ParamAccess.list);
            pManager.AddNumberParameter("Axial Forces", "F", "Axial force per bar element (empty for solid)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Strains", "ε", "Axial strain per bar element (empty for solid)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Utilization", "Util", "Utilization per element", GH_ParamAccess.list);
            pManager.AddPointParameter("Deformed Nodes", "DN", "Deformed node positions", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Iterations", "Iter", "Solver iterations", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Converged", "Conv", "Whether the solver converged", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            FEA_Model? model = null;
            bool selfWeight = false;
            FeaOptimizationConfig? optConfig = null;

            if (!DA.GetData(0, ref model)) return;
            DA.GetData(1, ref selfWeight);
            DA.GetData(2, ref optConfig);

            if (model == null || !model.Valid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                    model?.ValidationMessage ?? "Invalid FEA model.");
                return;
            }

            var loads = model.Loads.Select(load => load.Force).ToList();
            var loadNodeIdx = model.Loads.Select(load => load.NodeIndex).ToList();

            int constrainedLoadCount = 0;
            foreach (var load in model.Loads)
            {
                var constrained = model.GetConstrainedDofs(load.NodeIndex);
                if ((constrained[0] && Math.Abs(load.Force.X) > 1e-12) ||
                    (constrained[1] && Math.Abs(load.Force.Y) > 1e-12) ||
                    (constrained[2] && Math.Abs(load.Force.Z) > 1e-12))
                {
                    constrainedLoadCount++;
                }
            }

            if (constrainedLoadCount > 0)
            {
                AddRuntimeMessage(
                    GH_RuntimeMessageLevel.Warning,
                    $"{constrainedLoadCount} load(s) act on constrained DOFs; those components will appear as reactions rather than free displacement.");
            }

            if (loads.Count == 0 && !selfWeight)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No loads applied.");
                return;
            }

            if (model.IsCoupled)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                    "Coupled bar+solid models are not yet supported.");
                return;
            }

            try
            {
                FeaResult result;

                if (model.HasBarElements)
                {
                    result = SolveBar(model, loads, loadNodeIdx, selfWeight, optConfig);
                }
                else
                {
                    result = SolveSolid(model, loads, loadNodeIdx, selfWeight);
                }

                DA.SetData(0, result);
                DA.SetDataList(1, result.Displacements);
                DA.SetDataList(2, result.Reactions);

                if (model.HasBarElements)
                {
                    var reactionNodes = new List<Node>();
                    foreach (int ni in result.ReactionNodeIndices)
                        reactionNodes.Add(model.BarGraph!.Nodes[ni]);
                    DA.SetDataList(3, reactionNodes);
                    DA.SetDataList(4, result.BarStresses ?? []);
                    DA.SetDataList(5, result.AxialForces ?? []);
                    DA.SetDataList(6, result.BarStrains ?? []);
                }
                else
                {
                    var reactionNodes = new List<Node>();
                    var mesh = model.SolidMesh!;
                    foreach (int ni in result.ReactionNodeIndices)
                    {
                        reactionNodes.Add(new Node
                        {
                            Value = mesh.TetNodes[ni],
                            Index = ni
                        });
                    }
                    DA.SetDataList(3, reactionNodes);
                    DA.SetDataList(4, result.VonMises ?? []);
                    DA.SetDataList(5, new List<double>());
                    DA.SetDataList(6, new List<double>());
                }

                DA.SetDataList(7, result.Utilization);
                DA.SetDataList(8, result.DeformedNodes);
                DA.SetData(9, result.Iterations);
                DA.SetData(10, result.Converged);
            }
            catch (Exception ex)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"FEA solve failed: {ex.Message}");
            }
        }

        private static FeaResult SolveBar(
            FEA_Model model,
            List<Vector3d> loads,
            List<int>? loadNodeIdx,
            bool selfWeight,
            FeaOptimizationConfig? optConfig)
        {
            var network = model.ToBarNetwork();

            BarSolveData barResult;
            if (optConfig != null && optConfig.Run)
                barResult = FeaSolverService.Solve(network, loads, loadNodeIdx, selfWeight, optConfig);
            else
                barResult = FeaSolverService.SolveForward(network, loads, loadNodeIdx, selfWeight);

            var deformedModel = new FEA_Model(model);
            for (int i = 0; i < barResult.DeformedNodes.Length; i++)
                deformedModel.BarGraph!.Nodes[i].Value = barResult.DeformedNodes[i];
            for (int e = 0; e < barResult.DeformedEdges.Length; e++)
                deformedModel.BarGraph!.Edges[e].Value = barResult.DeformedEdges[e];

            int ne = model.BarGraph!.Ne;
            var utilization = barResult.Utilization ?? new double[ne];

            return new FeaResult
            {
                Model = model,
                DeformedModel = deformedModel,
                Displacements = barResult.Displacements,
                Reactions = barResult.Reactions,
                ReactionNodeIndices = barResult.ReactionNodeIndices,
                DeformedNodes = barResult.DeformedNodes,
                AxialForces = barResult.AxialForces,
                BarStresses = barResult.Stresses,
                BarStrains = barResult.Strains,
                Utilization = utilization,
                Iterations = barResult.Iterations,
                Converged = barResult.Converged,
            };
        }

        private static FeaResult SolveSolid(
            FEA_Model model,
            List<Vector3d> loads,
            List<int>? loadNodeIdx,
            bool selfWeight)
        {
            var network = model.ToSolidNetwork();
            var solidResult = SolidSolverService.SolveForward(network, loads, loadNodeIdx, selfWeight);

            var deformedModel = new FEA_Model(model);
            for (int i = 0; i < solidResult.DeformedNodes.Length; i++)
                deformedModel.SolidMesh!.TetNodes[i] = solidResult.DeformedNodes[i];

            int ne = model.SolidMesh!.Ne;
            var utilization = solidResult.VonMises.Select((vm, e) =>
            {
                int matIdx = model.MaterialAssignment[e];
                double yieldStress = model.Materials[matIdx].YieldStress;
                return yieldStress > 0 ? vm / yieldStress : 0.0;
            }).ToArray();

            return new FeaResult
            {
                Model = model,
                DeformedModel = deformedModel,
                Displacements = solidResult.Displacements,
                Reactions = solidResult.Reactions,
                ReactionNodeIndices = solidResult.ReactionNodeIndices,
                DeformedNodes = solidResult.DeformedNodes,
                SolidStresses = solidResult.Stresses,
                VonMises = solidResult.VonMises,
                Utilization = utilization,
                Iterations = 1,
                Converged = true,
            };
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-200000000003");
    }
}
