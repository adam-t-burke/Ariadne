using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Grasshopper.Kernel;
using Ariadne.Graphs;
using Rhino.Geometry;
using Theseus.Interop;

namespace Ariadne.FEA.Objectives
{
    public class FeaTargetDisplacementObj : FeaNodeObjective
    {
        public List<Vector3d> Targets { get; init; } = [];

        public override void ApplyTo(FeaSolver solver, FeaSolverContext context)
        {
            var indices = context.ResolveNodeIndices(TargetNodes);
            var flat = new double[indices.Length * 3];
            for (int i = 0; i < indices.Length; i++)
            {
                var t = i < Targets.Count ? Targets[i] : Targets[^1];
                flat[i * 3] = t.X;
                flat[i * 3 + 1] = t.Y;
                flat[i * 3 + 2] = t.Z;
            }
            solver.AddTargetDisplacement(Weight, indices, flat);
        }
    }

    public class TargetDisplacementComponent : GH_Component
    {
        public TargetDisplacementComponent()
            : base("FEA Target Displacement", "TgtDisp",
                "Target specific displacements at nodes",
                "Theseus-FEA", "Objectives")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Nodes", "N", "Target nodes (optional)", GH_ParamAccess.list);
            pManager.AddVectorParameter("Targets", "T", "Target displacement vectors", GH_ParamAccess.list);
            pManager.AddNumberParameter("Weight", "W", "Objective weight", GH_ParamAccess.item, 1.0);
            pManager[0].Optional = true;
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Objective", "Obj", "FEA target displacement objective", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            var nodes = new List<Node>();
            var targets = new List<Vector3d>();
            double weight = 1.0;
            DA.GetDataList(0, nodes);
            if (!DA.GetDataList(1, targets)) return;
            DA.GetData(2, ref weight);

            DA.SetData(0, new FeaTargetDisplacementObj
            {
                Weight = weight,
                TargetNodes = nodes.Count > 0 ? nodes : null,
                Targets = targets,
            });
        }

        protected override Bitmap Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-000000000012");
    }
}
