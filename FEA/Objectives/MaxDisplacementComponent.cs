using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Grasshopper.Kernel;
using Ariadne.Graphs;
using Theseus.Interop;

namespace Ariadne.FEA.Objectives
{
    public class FeaMaxDisplacementObj : FeaNodeObjective
    {
        public override void ApplyTo(FeaSolver solver, FeaSolverContext context)
        {
            var indices = context.ResolveNodeIndices(TargetNodes);
            solver.AddMaxDisplacement(Weight, indices);
        }
    }

    public class MaxDisplacementComponent : GH_Component
    {
        public MaxDisplacementComponent()
            : base("FEA Max Displacement", "MaxDisp",
                "Minimize maximum nodal displacement",
                "Theseus-FEA", "Objectives")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Nodes", "N", "Target nodes (optional, default all)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Weight", "W", "Objective weight", GH_ParamAccess.item, 1.0);
            pManager[0].Optional = true;
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Objective", "Obj", "FEA max displacement objective", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            var nodes = new List<Node>();
            double weight = 1.0;
            DA.GetDataList(0, nodes);
            DA.GetData(1, ref weight);

            DA.SetData(0, new FeaMaxDisplacementObj
            {
                Weight = weight,
                TargetNodes = nodes.Count > 0 ? nodes : null,
            });
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-000000000011");
    }
}
