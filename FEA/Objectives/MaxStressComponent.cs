using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Grasshopper.Kernel;
using Ariadne.Graphs;
using Theseus.Interop;

namespace Ariadne.FEA.Objectives
{
    public class FeaMaxStressObj : FeaEdgeObjective
    {
        public List<double> Thresholds { get; init; } = [];

        public override void ApplyTo(FeaSolver solver, FeaSolverContext context)
        {
            var indices = context.ResolveEdgeIndices(TargetEdges);
            var thresholds = new double[indices.Length];
            for (int i = 0; i < indices.Length; i++)
                thresholds[i] = i < Thresholds.Count ? Thresholds[i] : Thresholds[^1];
            solver.AddMaxStress(Weight, indices, thresholds, Sharpness);
        }
    }

    public class MaxStressComponent : GH_Component
    {
        public MaxStressComponent()
            : base("FEA Max Stress", "MaxStress",
                "Penalize elements exceeding stress threshold",
                "Theseus-FEA", "Objectives")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Edges", "E", "Target edges (optional, default all)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Threshold", "T", "Stress threshold (Pa)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Weight", "W", "Objective weight", GH_ParamAccess.item, 1.0);
            pManager.AddNumberParameter("Sharpness", "S", "Barrier sharpness", GH_ParamAccess.item, 20.0);
            pManager[0].Optional = true;
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Objective", "Obj", "FEA max stress objective", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            var edges = new List<Edge>();
            var thresholds = new List<double>();
            double weight = 1.0, sharpness = 20.0;

            DA.GetDataList(0, edges);
            if (!DA.GetDataList(1, thresholds)) return;
            DA.GetData(2, ref weight);
            DA.GetData(3, ref sharpness);

            if (thresholds.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "At least one stress threshold is required.");
                return;
            }

            DA.SetData(0, new FeaMaxStressObj
            {
                Weight = weight,
                TargetEdges = edges.Count > 0 ? edges : null,
                Thresholds = thresholds,
                Sharpness = sharpness,
            });
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-000000000014");
    }
}
