using System;
using System.Drawing;
using Grasshopper.Kernel;
using Theseus.Interop;

namespace Ariadne.FEA.Objectives
{
    public class FeaMinWeightObj : FeaObjective
    {
        public override void ApplyTo(FeaSolver solver, FeaSolverContext context)
        {
            solver.AddMinWeight(Weight);
        }
    }

    public class MinWeightComponent : GH_Component
    {
        public MinWeightComponent()
            : base("FEA Min Weight", "MinWeight",
                "Minimize total structural weight",
                "Theseus-FEA", "Objectives")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Weight", "W", "Objective weight", GH_ParamAccess.item, 1.0);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Objective", "Obj", "FEA min weight objective", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double weight = 1.0;
            DA.GetData(0, ref weight);
            DA.SetData(0, new FeaMinWeightObj { Weight = weight });
        }

        protected override Bitmap Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-000000000013");
    }
}
