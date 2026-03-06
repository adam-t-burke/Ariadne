using System;
using System.Drawing;
using Grasshopper.Kernel;

namespace Ariadne.FEA.Components
{
    public class SectionComponent : GH_Component
    {
        public SectionComponent()
            : base("FEA Section", "Section",
                "Define cross-section properties",
                "Theseus-FEA", "Setup")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Area", "A", "Cross-sectional area (m²)", GH_ParamAccess.item, 0.01);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Section", "S", "FEA section", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double area = 0.01;
            DA.GetData(0, ref area);

            if (area <= 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Area must be positive."); return; }

            DA.SetData(0, new FeaSection { Area = area });
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-000000000003");
    }
}
