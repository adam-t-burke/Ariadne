using System;
using System.Drawing;
using Grasshopper.Kernel;

namespace Ariadne.FEA.Components
{
    public class MaterialComponent : GH_Component
    {
        public MaterialComponent()
            : base("FEA Material", "Material",
                "Define material properties (defaults to steel)",
                "Theseus-FEA", "Design")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("E", "E", "Young's modulus (Pa)", GH_ParamAccess.item, 210e9);
            pManager.AddNumberParameter("Density", "ρ", "Density (kg/m³)", GH_ParamAccess.item, 7850.0);
            pManager.AddNumberParameter("Yield Stress", "σy", "Yield stress (Pa)", GH_ParamAccess.item, 250e6);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Material", "M", "FEA material", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double e = 210e9, density = 7850.0, yieldStress = 250e6;
            DA.GetData(0, ref e);
            DA.GetData(1, ref density);
            DA.GetData(2, ref yieldStress);

            if (e <= 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "E must be positive."); return; }
            if (density <= 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Density must be positive."); return; }

            DA.SetData(0, new FeaMaterial { E = e, Density = density, YieldStress = yieldStress });
        }

        protected override Bitmap Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-000000000002");
    }
}
