using System;
using System.Drawing;
using Grasshopper.Kernel;

namespace Ariadne.FEA.Components
{
    public class MaterialComponent : GH_Component
    {
        public MaterialComponent()
            : base("FEA Material", "Material",
                "Define material properties (defaults to steel). Use Preset to select a built-in material.",
                "Theseus-FEA", "Setup")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddIntegerParameter("Preset", "P", "Material preset: 0=Custom, 1=Steel, 2=Concrete, 3=Aluminum, 4=Timber", GH_ParamAccess.item, 0);
            pManager.AddNumberParameter("E", "E", "Young's modulus (Pa)", GH_ParamAccess.item, 210e9);
            pManager.AddNumberParameter("Nu", "ν", "Poisson's ratio", GH_ParamAccess.item, 0.3);
            pManager.AddNumberParameter("Density", "ρ", "Density (kg/m³)", GH_ParamAccess.item, 7850.0);
            pManager.AddNumberParameter("Yield Stress", "σy", "Yield stress (Pa)", GH_ParamAccess.item, 250e6);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Material", "M", "FEA material", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            int preset = 0;
            double e = 210e9, nu = 0.3, density = 7850.0, yieldStress = 250e6;
            DA.GetData(0, ref preset);
            DA.GetData(1, ref e);
            DA.GetData(2, ref nu);
            DA.GetData(3, ref density);
            DA.GetData(4, ref yieldStress);

            if (preset > 0)
            {
                var mat = preset switch
                {
                    1 => FeaMaterial.Steel(),
                    2 => FeaMaterial.Concrete(),
                    3 => FeaMaterial.Aluminum(),
                    4 => FeaMaterial.Timber(),
                    _ => null as FeaMaterial,
                };

                if (mat == null)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Unknown preset index {preset}. Use 0-4.");
                    return;
                }

                DA.SetData(0, mat);
                return;
            }

            if (e <= 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "E must be positive."); return; }
            if (nu <= -1.0 || nu >= 0.5) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Poisson's ratio must be in (-1, 0.5)."); return; }
            if (density <= 0) { AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Density must be positive."); return; }

            DA.SetData(0, new FeaMaterial { E = e, Nu = nu, Density = density, YieldStress = yieldStress });
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-000000000002");
    }
}
