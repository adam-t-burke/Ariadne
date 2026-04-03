using System;
using System.Drawing;
using Grasshopper.Kernel;

namespace Ariadne.FEA.Components
{
    public class SectionComponent : GH_Component
    {
        public SectionComponent()
            : base("FEA Section", "Section",
                "Define cross-section properties for bar/beam/shell elements",
                "Theseus-FEA", "Setup")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Area", "A", "Cross-sectional area (m²)", GH_ParamAccess.item, 0.01);
            pManager.AddNumberParameter("Iy", "Iy", "Second moment about local y for bending in x-z plane (m⁴)", GH_ParamAccess.item, 0.0);
            pManager.AddNumberParameter("Iz", "Iz", "Second moment about local z for bending in x-y plane (m⁴)", GH_ParamAccess.item, 0.0);
            pManager.AddNumberParameter("J", "J", "Torsion constant (m⁴)", GH_ParamAccess.item, 0.0);
            pManager.AddNumberParameter("Asy", "Asy", "Shear area in y (Timoshenko only; 0 = auto)", GH_ParamAccess.item, 0.0);
            pManager.AddNumberParameter("Asz", "Asz", "Shear area in z (Timoshenko only; 0 = auto)", GH_ParamAccess.item, 0.0);
            pManager.AddNumberParameter("Offset", "Off", "Shell mid-surface offset (m)", GH_ParamAccess.item, 0.0);

            for (int i = 1; i <= 6; i++)
                pManager[i].Optional = true;
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Section", "S", "FEA section", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double area = 0.01, iy = 0, iz = 0, j = 0, asy = 0, asz = 0, offset = 0;

            DA.GetData(0, ref area);
            DA.GetData(1, ref iy);
            DA.GetData(2, ref iz);
            DA.GetData(3, ref j);
            DA.GetData(4, ref asy);
            DA.GetData(5, ref asz);
            DA.GetData(6, ref offset);

            if (area <= 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Area must be positive.");
                return;
            }

            DA.SetData(0, new FeaSection
            {
                Area = area,
                Iy = iy,
                Iz = iz,
                J = j,
                Asy = asy,
                Asz = asz,
                Offset = offset,
            });
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-000000000003");
    }
}
