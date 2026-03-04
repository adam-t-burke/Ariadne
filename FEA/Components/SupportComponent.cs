using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    public class SupportComponent : GH_Component
    {
        public SupportComponent()
            : base("FEA Support", "Support",
                "Define a support boundary condition with per-DOF constraints",
                "Theseus-FEA", "Design")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Point", "P", "Support location", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Fix X", "X", "Constrain X translation", GH_ParamAccess.item, true);
            pManager.AddBooleanParameter("Fix Y", "Y", "Constrain Y translation", GH_ParamAccess.item, true);
            pManager.AddBooleanParameter("Fix Z", "Z", "Constrain Z translation", GH_ParamAccess.item, true);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Support", "S", "FEA support", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Point3d pt = Point3d.Origin;
            bool fixX = true, fixY = true, fixZ = true;

            if (!DA.GetData(0, ref pt)) return;
            DA.GetData(1, ref fixX);
            DA.GetData(2, ref fixY);
            DA.GetData(3, ref fixZ);

            if (!fixX && !fixY && !fixZ)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Support has no constrained DOFs.");
                return;
            }

            DA.SetData(0, new FeaSupport
            {
                Location = pt,
                FixX = fixX,
                FixY = fixY,
                FixZ = fixZ,
            });
        }

        protected override Bitmap Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-000000000001");
    }
}
