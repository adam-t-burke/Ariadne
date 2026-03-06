using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    /// <summary>
    /// Computes per-element quality metrics for a SolidMesh.
    /// Quality = 3 * r_in / r_circum, normalized so a regular tet = 1.
    /// </summary>
    public class MeshQualityComponent : GH_Component
    {
        public MeshQualityComponent()
            : base("Mesh Quality", "MeshQ",
                "Compute element quality metrics for a tetrahedral SolidMesh",
                "Theseus-FEA", "Utilities")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("SolidMesh", "SM", "Tetrahedral solid mesh", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Aspect Ratios", "Q", "Quality metric per element (0=degenerate, 1=regular tet)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Min Quality", "Min", "Minimum quality across all elements", GH_ParamAccess.item);
            pManager.AddNumberParameter("Max Quality", "Max", "Maximum quality across all elements", GH_ParamAccess.item);
            pManager.AddNumberParameter("Avg Quality", "Avg", "Average quality across all elements", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            SolidMesh? mesh = null;
            if (!DA.GetData(0, ref mesh)) return;
            if (mesh == null || mesh.Ne == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid or empty SolidMesh.");
                return;
            }

            var qualities = new List<double>(mesh.Ne);
            for (int i = 0; i < mesh.Ne; i++)
            {
                var e = mesh.Elements[i];
                double q = TetQuality(mesh.TetNodes[e[0]], mesh.TetNodes[e[1]], mesh.TetNodes[e[2]], mesh.TetNodes[e[3]]);
                qualities.Add(q);
            }

            DA.SetDataList(0, qualities);
            DA.SetData(1, qualities.Min());
            DA.SetData(2, qualities.Max());
            DA.SetData(3, qualities.Average());
        }

        private static double TetQuality(Point3d p0, Point3d p1, Point3d p2, Point3d p3)
            => TetMeshFromMeshComponent.TetQuality(p0, p1, p2, p3);

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-100000000022");
    }
}
