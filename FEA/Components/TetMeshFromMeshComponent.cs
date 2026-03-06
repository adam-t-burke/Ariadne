using System;
using System.Drawing;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    /// <summary>
    /// Generates a tet4 volume mesh from a closed triangulated surface mesh
    /// using TetGen (constrained Delaunay tetrahedralization).
    /// </summary>
    public class TetMeshFromMeshComponent : GH_Component
    {
        public TetMeshFromMeshComponent()
            : base("Tet Mesh From Mesh", "TetMesh",
                "Generate a tet4 volume mesh from a closed surface mesh via TetGen",
                "Theseus-FEA", "Meshing")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "M", "Closed triangulated surface mesh", GH_ParamAccess.item);
            pManager.AddNumberParameter("Max Volume", "V", "Maximum tetrahedron volume (0 = no constraint)", GH_ParamAccess.item, 0.0);
            pManager.AddNumberParameter("Quality", "Q", "Minimum radius-edge ratio (0 = no quality, typical 1.2-2.0)", GH_ParamAccess.item, 2.0);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("SolidMesh", "SM", "Tetrahedral volume mesh", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Mesh? mesh = null;
            double maxVolume = 0.0;
            double quality = 2.0;

            if (!DA.GetData(0, ref mesh)) return;
            DA.GetData(1, ref maxVolume);
            DA.GetData(2, ref quality);

            if (mesh == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Null mesh input.");
                return;
            }

            if (!mesh.IsClosed)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Input mesh must be closed.");
                return;
            }

            maxVolume = Math.Max(0.0, maxVolume);
            quality = Math.Max(0.0, quality);

            try
            {
                string? error = Ariadne.TetGen.TetGenMesher.Tetrahedralize(
                    mesh, maxVolume, quality,
                    out var nodes, out var elements, out var boundaryFaces);

                if (error != null)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"TetGen failed: {error}");
                    return;
                }

                if (elements.Count == 0)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                        "TetGen produced no elements. Check that the mesh is valid and closed.");
                    return;
                }

                var solidMesh = new SolidMesh(nodes, elements, boundaryFaces);

                AddRuntimeMessage(GH_RuntimeMessageLevel.Remark,
                    $"Generated {solidMesh.Ne} tets with {solidMesh.Nn} nodes.");

                DA.SetData(0, solidMesh);
            }
            catch (Exception ex)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, $"Tet meshing failed: {ex.Message}");
            }
        }

        /// <summary>
        /// Tet quality: ratio of inscribed to circumscribed sphere radii,
        /// normalized so a regular tet = 1. Uses the correct formula:
        /// q = 3 * r_in / R_circum where r_in = 3V/A_surface.
        /// </summary>
        internal static double TetQuality(Point3d p0, Point3d p1, Point3d p2, Point3d p3)
        {
            var e01 = p1 - p0; var e02 = p2 - p0; var e03 = p3 - p0;
            var e12 = p2 - p1; var e13 = p3 - p1; var e23 = p3 - p2;

            double vol6 = Math.Abs(
                e01.X * (e02.Y * e03.Z - e02.Z * e03.Y) -
                e01.Y * (e02.X * e03.Z - e02.Z * e03.X) +
                e01.Z * (e02.X * e03.Y - e02.Y * e03.X));

            if (vol6 < 1e-30) return 0.0;

            double volume = vol6 / 6.0;

            double a0 = 0.5 * Vector3d.CrossProduct(e12, e13).Length;
            double a1 = 0.5 * Vector3d.CrossProduct(e02, e03).Length;
            double a2 = 0.5 * Vector3d.CrossProduct(e01, e03).Length;
            double a3 = 0.5 * Vector3d.CrossProduct(e01, e02).Length;

            double surfaceArea = a0 + a1 + a2 + a3;
            if (surfaceArea < 1e-30) return 0.0;

            double rIn = 3.0 * volume / surfaceArea;

            var n1 = Vector3d.CrossProduct(e02, e03);
            var n2 = Vector3d.CrossProduct(e03, e01);
            var n3 = Vector3d.CrossProduct(e01, e02);

            double denom = 2.0 * (e01.X * n1.X + e01.Y * n1.Y + e01.Z * n1.Z);
            if (Math.Abs(denom) < 1e-30) return 0.0;

            double l01s = e01.SquareLength, l02s = e02.SquareLength, l03s = e03.SquareLength;
            var circumOffset = (l01s * n1 + l02s * n2 + l03s * n3) / denom;
            double rCircum = circumOffset.Length;

            if (rCircum < 1e-30) return 0.0;

            double q = 3.0 * rIn / rCircum;
            return Math.Max(0.0, Math.Min(1.0, q));
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-100000000020");
    }
}
