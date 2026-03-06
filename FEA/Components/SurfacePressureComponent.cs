using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Grasshopper.Kernel;
using Rhino;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    /// <summary>
    /// Computes consistent nodal forces from a uniform surface pressure applied to
    /// boundary faces of a solid mesh. For each loaded triangular face with area A
    /// and outward normal n, the consistent nodal force is f_i = (p * A / 3) * n
    /// for each of the 3 face nodes. This is exact for linear (tet4) elements.
    /// (Zienkiewicz and Taylor Vol. 1, Ch. 3.7; Bathe Ch. 4.3.2)
    /// </summary>
    public class SurfacePressureComponent : GH_Component
    {
        public SurfacePressureComponent()
            : base("Surface Pressure", "Pressure",
                "Compute consistent nodal forces from uniform surface pressure on boundary faces",
                "Theseus-FEA", "Setup")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Solid Mesh", "SM", "Solid mesh with boundary faces", GH_ParamAccess.item);
            pManager.AddNumberParameter("Pressure", "P", "Pressure magnitude (Pa, positive = compression into surface)", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Face Indices", "FI", "Boundary face indices to load (optional, default = all)", GH_ParamAccess.list);
            pManager[2].Optional = true;
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Node Loads", "NL", "Canonical nodal loads for the selected boundary faces", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            SolidMesh? mesh = null;
            double pressure = 0;
            var faceIndices = new List<int>();

            if (!DA.GetData(0, ref mesh)) return;
            if (!DA.GetData(1, ref pressure)) return;
            DA.GetDataList(2, faceIndices);

            if (mesh == null || mesh.Ne == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid solid mesh.");
                return;
            }

            if (mesh.BoundaryFaces.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Mesh has no boundary faces.");
                return;
            }

            // Determine which faces to load
            List<int[]> facesToLoad;
            if (faceIndices.Count > 0)
            {
                facesToLoad = new List<int[]>(faceIndices.Count);
                foreach (int fi in faceIndices)
                {
                    if (fi >= 0 && fi < mesh.BoundaryFaces.Count)
                        facesToLoad.Add(mesh.BoundaryFaces[fi]);
                    else
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, $"Face index {fi} out of range.");
                }
            }
            else
            {
                facesToLoad = mesh.BoundaryFaces;
            }

            // Accumulate consistent nodal forces: f_i = (p * A / 3) * n per face node
            // Pressure convention: positive pressure acts inward (compression),
            // so force direction is -n (inward normal)
            var nodeForces = new Dictionary<int, Vector3d>();
            double toMeters = RhinoDoc.ActiveDoc == null
                ? 1.0
                : RhinoMath.UnitScale(RhinoDoc.ActiveDoc.ModelUnitSystem, UnitSystem.Meters);

            foreach (int[] face in facesToLoad)
            {
                int faceIndex = mesh.BoundaryFaces.IndexOf(face);
                if (!mesh.TryGetBoundaryFaceLoadData(faceIndex, out var outwardNormal, out double areaModelUnits))
                    continue;

                // Consistent nodal force: (p * A / 3) * (-n) for each of the 3 nodes
                // Negative sign: positive pressure = compression = force inward = -n
                double areaSquareMeters = areaModelUnits * toMeters * toMeters;
                var force = -pressure * areaSquareMeters / 3.0 * outwardNormal;

                for (int i = 0; i < 3; i++)
                {
                    int ni = face[i];
                    if (nodeForces.ContainsKey(ni))
                        nodeForces[ni] += force;
                    else
                        nodeForces[ni] = force;
                }
            }

            var loads = new List<FeaLoad>(nodeForces.Count);

            foreach (var kvp in nodeForces.OrderBy(k => k.Key))
            {
                loads.Add(new FeaLoad
                {
                    NodeIndex = kvp.Key,
                    Force = kvp.Value,
                });
            }

            DA.SetDataList(0, loads);
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-200000000010");
    }
}
