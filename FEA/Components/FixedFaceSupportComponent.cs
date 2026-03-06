using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;

namespace Ariadne.FEA.Components
{
    public class FixedFaceSupportComponent : GH_Component
    {
        public FixedFaceSupportComponent()
            : base("Fixed Face Support", "FaceFix",
                "Create a fully fixed support region on selected solid boundary faces",
                "Theseus-FEA", "Setup")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Solid Mesh", "SM", "Solid mesh with boundary faces", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Face Indices", "FI", "Boundary face indices to constrain (optional, default = all)", GH_ParamAccess.list);
            pManager[1].Optional = true;
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Support", "S", "Solid face support region", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            SolidMesh? mesh = null;
            var faceIndices = new List<int>();

            if (!DA.GetData(0, ref mesh)) return;
            DA.GetDataList(1, faceIndices);

            if (mesh == null || mesh.BoundaryFaces.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Solid mesh must contain boundary faces.");
                return;
            }

            var selectedFaces = faceIndices.Count > 0
                ? faceIndices
                : new List<int>(System.Linq.Enumerable.Range(0, mesh.BoundaryFaces.Count));
            var nodeIndices = mesh.GetBoundaryFaceNodeIndices(selectedFaces);

            if (nodeIndices.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "No valid boundary faces were selected.");
                return;
            }

            DA.SetData(0, new FeaSolidSupportRegion
            {
                NodeIndices = nodeIndices,
                FaceIndices = selectedFaces,
                FixX = true,
                FixY = true,
                FixZ = true,
            });
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-200000000012");
    }
}
