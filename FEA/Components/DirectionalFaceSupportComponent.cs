using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;

namespace Ariadne.FEA.Components
{
    public class DirectionalFaceSupportComponent : GH_Component
    {
        public DirectionalFaceSupportComponent()
            : base("Directional Face Support", "FaceSup",
                "Create a directional support region on selected solid boundary faces",
                "Theseus-FEA", "Setup")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Solid Mesh", "SM", "Solid mesh with boundary faces", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Face Indices", "FI", "Boundary face indices to constrain (optional, default = all)", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Fix X", "X", "Constrain X translation", GH_ParamAccess.item, false);
            pManager.AddBooleanParameter("Fix Y", "Y", "Constrain Y translation", GH_ParamAccess.item, false);
            pManager.AddBooleanParameter("Fix Z", "Z", "Constrain Z translation", GH_ParamAccess.item, true);
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
            bool fixX = false, fixY = false, fixZ = true;

            if (!DA.GetData(0, ref mesh)) return;
            DA.GetDataList(1, faceIndices);
            DA.GetData(2, ref fixX);
            DA.GetData(3, ref fixY);
            DA.GetData(4, ref fixZ);

            if (mesh == null || mesh.BoundaryFaces.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Solid mesh must contain boundary faces.");
                return;
            }

            if (!fixX && !fixY && !fixZ)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "At least one translational DOF must be constrained.");
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
                FixX = fixX,
                FixY = fixY,
                FixZ = fixZ,
            });
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-200000000013");
    }
}
