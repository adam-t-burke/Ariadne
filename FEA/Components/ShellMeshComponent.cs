using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    public class ShellMeshComponent : GH_Component
    {
        public ShellMeshComponent()
            : base("Shell Mesh", "ShellMesh",
                "Prepare a shell mesh from a Rhino mesh: triangulate quads, assign per-node thickness, and expose nodes for load/support wiring",
                "Theseus-FEA", "Meshing")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("Mesh", "M", "Rhino mesh (triangles and/or quads)", GH_ParamAccess.item);
            pManager.AddNumberParameter("Thickness", "T", "Shell thickness: single value or per-vertex list", GH_ParamAccess.list);
            pManager.AddGenericParameter("Material", "Mat", "Material properties (default: Steel)", GH_ParamAccess.item);
            pManager[2].Optional = true;
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Shell Mesh", "SM", "Prepared shell mesh for FEA", GH_ParamAccess.item);
            pManager.AddGenericParameter("Nodes", "N", "Node objects with Index (wire into Point Load or Support)", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Mesh? mesh = null;
            var thicknesses = new List<double>();

            if (!DA.GetData(0, ref mesh)) return;
            if (!DA.GetDataList(1, thicknesses) || thicknesses.Count == 0) return;

            if (mesh == null || mesh.Faces.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Input mesh is null or has no faces.");
                return;
            }

            if (thicknesses.Any(t => t <= 0))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Thickness values must be positive.");
                return;
            }

            int nn = mesh.Vertices.Count;
            double[] nodeThicknesses;

            if (thicknesses.Count == 1)
            {
                nodeThicknesses = Enumerable.Repeat(thicknesses[0], nn).ToArray();
            }
            else if (thicknesses.Count == nn)
            {
                nodeThicknesses = thicknesses.ToArray();
            }
            else
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                    $"Thickness list has {thicknesses.Count} values but mesh has {nn} vertices. Provide 1 value or {nn} values.");
                return;
            }

            FeaMaterial material = FeaMaterial.Steel();
            {
                object? matInput = null;
                if (DA.GetData(2, ref matInput) && matInput != null)
                {
                    if (matInput is GH_ObjectWrapper w) matInput = w.Value;
                    if (matInput is FeaMaterial m) material = m;
                }
            }

            var shellMesh = new ShellMesh(mesh, nodeThicknesses, material);

            var nodes = new List<Graphs.Node>(nn);
            for (int i = 0; i < nn; i++)
            {
                nodes.Add(new Graphs.Node
                {
                    Value = shellMesh.Nodes[i],
                    Index = i
                });
            }

            AddRuntimeMessage(GH_RuntimeMessageLevel.Remark,
                $"Shell mesh: {shellMesh.Nn} nodes, {shellMesh.Ne} triangles" +
                (shellMesh.Ne != shellMesh.OriginalFaceCount
                    ? $" ({shellMesh.OriginalFaceCount} original faces, {shellMesh.Ne - shellMesh.OriginalFaceCount} quads split)"
                    : ""));

            DA.SetData(0, shellMesh);
            DA.SetDataList(1, nodes);
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-100000000050");
    }
}
