using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Display;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    public class DeconstructSolidMeshComponent : GH_Component
    {
        private SolidMesh? _mesh;
        private bool _showEdges = true;

        public DeconstructSolidMeshComponent()
            : base("Deconstruct Solid Mesh", "Decon Solid",
                "Decompose solid mesh into nodes, boundary mesh, and elements. Displays mesh with optional edges.",
                "Theseus-FEA", "Meshing")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Solid Mesh", "SM", "Tetrahedral solid mesh", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Preview Edges", "PE", "Show element edges in viewport", GH_ParamAccess.item, true);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Nodes", "N", "Node objects (GH_Point) with Index for FEA setup", GH_ParamAccess.list);
            pManager.AddMeshParameter("Boundary Mesh", "BM", "Surface mesh for native Deconstruct Mesh, face selection", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Elements", "E", "Tet connectivity (4 node indices per element, tree)", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            SolidMesh? solidMesh = null;
            bool previewEdges = true;

            if (!DA.GetData(0, ref solidMesh)) return;
            DA.GetData(1, ref previewEdges);

            if (solidMesh == null || solidMesh.Ne == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid or empty solid mesh.");
                return;
            }

            _mesh = solidMesh;
            _showEdges = previewEdges;

            var nodes = new List<Graphs.Node>();
            for (int i = 0; i < solidMesh.TetNodes.Count; i++)
            {
                nodes.Add(new Graphs.Node
                {
                    Value = solidMesh.TetNodes[i],
                    Index = i
                });
            }

            var boundaryMesh = solidMesh.Value;
            var elementTree = new GH_Structure<GH_Integer>();
            for (int e = 0; e < solidMesh.Elements.Count; e++)
            {
                var path = new GH_Path(e);
                foreach (int idx in solidMesh.Elements[e])
                    elementTree.Append(new GH_Integer(idx), path);
            }

            DA.SetDataList(0, nodes);
            DA.SetData(1, boundaryMesh);
            DA.SetDataTree(2, elementTree);
        }

        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            if (_mesh == null) return;

            var rhinoMesh = new Mesh();
            rhinoMesh.Vertices.AddVertices(_mesh.TetNodes);
            foreach (var face in _mesh.BoundaryFaces)
                rhinoMesh.Faces.AddFace(face[0], face[1], face[2]);
            rhinoMesh.Normals.ComputeNormals();

            var mat = new DisplayMaterial(System.Drawing.Color.White, 0.5);
            args.Display.DrawMeshShaded(rhinoMesh, mat);
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            if (_mesh == null || !_showEdges) return;

            var positions = _mesh.TetNodes;
            var edgeColor = System.Drawing.Color.FromArgb(100, 100, 100);

            foreach (var elem in _mesh.Elements)
            {
                var p0 = positions[elem[0]];
                var p1 = positions[elem[1]];
                var p2 = positions[elem[2]];
                var p3 = positions[elem[3]];

                args.Display.DrawLine(p0, p1, edgeColor, 1);
                args.Display.DrawLine(p0, p2, edgeColor, 1);
                args.Display.DrawLine(p0, p3, edgeColor, 1);
                args.Display.DrawLine(p1, p2, edgeColor, 1);
                args.Display.DrawLine(p1, p3, edgeColor, 1);
                args.Display.DrawLine(p2, p3, edgeColor, 1);
            }
        }

        public override BoundingBox ClippingBox
        {
            get
            {
                if (_mesh == null) return BoundingBox.Empty;
                var bb = new BoundingBox(_mesh.TetNodes);
                bb.Inflate(1.0);
                return bb;
            }
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-100000000030");
    }
}
