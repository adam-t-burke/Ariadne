using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Display;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    public class DeconstructShellMeshComponent : GH_Component
    {
        private ShellMesh? _mesh;
        private bool _showEdges = true;

        public DeconstructShellMeshComponent()
            : base("Deconstruct Shell Mesh", "Decon Shell",
                "Decompose a shell mesh into nodes, elements, boundary edges, and thicknesses. Displays mesh with optional edges.",
                "Theseus-FEA", "Meshing")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Shell Mesh", "SM", "Prepared shell mesh", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Preview Edges", "PE", "Show element edges in viewport", GH_ParamAccess.item, true);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Nodes", "N", "Node objects with Index for FEA setup", GH_ParamAccess.list);
            pManager.AddMeshParameter("Mesh", "M", "Triangulated mesh for visualization", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Elements", "E", "Triangle connectivity (3 node indices per element, tree)", GH_ParamAccess.tree);
            pManager.AddTextParameter("Boundary Edges", "BE", "Naked edge node pairs (i-j)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Thicknesses", "T", "Per-node thickness values", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            ShellMesh? shellMesh = null;
            bool previewEdges = true;

            if (!DA.GetData(0, ref shellMesh)) return;
            DA.GetData(1, ref previewEdges);

            if (shellMesh == null || shellMesh.Ne == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid or empty shell mesh.");
                return;
            }

            _mesh = shellMesh;
            _showEdges = previewEdges;

            var nodes = new List<Graphs.Node>(shellMesh.Nn);
            for (int i = 0; i < shellMesh.Nn; i++)
            {
                nodes.Add(new Graphs.Node
                {
                    Value = shellMesh.Nodes[i],
                    Index = i
                });
            }

            var triMesh = shellMesh.BuildTriangulatedMesh();

            var elementTree = new GH_Structure<GH_Integer>();
            for (int e = 0; e < shellMesh.Ne; e++)
            {
                var path = new GH_Path(e);
                foreach (int idx in shellMesh.Elements[e])
                    elementTree.Append(new GH_Integer(idx), path);
            }

            var boundaryEdgeStrings = shellMesh.BoundaryEdges
                .Select(e => $"{e.Item1}-{e.Item2}")
                .ToList();

            DA.SetDataList(0, nodes);
            DA.SetData(1, triMesh);
            DA.SetDataTree(2, elementTree);
            DA.SetDataList(3, boundaryEdgeStrings);
            DA.SetDataList(4, shellMesh.NodeThicknesses);
        }

        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            if (_mesh == null) return;
            var rhinoMesh = _mesh.BuildTriangulatedMesh();
            var mat = new DisplayMaterial(Color.White, 0.5);
            args.Display.DrawMeshShaded(rhinoMesh, mat);
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            if (_mesh == null || !_showEdges) return;

            var edgeColor = Color.FromArgb(100, 100, 100);
            foreach (var tri in _mesh.Elements)
            {
                var p0 = _mesh.Nodes[tri[0]];
                var p1 = _mesh.Nodes[tri[1]];
                var p2 = _mesh.Nodes[tri[2]];
                args.Display.DrawLine(p0, p1, edgeColor, 1);
                args.Display.DrawLine(p1, p2, edgeColor, 1);
                args.Display.DrawLine(p2, p0, edgeColor, 1);
            }
        }

        public override BoundingBox ClippingBox
        {
            get
            {
                if (_mesh == null || _mesh.Nn == 0) return BoundingBox.Empty;
                var bb = new BoundingBox(_mesh.Nodes);
                bb.Inflate(1.0);
                return bb;
            }
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-100000000051");
    }
}
