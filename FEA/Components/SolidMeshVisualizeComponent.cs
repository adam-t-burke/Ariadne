using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Grasshopper.Kernel.Types;
using Grasshopper.Kernel;
using Rhino.Display;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    public class SolidMeshVisualizeComponent : GH_Component
    {
        private SolidMesh? _mesh;
        private FeaResult? _result;
        private int _colorProp;
        private bool _showEdges = true;
        private bool _showFaces = true;
        private double _deformScale = 1.0;
        private Color _colorMin = Color.Blue;
        private Color _colorMed = Color.White;
        private Color _colorMax = Color.Red;
        private Plane _clipPlane = Plane.Unset;
        private bool _useClip;

        public SolidMeshVisualizeComponent()
            : base("Solid Mesh Visualize", "SolidViz",
                "Visualize solid mesh with stress coloring and deformation",
                "Theseus-FEA", "Utilities")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddBooleanParameter("Show", "Show", "Enable visualization", GH_ParamAccess.item, true);
            pManager.AddGenericParameter("Mesh", "M", "Solid mesh or FEA result", GH_ParamAccess.item);
            pManager.AddGenericParameter("Results", "Res", "Solid FEA result (optional)", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Color Property", "CP", "0=None, 1=VonMises, 2=Displacement", GH_ParamAccess.item, 0);
            pManager.AddBooleanParameter("Show Edges", "SE", "Show element edges", GH_ParamAccess.item, true);
            pManager.AddBooleanParameter("Show Faces", "SF", "Show boundary faces", GH_ParamAccess.item, true);
            pManager.AddNumberParameter("Deform Scale", "DS", "Deformation scale factor", GH_ParamAccess.item, 1.0);
            pManager.AddColourParameter("Color Min", "CMin", "Color for minimum values", GH_ParamAccess.item, Color.Blue);
            pManager.AddColourParameter("Color Med", "CMed", "Color for median values", GH_ParamAccess.item, Color.White);
            pManager.AddColourParameter("Color Max", "CMax", "Color for maximum values", GH_ParamAccess.item, Color.Red);
            pManager.AddPlaneParameter("Clip Plane", "Clip", "Clipping plane (optional)", GH_ParamAccess.item);

            pManager[2].Optional = true;
            pManager[10].Optional = true;
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager) { }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            bool show = true;
            DA.GetData(0, ref show);
            if (!show) { _mesh = null; _result = null; return; }

            object? meshObj = null;
            if (!DA.GetData(1, ref meshObj)) return;
            if (meshObj is GH_ObjectWrapper meshWrapper)
                meshObj = meshWrapper.Value;

            _mesh = meshObj as SolidMesh;
            _result = null;

            if (_mesh == null && meshObj is FeaResult meshResult)
            {
                if (!meshResult.Model.HasSolidElements)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "FeaResult has no solid mesh.");
                    return;
                }

                _mesh = meshResult.Model.SolidMesh;
                _result = meshResult;
            }

            object? resultObj = null;
            if (DA.GetData(2, ref resultObj))
            {
                if (resultObj is GH_ObjectWrapper resultWrapper)
                    resultObj = resultWrapper.Value;

                if (resultObj is FeaResult result)
                {
                    if (!result.Model.HasSolidElements)
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Results input must contain a solid FeaResult.");
                        return;
                    }

                    _result = result;
                    _mesh = result.Model.SolidMesh;
                }
            }

            DA.GetData(3, ref _colorProp);
            DA.GetData(4, ref _showEdges);
            DA.GetData(5, ref _showFaces);
            DA.GetData(6, ref _deformScale);
            DA.GetData(7, ref _colorMin);
            DA.GetData(8, ref _colorMed);
            DA.GetData(9, ref _colorMax);

            _clipPlane = Plane.Unset;
            _useClip = DA.GetData(10, ref _clipPlane);
        }

        public override void DrawViewportMeshes(IGH_PreviewArgs args)
        {
            if (_mesh == null || !_showFaces) return;

            var displayMesh = BuildDisplayMesh();
            if (displayMesh == null) return;

            var mat = new DisplayMaterial(_colorMed, 0.5);
            args.Display.DrawMeshShaded(displayMesh, mat);
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            if (_mesh == null) return;

            if (_showEdges)
            {
                var positions = GetDisplayPositions();
                var edgeColor = Color.FromArgb(100, 100, 100);

                foreach (var elem in _mesh.Elements)
                {
                    var p = new Point3d[4];
                    for (int i = 0; i < 4; i++) p[i] = positions[elem[i]];

                    if (_useClip && !PassesClip(p)) continue;

                    for (int i = 0; i < 4; i++)
                        for (int j = i + 1; j < 4; j++)
                            args.Display.DrawLine(p[i], p[j], edgeColor, 1);
                }
            }
        }

        private Mesh BuildDisplayMesh()
        {
            var positions = GetDisplayPositions();
            var rhinoMesh = new Mesh();
            rhinoMesh.Vertices.AddVertices(positions);

            double[]? values = null;
            if (_result != null && _colorProp > 0)
            {
                values = _colorProp switch
                {
                    1 => _result.VonMises,
                    2 => positions.Select((_, i) => _result.Displacements[i].Length).ToArray(),
                    _ => null,
                };
            }

            if (_mesh == null) return rhinoMesh;

            foreach (var face in _mesh.BoundaryFaces)
            {
                if (_useClip)
                {
                    var centroid = (positions[face[0]] + positions[face[1]] + positions[face[2]]) / 3.0;
                    if (_clipPlane.IsValid && _clipPlane.DistanceTo((Point3d)centroid) < 0) continue;
                }
                rhinoMesh.Faces.AddFace(face[0], face[1], face[2]);
            }

            if (values != null && values.Length > 0)
            {
                double[] nodeValues;
                if (values.Length == _mesh.Ne)
                {
                    nodeValues = new double[_mesh.Nn];
                    var counts = new int[_mesh.Nn];
                    for (int e = 0; e < _mesh.Ne; e++)
                    {
                        var elem = _mesh.Elements[e];
                        for (int j = 0; j < 4; j++)
                        {
                            nodeValues[elem[j]] += values[e];
                            counts[elem[j]]++;
                        }
                    }
                    for (int n = 0; n < _mesh.Nn; n++)
                        if (counts[n] > 0) nodeValues[n] /= counts[n];
                }
                else
                {
                    nodeValues = values;
                }

                double vMin = nodeValues.Min();
                double vMax = nodeValues.Max();
                double range = Math.Max(Math.Abs(vMax - vMin), 1e-12);

                for (int n = 0; n < positions.Length; n++)
                {
                    double t = (nodeValues[n] - vMin) / range;
                    Color c = t < 0.5
                        ? Lerp(_colorMin, _colorMed, t * 2.0)
                        : Lerp(_colorMed, _colorMax, (t - 0.5) * 2.0);
                    rhinoMesh.VertexColors.Add(c);
                }
            }

            rhinoMesh.Normals.ComputeNormals();
            return rhinoMesh;
        }

        private Point3d[] GetDisplayPositions()
        {
            if (_mesh == null) return [];
            if (_result?.DeformedNodes != null && _deformScale == 1.0)
                return _result.DeformedNodes.ToArray();

            var positions = _mesh.TetNodes.ToArray();
            if (_result?.Displacements != null && _deformScale != 0)
            {
                for (int i = 0; i < positions.Length && i < _result.Displacements.Length; i++)
                    positions[i] += _result.Displacements[i] * _deformScale;
            }
            return positions;
        }

        private bool PassesClip(Point3d[] pts)
        {
            if (!_clipPlane.IsValid) return true;
            return pts.Any(p => _clipPlane.DistanceTo(p) >= 0);
        }

        private static Color Lerp(Color a, Color b, double t)
        {
            t = Math.Max(0, Math.Min(1, t));
            return Color.FromArgb(
                (int)(a.R + (b.R - a.R) * t),
                (int)(a.G + (b.G - a.G) * t),
                (int)(a.B + (b.B - a.B) * t));
        }

        public override BoundingBox ClippingBox
        {
            get
            {
                if (_mesh == null) return BoundingBox.Empty;
                var bb = new BoundingBox(GetDisplayPositions());
                bb.Inflate(1.0);
                return bb;
            }
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-100000000023");
    }
}
