using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Ariadne.FEA;
using Ariadne.Graphs;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Display;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    public class FeaVisualizeComponent : GH_Component
    {
        private Graph? _graph;
        private Vector3d[]? _displacements;
        private double[]? _axialForces;
        private double[]? _stresses;
        private double[]? _utilization;
        private Point3d[]? _deformedNodes;
        private int _colorProperty;
        private bool _showDeformed = true;
        private double _deformScale = 1.0;
        private Color _colorMin = Color.Blue;
        private Color _colorMed = Color.White;
        private Color _colorMax = Color.Red;
        private double _lineThickness = 2.0;

        public FeaVisualizeComponent()
            : base("FEA Visualize", "FEA Viz",
                "Visualize FEA results with color-coded bar elements",
                "Theseus-FEA", "Utilities")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddBooleanParameter("Show", "Show", "Enable visualization", GH_ParamAccess.item, true);
            pManager.AddGenericParameter("Results", "R", "FEA solve result (bar only)", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Color Property", "CP", "0=AxialForce, 1=Stress, 2=Utilization, 3=Displacement", GH_ParamAccess.item, 2);
            pManager.AddBooleanParameter("Show Deformed", "Def", "Show deformed geometry", GH_ParamAccess.item, true);
            pManager.AddNumberParameter("Deform Scale", "DS", "Deformation scale factor", GH_ParamAccess.item, 1.0);
            pManager.AddColourParameter("Color Min", "CMin", "Color for minimum values", GH_ParamAccess.item, Color.Blue);
            pManager.AddColourParameter("Color Med", "CMed", "Color for zero/median values", GH_ParamAccess.item, Color.White);
            pManager.AddColourParameter("Color Max", "CMax", "Color for maximum values", GH_ParamAccess.item, Color.Red);
            pManager.AddNumberParameter("Line Thickness", "LT", "Display line thickness", GH_ParamAccess.item, 2.0);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            bool show = true;
            DA.GetData(0, ref show);
            if (!show)
            {
                _graph = null;
                _displacements = null;
                _axialForces = null;
                _stresses = null;
                _utilization = null;
                _deformedNodes = null;
                return;
            }

            object? input = null;
            if (!DA.GetData(1, ref input)) return;

            // Unwrap GH_ObjectWrapper when Grasshopper wraps custom types for generic parameters
            if (input is GH_ObjectWrapper wrapper)
                input = wrapper.Value;

            if (input is FeaResult result)
            {
                if (!result.Model.HasBarElements)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "FeaResult has no bar elements. Use Solid Mesh Visualize for solid results.");
                    _graph = null;
                    return;
                }
                _graph = _showDeformed ? result.DeformedModel.BarGraph! : result.Model.BarGraph!;
                _displacements = result.Displacements;
                _axialForces = result.AxialForces ?? [];
                _stresses = result.BarStresses ?? [];
                _utilization = result.Utilization;
                _deformedNodes = result.DeformedNodes;
            }
            else
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Input must be a bar FeaResult.");
                _graph = null;
                return;
            }

            DA.GetData(2, ref _colorProperty);
            DA.GetData(3, ref _showDeformed);
            DA.GetData(4, ref _deformScale);
            DA.GetData(5, ref _colorMin);
            DA.GetData(6, ref _colorMed);
            DA.GetData(7, ref _colorMax);
            DA.GetData(8, ref _lineThickness);

            // Re-fetch graph with correct deformed state after _showDeformed is read
            if (input is FeaResult u)
            {
                _graph = _showDeformed ? u.DeformedModel.BarGraph! : u.Model.BarGraph!;
            }
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            if (_graph == null || _utilization == null) return;

            int ne = _graph.Ne;
            if (ne == 0) return;

            double[] values = _colorProperty switch
            {
                0 => _axialForces ?? [],
                1 => _stresses ?? [],
                2 => _utilization,
                3 => GetDisplacementMagnitudes(),
                _ => _utilization,
            };

            if (values.Length == 0) return;

            double vMin = values.Min();
            double vMax = values.Max();
            double range = Math.Max(Math.Abs(vMax - vMin), 1e-12);

            for (int e = 0; e < ne; e++)
            {
                var edge = _graph.Edges[e];
                Point3d p0, p1;

                if (_showDeformed && _deformedNodes != null && _deformedNodes.Length > 0)
                {
                    int ni = edge.Start.Index, nj = edge.End.Index;
                    if (ni < _deformedNodes.Length && nj < _deformedNodes.Length)
                    {
                        p0 = _deformedNodes[ni];
                        p1 = _deformedNodes[nj];
                    }
                    else
                    {
                        p0 = _graph.Nodes[ni].Value;
                        p1 = _graph.Nodes[nj].Value;
                    }
                }
                else
                {
                    p0 = _graph.Nodes[edge.Start.Index].Value;
                    p1 = _graph.Nodes[edge.End.Index].Value;
                }

                double t = (values[e] - vMin) / range;
                Color color = t < 0.5
                    ? InterpolateColor(_colorMin, _colorMed, t * 2.0)
                    : InterpolateColor(_colorMed, _colorMax, (t - 0.5) * 2.0);

                args.Display.DrawLine(new Line(p0, p1), color, (int)_lineThickness);
            }
        }

        private double[] GetDisplacementMagnitudes()
        {
            if (_displacements == null || _graph == null) return [];
            int nn = _graph.Nn;
            var mags = new double[nn];
            for (int i = 0; i < nn && i < _displacements.Length; i++)
                mags[i] = _displacements[i].Length;

            int ne = _graph.Ne;
            var edgeMags = new double[ne];
            for (int e = 0; e < ne; e++)
            {
                var edge = _graph.Edges[e];
                int i = edge.Start.Index, j = edge.End.Index;
                edgeMags[e] = (mags[i] + mags[j]) * 0.5;
            }
            return edgeMags;
        }

        private static Color InterpolateColor(Color a, Color b, double t)
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
                if (_deformedNodes == null || _deformedNodes.Length == 0) return BoundingBox.Empty;
                var bb = new BoundingBox(_deformedNodes);
                bb.Inflate(1.0);
                return bb;
            }
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-000000000021");
    }
}
