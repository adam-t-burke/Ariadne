using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Grasshopper.Kernel;
using Rhino.Display;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    public class FeaVisualizeComponent : GH_Component
    {
        private FeaSolveResult _result;
        private int _colorProperty;
        private bool _showDeformed = true;
        private double _deformScale = 1.0;
        private Color _colorMin = Color.Blue;
        private Color _colorMed = Color.White;
        private Color _colorMax = Color.Red;
        private double _lineThickness = 2.0;

        public FeaVisualizeComponent()
            : base("FEA Visualize", "FEA Viz",
                "Visualize FEA results with color-coded elements",
                "Theseus-FEA", "Utilities")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddBooleanParameter("Show", "Show", "Enable visualization", GH_ParamAccess.item, true);
            pManager.AddGenericParameter("Results", "R", "FEA solve result", GH_ParamAccess.item);
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
            if (!show) { _result = null; return; }

            FeaSolveResult result = null;
            if (!DA.GetData(1, ref result)) return;

            _result = result;
            DA.GetData(2, ref _colorProperty);
            DA.GetData(3, ref _showDeformed);
            DA.GetData(4, ref _deformScale);
            DA.GetData(5, ref _colorMin);
            DA.GetData(6, ref _colorMed);
            DA.GetData(7, ref _colorMax);
            DA.GetData(8, ref _lineThickness);
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            if (_result == null) return;

            var graph = _result.Network.Graph;
            int ne = graph.Ne;

            double[] values = _colorProperty switch
            {
                0 => _result.AxialForces,
                1 => _result.Stresses,
                2 => _result.Utilization,
                3 => GetDisplacementMagnitudes(),
                _ => _result.Utilization,
            };

            if (values == null || values.Length == 0) return;

            double vMin = values.Min();
            double vMax = values.Max();
            double range = Math.Max(Math.Abs(vMax - vMin), 1e-12);

            for (int e = 0; e < ne; e++)
            {
                var edge = graph.Edges[e];
                Point3d p0, p1;

                if (_showDeformed && _result.DeformedNodes != null)
                {
                    int ni = edge.Start.Index, nj = edge.End.Index;
                    var orig0 = graph.Nodes[ni].Value;
                    var orig1 = graph.Nodes[nj].Value;
                    var def0 = _result.DeformedNodes[ni];
                    var def1 = _result.DeformedNodes[nj];
                    p0 = orig0 + (def0 - orig0) * _deformScale;
                    p1 = orig1 + (def1 - orig1) * _deformScale;
                }
                else
                {
                    p0 = graph.Nodes[edge.Start.Index].Value;
                    p1 = graph.Nodes[edge.End.Index].Value;
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
            if (_result?.Displacements == null) return [];
            int nn = _result.Network.Graph.Nn;
            var mags = new double[nn];
            for (int i = 0; i < nn; i++)
            {
                double ux = _result.Displacements[i * 3];
                double uy = _result.Displacements[i * 3 + 1];
                double uz = _result.Displacements[i * 3 + 2];
                mags[i] = Math.Sqrt(ux * ux + uy * uy + uz * uz);
            }
            // Map to edges by averaging node values
            int ne = _result.Network.Graph.Ne;
            var edgeMags = new double[ne];
            for (int e = 0; e < ne; e++)
            {
                var edge = _result.Network.Graph.Edges[e];
                edgeMags[e] = (mags[edge.Start.Index] + mags[edge.End.Index]) * 0.5;
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
                if (_result?.DeformedNodes == null) return BoundingBox.Empty;
                var bb = new BoundingBox(_result.DeformedNodes);
                bb.Inflate(1.0);
                return bb;
            }
        }

        protected override Bitmap Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-000000000021");
    }
}
