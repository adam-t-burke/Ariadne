using System;
using System.Collections.Generic;
using System.Linq;
using Grasshopper.Kernel;
using Rhino.Geometry;
using Grasshopper.Kernel.Parameters;
using Grasshopper.GUI.Gradient;
using Ariadne.FDM;
using Ariadne.Graphs;
using System.Drawing;

namespace Ariadne.Utilities
{
    /// <summary>
    /// Visualizes an FDM network with optional property coloring (e.g. force, length), loads, and reactions.
    /// </summary>
    internal class Visualize : GH_Component
    {
        private FDM_Network _network;
        private Line[] _edges;
        private List<double> _property;
        private Line[] _externalForces;
        private Line[] _reactionForces;
        private Color _c0;
        private Color _cMed;
        private Color _c1;
        private GH_Gradient _grad;
        private int _thickness;
        private Color _cLoad;
        private Color _cReact;
        private bool _showLoads;
        private bool _showReactions;
        private int _prop;
        private bool _show;

        private static readonly Color DefaultLightGray = Color.FromArgb(230, 231, 232);
        private static readonly Color DefaultBlue = Color.FromArgb(62, 168, 222);
        private static readonly Color DefaultPink = Color.FromArgb(255, 123, 172);
        private static readonly Color DefaultGreen = Color.FromArgb(71, 181, 116);
        private static readonly Color DefaultRed = Color.FromArgb(150, 235, 52, 73);

        public Visualize()
          : base("VisualizeNetwork", "Visualize",
              "Visualize a FDM network",
              "Ariadne", "Utilities")
        {
        }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddBooleanParameter("Show", "Show", "Active status of component", GH_ParamAccess.item, true);
            pManager.AddGenericParameter("Network", "Network", "Network to visualize", GH_ParamAccess.item);
            pManager.AddVectorParameter("Loads", "P", "Applied loads", GH_ParamAccess.list, new Vector3d(0, 0, 0));
            pManager.AddNumberParameter("Load Scale", "Pscale", "Scale factor for length of arrows", GH_ParamAccess.item, 100);
            pManager.AddColourParameter("ColourMin", "Cmin", "Colour for minimum value", GH_ParamAccess.item, DefaultPink);
            pManager.AddColourParameter("ColourMed", "Cmed", "Colour for neutral value", GH_ParamAccess.item, DefaultLightGray);
            pManager.AddColourParameter("ColourMax", "Cmax", "Colour for maximum value", GH_ParamAccess.item, DefaultBlue);
            pManager.AddIntegerParameter("Color Property", "Property", "Property displayed by colour gradient", GH_ParamAccess.item, 0);
            pManager.AddIntegerParameter("Line Thickness", "Thickness", "Thickness of preview lines", GH_ParamAccess.item, 8);
            pManager.AddColourParameter("Load Colour", "Cload", "Colour for applied loads", GH_ParamAccess.item, DefaultRed);
            pManager.AddBooleanParameter("Show Loads", "Load", "Show external loads in preview", GH_ParamAccess.item, true);
            pManager.AddColourParameter("Reaction Colour", "Creaction", "Colour for support reactions", GH_ParamAccess.item, DefaultGreen);
            pManager.AddBooleanParameter("Show Reactions", "Reaction", "Show anchor reactions in preview", GH_ParamAccess.item, true);

            Param_Integer param = pManager[7] as Param_Integer;
            param.AddNamedValue("None", -1);
            param.AddNamedValue("Force", 0);
            param.AddNamedValue("Q", 1);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            ClearData();

            FDM_Network network = new();
            List<Vector3d> loads = new();
            double scale = 1.0;
            _c0 = DefaultPink;
            _cMed = DefaultLightGray;
            _c1 = DefaultBlue;
            _cLoad = DefaultPink;
            _showLoads = true;
            _cReact = DefaultGreen;
            _showReactions = false;

            DA.GetData(0, ref _show);
            if (!DA.GetData(1, ref network)) return;
            DA.GetDataList(2, loads);
            DA.GetData(3, ref scale);
            DA.GetData(4, ref _c0);
            DA.GetData(5, ref _cMed);
            DA.GetData(6, ref _c1);
            DA.GetData(7, ref _prop);
            DA.GetData(8, ref _thickness);
            DA.GetData(9, ref _cLoad);
            DA.GetData(10, ref _showLoads);
            DA.GetData(11, ref _cReact);
            DA.GetData(12, ref _showReactions);

            _network = network;

            _edges = EdgesToLines(network.Graph.Edges);
            _externalForces = LoadMaker(network.Graph.Nodes, network.FreeNodes, loads, scale);

            List<double> lengths = UtilityFunctions.GetLengths(network.Graph.Edges);
            List<double> forces = UtilityFunctions.GetForces(lengths, network.Graph.Edges);

            if (_showReactions)
            {
                List<Vector3d> reactions = UtilityFunctions.GetReactions(forces, network);
                List<Point3d> nodePoints = network.Graph.Nodes.Select(n => n.Value).ToList();
                _reactionForces = ReactionMaker(reactions, nodePoints, network.FixedNodes, scale);
            }
            else
            {
                _reactionForces = Array.Empty<Line>();
            }

            if (_prop == 0)
            {
                _property = forces;
                GradientMaker(_property);
            }
            else if (_prop == 1)
            {
                _property = network.Graph.Edges.Select(x => x.Q).ToList();
                GradientMaker(_property);
            }
        }

        public override BoundingBox ClippingBox
        {
            get
            {
                if (_network == null || _edges == null || _edges.Length == 0)
                    return BoundingBox.Empty;

                BoundingBox bb = new BoundingBox(_network.Graph.Nodes.Select(node => node.Value));
                if (_externalForces != null)
                    for (int i = 0; i < _externalForces.Length; i++)
                        bb.Union(_externalForces[i].BoundingBox);
                if (_reactionForces != null)
                    for (int i = 0; i < _reactionForces.Length; i++)
                        bb.Union(_reactionForces[i].BoundingBox);

                return bb;
            }
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            base.DrawViewportWires(args);

            if (!_show || _edges == null) return;

            if (_showLoads && _externalForces != null && _externalForces.Length > 0)
                args.Display.DrawArrows(_externalForces, _cLoad);

            if (_showReactions && _reactionForces != null && _reactionForces.Length > 0)
                args.Display.DrawArrows(_reactionForces, _cReact);

            if (_prop == -1)
            {
                args.Display.DrawLines(_edges, _c1, _thickness);
            }
            else if (_property != null && _grad != null)
            {
                for (int i = 0; i < _edges.Length; i++)
                {
                    args.Display.DrawLine(_edges[i], _grad.ColourAt(_property[i]), _thickness);
                }
            }
        }

        private static Line[] EdgesToLines(List<Edge> edges)
        {
            Line[] lines = new Line[edges.Count];
            for (int i = 0; i < edges.Count; i++)
            {
                lines[i] = new Line(edges[i].Start.Value, edges[i].End.Value);
            }
            return lines;
        }

        private void GradientMaker(List<double> property)
        {
            double minprop = property.Min();
            double maxprop = property.Max();

            int signmin = Math.Sign(minprop);
            int signmax = Math.Sign(maxprop);

            if (signmin <= 0 && signmax <= 0)
            {
                _grad = new GH_Gradient();
                _grad.AddGrip(minprop, _c0);
                _grad.AddGrip(0, _cMed);
            }
            else if (signmin <= 0 && signmax >= 0)
            {
                _grad = new GH_Gradient();
                _grad.AddGrip(minprop, _c0);
                _grad.AddGrip(0, _cMed);
                _grad.AddGrip(maxprop, _c1);
            }
            else
            {
                _grad = new GH_Gradient();
                _grad.AddGrip(0, _cMed);
                _grad.AddGrip(maxprop, _c1);
            }
        }

        private static Line[] ReactionMaker(List<Vector3d> anchorforces, List<Point3d> points, List<int> fixedIndices, double scale)
        {
            if (anchorforces.Count == 0) return Array.Empty<Line>();

            var mags = anchorforces.Select(p => p.Length).ToList();
            var normalizer = mags.Max();
            if (normalizer < 1e-12) return Array.Empty<Line>();

            List<Line> reactions = new(fixedIndices.Count);
            for (int i = 0; i < fixedIndices.Count; i++)
            {
                var index = fixedIndices[i];
                reactions.Add(new Line(points[index], anchorforces[i] * 3 * scale / normalizer));
            }

            return reactions.ToArray();
        }

        private static Line[] LoadMaker(List<Node> nodes, List<int> N, List<Vector3d> loads, double scale)
        {
            List<Line> loadvectors = new();

            if (N.Count != loads.Count && loads.Count != 1)
                throw new ArgumentException("Length of force vectors must be 1 or match length of free nodes.");

            if (loads.Count == 1)
            {
                double len = loads[0].Length;
                if (len < 1e-12) return Array.Empty<Line>();

                for (int i = 0; i < N.Count; i++)
                {
                    Point3d p = nodes[N[i]].Value;
                    loadvectors.Add(new Line(p, loads[0] / len * scale));
                }
            }
            else
            {
                var lns = loads.Select(p => p.Length).ToList();
                var normalizer = lns.Max();
                if (normalizer < 1e-12) return Array.Empty<Line>();

                for (int i = 0; i < N.Count; i++)
                {
                    Point3d p = nodes[N[i]].Value;
                    Vector3d l = loads[i];

                    if (l.Length < 0.1) continue;

                    loadvectors.Add(new Line(p, l * scale / normalizer));
                }
            }

            return loadvectors.ToArray();
        }

        protected override Bitmap Icon => Properties.Resources.visualize;

        public override Guid ComponentGuid => new Guid("5B9176B3-B940-4C2C-AFFE-BF4532FB2111");
    }
}
