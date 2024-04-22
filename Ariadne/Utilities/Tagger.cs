using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Xml.Linq;
using Rhino.Collections;
using Ariadne.FDM;
using Ariadne.Graphs;
using System.Drawing;

namespace Ariadne.Utilities
{
    public class Tagger : GH_Component
    {
        FDM_Problem network;
        List<string> tags;
        List<Point3d> points;
        bool show;
        double size;
        System.Drawing.Color col;

        /// <summary>
        /// Initializes a new instance of the Tagger class.
        /// </summary>
        public Tagger()
          : base("ElementTagger", "Tagger",
              "Tagged information about elements",
              "Ariadne", "Utilities")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("FDM Network", "Network", "FDM network class", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Show", "Show", "Show tags", GH_ParamAccess.item, true);
            pManager.AddNumberParameter("Tag Size", "Size", "Size of text tags", GH_ParamAccess.item, 2);
            pManager.AddColourParameter("Tag Colour", "Colour", "Colour of text tags", GH_ParamAccess.item, System.Drawing.Color.Black);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            network = new FDM_Problem();
            show = true;
            size = 2;
            col = System.Drawing.Color.Black;

            if (!DA.GetData(0, ref network)) return;
            DA.GetData(1, ref show);
            DA.GetData(2, ref size);
            DA.GetData(3, ref col);

            List<double> lengths = UtilityFunctions.GetLengths(network.Network.Graph.Edges);
            List<double> forces = UtilityFunctions.GetForces(lengths, network.Q);

            GetPoints(network);
            GetText(network, forces);


        }

        private void GetPoints(FDM_Problem network)
        {
            points = new List<Point3d>();

            foreach (Edge edge in network.Network.Graph.Edges)
            {
                Curve curve = edge.Value;
                points.Add(curve.PointAtNormalizedLength(0.5));
            }
        }

        private void GetText(FDM_Problem network, List<double> forces)
        {
            tags = new List<string>();

            for (int i = 0; i < network.Network.Graph.Ne; i++)
            {
                Curve curve = network.Network.Graph.Edges[i].Value;
                string length = string.Format("{0:0.0}", curve.GetLength());
                string density = string.Format("{0:0.0}", network.Q[i]);
                string force = string.Format("{0:0.00}", forces[i]);

                string info = $@"L = {length}
                                 q = {density}
                                 F = {force}";

                tags.Add(info);
            }
        }

        public override BoundingBox ClippingBox
        {
            get
            {
                Point3dList pts = new Point3dList(points);
                return pts.BoundingBox;
            }
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            base.DrawViewportWires(args);

            Plane plane;
            args.Viewport.GetFrustumFarPlane(out plane);

            if (show)
            {
                for (int i = 0; i < tags.Count; i++)
                {
                    string text = tags[i];
                    Point3d point = points[i];
                    plane.Origin = point;

                    Rhino.Display.Text3d drawText = new Rhino.Display.Text3d(text, plane, size);
                    args.Display.Draw3dText(drawText, col);
                    drawText.Dispose();
                }
            }
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override Bitmap Icon => Properties.Resources.tagger;

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("69EEC54F-B955-4840-8A46-4E1560298BE7"); }
        }
    }
}