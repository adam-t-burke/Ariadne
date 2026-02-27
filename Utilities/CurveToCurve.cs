using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Ariadne.FDM;
using System.Drawing;

namespace Ariadne.Utilities
{
    /// <summary>
    /// Draws lines connecting initial and solved edge geometry to show deformation or mapping.
    /// </summary>
    public class CurveToCurve : GH_Component
    {
        Line[] lines = null!;
        System.Drawing.Color col;
        int thickness;
        bool show;
        int index;

        /// <summary>
        /// Initializes a new instance of the CurveToCurve class.
        /// </summary>
        public CurveToCurve()
          : base("CurvetoCurve", "CurvePairs",
              "Connecting the starting Value of edges to equilibrium Values",
              "Ariadne", "Utilities")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Initial Network", "Network1", "Network before analysis", GH_ParamAccess.item);
            pManager.AddGenericParameter("Solved Network", "Network2", "Network after analysis", GH_ParamAccess.item);
            pManager.AddColourParameter("Line Colour", "Colour", "Colour of line", GH_ParamAccess.item, System.Drawing.Color.SlateGray);
            pManager.AddIntegerParameter("Line Thickness", "Weight", "Weight of lines", GH_ParamAccess.item, 2);
            pManager.AddBooleanParameter("Show", "Show", "Show assignments", GH_ParamAccess.item, true);
            pManager.AddIntegerParameter("Indexer", "Indexer", "Specific curve index", GH_ParamAccess.item, -1);
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
            FDM_Network fdm1 = new();
            FDM_Network fdm2 = new();
            col = System.Drawing.Color.SlateGray;
            thickness = 2;
            show = true;

            if (!DA.GetData(0, ref fdm1)) return;
            if (!DA.GetData(1, ref fdm2)) return;

            if (fdm1.Graph.Ne != fdm2.Graph.Ne) AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Input and output networks must have same number of elements");

            lines = GetLines(fdm1, fdm2);

            DA.GetData(2, ref col);
            DA.GetData(3, ref thickness);
            DA.GetData(4, ref show);
            DA.GetData(5, ref index);

            if (index >= fdm1.Graph.Ne)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Index greater than total number of elements");
            }
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            base.DrawViewportWires(args);

            if (show)
            {
                if (index == -1)
                {
                    args.Display.DrawLines(lines, col, thickness);
                }
                else
                {
                    args.Display.DrawLine(lines[index], col, thickness);
                }
            }

        }

        private Line[] GetLines(FDM_Network n1, FDM_Network n2)
        {
            Line[] ls = new Line[n1.Graph.Edges.Count];

            for (int i = 0; i < n1.Graph.Edges.Count; i++)
            {
                Curve c1 = n1.Graph.Edges[i].Value;
                Curve c2 = n2.Graph.Edges[i].Value;

                ls[i] = new Line(c1.PointAtNormalizedLength(0.5), c2.PointAtNormalizedLength(0.5));
            }

            return ls;
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override Bitmap Icon => Properties.Resources.curve2curve;

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("4E9B60C2-E43A-4680-800A-B66871229F3B"); }
        }
    }
}