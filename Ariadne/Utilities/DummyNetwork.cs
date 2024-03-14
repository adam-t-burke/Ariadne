using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Ariadne.FDM;
using Grasshopper.Kernel.Data;
using System.Drawing;

namespace Ariadne.Utilities
{
    public class DummyNetwork : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the DummyNetwork class.
        /// </summary>
        public DummyNetwork()
          : base("DummyNetwork", "Dummy",
              "Dummy network to test downstream functionality",
              "Ariadne", "Utilities")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("FDM Network", "Network", "FDM Network", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Point3d> anchors = new List<Point3d>
            {
                new Point3d(0, 0, 0),
                new Point3d(1, 0, 0),
                new Point3d(0, 1, 0),
                new Point3d(1, 1, 0)
            };

            List<Curve> edges = new List<Curve>
            {
                new LineCurve(new Point3d(0, 0, 0), new Point3d(0.25, 0.25, 0)),
                new LineCurve(new Point3d(0.25, 0.75, 0), new Point3d(0, 1, 0)),
                new LineCurve(new Point3d(0.75, 0.75, 0), new Point3d(1, 1, 0)),
                new LineCurve(new Point3d(0.75, 0.25, 0), new Point3d(1, 0, 0)),
                new LineCurve(new Point3d(0.25, 0.25, 0), new Point3d(0.25, 0.75, 0)),
                new LineCurve(new Point3d(0.25, 0.75, 0), new Point3d(0.75, 0.75, 0)),
                new LineCurve(new Point3d(0.75, 0.75, 0), new Point3d(0.75, 0.25, 0)),
                new LineCurve(new Point3d(0.75, 0.25, 0), new Point3d(0.25, 0.25, 0))
            };

            double etol = 0.01;
            double atol = 0.01;

            FDM_Network Network = new FDM_Network(edges, etol, anchors, atol);

            DA.SetData(0, Network);


        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                return Properties.Resources.Dummy_Network;

            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("36E2CE42-D573-4C81-A226-94ED15B158F0"); }
        }
    }
}