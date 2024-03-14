using Grasshopper;
using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino.Collections;
using System;
using System.Drawing;
using System.Collections.Generic;
using Ariadne.Utilities;
using Grasshopper.Kernel.Data;

namespace Ariadne.FDM
{
    public class ConstructNetwork : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public ConstructNetwork()
          : base("Construct FDM Network", "FDM Network",
            "Create FDM network",
            "Ariadne", "Design")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Edges", "E", "Edges between nodes", GH_ParamAccess.list);
            pManager.AddNumberParameter("Edge Tolerence", "E Tol", "Geometric tolerance for connecting edges", GH_ParamAccess.item, 0.001);
            pManager.AddPointParameter("Anchors", "A", "Anchor points", GH_ParamAccess.list);
            pManager.AddNumberParameter("Anchor Tolerance", "A Tol", "Geometric tolerance for picking anchor nodes geometry", GH_ParamAccess.item, 0.001);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Network", "Network", "FDM network", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //initialize
            List<Curve> edges = new List<Curve>();
            double e_tol = 0.001;
            List<Point3d> anchors = new List<Point3d>();
            double a_tol = 0.001;


            //assign
            if (!DA.GetDataList(0, edges)) return;
            if (!DA.GetData(1, ref e_tol)) return;
            if (!DA.GetDataList(2, anchors)) return;
            if (!DA.GetData(3, ref a_tol)) return;

            //create network from graph
            FDM_Network network = new FDM_Network(edges, e_tol, anchors, a_tol);

            //Check sufficent anchor definitions
            if (!network.AnchorCheck())
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "For stability, define at least 2 anchor points.");
                return;
            }

            //check that fixed/free nodes were properly captured
            if (!network.NFCheck())
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "The number of fixed points does not match the number of input anchor points. Check Anchor Tolerance.");
                return;
            }

            //assign
            DA.SetData(0, network);

        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// You can add image files to your project resources and access them like this:
        /// return Resources.IconForThisComponent;
        /// </summary>
        protected override Bitmap Icon => Properties.Resources.Construct_Network;

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid => new Guid("5DBF67D7-9DCE-4B0B-9DE9-15422CD2B4FA");
    }
}