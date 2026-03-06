using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using Rhino.Collections;
using System;
using System.Drawing;
using System.Collections.Generic;
using Ariadne.Utilities;
using Ariadne.Graphs;

namespace Ariadne.FDM
{
    /// <summary>
    /// Builds an FDM network from a tree of edges and anchor points. Outputs a single Network for use with Theseus Solve.
    /// </summary>
    public class ConstructNetwork : GH_Component
    {
        FDM_Network oldNet = null!;
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
            pManager.AddCurveParameter("Edges", "E", "Edges between nodes", GH_ParamAccess.tree);
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
            DA.DisableGapLogic();

            if (oldNet != null)
            {
                if(oldNet.IsUpdating)
                {
                    return;
                }
            }
            
            //initialize
            GH_Structure<GH_Curve> edges = new GH_Structure<GH_Curve>();
            double e_tol = 0.001;
            List<Point3d> anchors = new List<Point3d>();
            double a_tol = 0.001;


            //assign
            if (!DA.GetDataTree(0, out edges)) return;
            if (!DA.GetData(1, ref e_tol)) return;
            if (!DA.GetDataList(2, anchors)) return;
            if (!DA.GetData(3, ref a_tol)) return;

            //create graph from edge tree structure
            Graph graph = new(edges, e_tol);
            //create network from graph and anchors
            FDM_Network network = new(graph, anchors, a_tol);

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


            oldNet = network;
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