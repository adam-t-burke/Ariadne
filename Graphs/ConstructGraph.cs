using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;
using Rhino.Collections;
using System;
using System.Collections.Generic;
using Ariadne.Graphs;
using System.Drawing;

namespace Ariadne.Graphs
{
    public class ConstructGraph : GH_Component
    {
        /// <summary>
        /// Each implementation of GH_Component must provide a public 
        /// constructor without any arguments.
        /// Category represents the Tab in which the component will appear, 
        /// Subcategory the panel. If you use non-existing tab or panel names, 
        /// new tabs/panels will automatically be created.
        /// </summary>
        public ConstructGraph()
          : base("ConstructGraph", "Construct Graph",
            "Create Graph from Curve Network",
            "Ariadne", "Graphs")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Curves", "C", "Curves to convert to graph", GH_ParamAccess.tree);
            pManager.AddNumberParameter("IntersectionTolerance", "tol", "Geometric tolerance for connecting geometry", GH_ParamAccess.item, 0.0001);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Graph", "Graph", "Graph representation of curve network", GH_ParamAccess.item);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object can be used to retrieve data from input parameters and 
        /// to store data in output parameters.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<GH_Curve> edgeTree = new();
            double tol = 1.0;

            if (!DA.GetDataTree(0, out edgeTree)) return;
            DA.GetData(1, ref tol);

            Graph curveGraph = new Graph(edgeTree, tol);
   
            DA.SetData(0, curveGraph);
        }

        /// <summary>
        /// Provides an Icon for every component that will be visible in the User Interface.
        /// Icons need to be 24x24 pixels.
        /// You can add image files to your project resources and access them like this:
        /// return Resources.IconForThisComponent;
        /// </summary>
        protected override Bitmap Icon => Properties.Resources.Construct_Graph;

        /// <summary>
        /// Each component must have a unique Guid to identify it. 
        /// It is vital this Guid doesn't change otherwise old ghx files 
        /// that use the old ID will partially fail during loading.
        /// </summary>
        public override Guid ComponentGuid => new Guid("3b2c9fff-0e90-4059-b65a-1668299d121d");
    }
}