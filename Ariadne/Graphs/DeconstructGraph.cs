using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Ariadne.Graphs;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Ariadne.Graphs
{
    public class DeconstructGraph : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the MyComponent1 class.
        /// </summary>
        public DeconstructGraph()
          : base("DeconstructGraph", "Deconstruct Graph",
            "Deconstruct Curve Graph",
            "Ariadne", "Graphs")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Graph", "G", "Graph representation of curve network", GH_ParamAccess.item);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Nodes", "Nodes", "Ordered list of nodes in graph", GH_ParamAccess.list);
            pManager.AddCurveParameter("Edge Curves", "Edges", "Curve for each graph edge", GH_ParamAccess.list);
            pManager.AddNumberParameter("Edge Indices", "Edge Indices", "Edges defined by start and end node indices", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Adjacency List", "Adjacency List", "Adjacency List representation of graph", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            Graph curveGraph = new Graph();
            if (!DA.GetData(0, ref curveGraph)) return;

            curveGraph.EdgeIndicesToTree();
            curveGraph.AdjacencyListToTree();

            DA.SetDataList(0, curveGraph.Nodes.Select(node => node.Position));
            DA.SetDataList(1, curveGraph.Edges.Select(edge => edge.Curve));
            DA.SetDataTree(2, curveGraph.IndicesTree);
            DA.SetDataTree(3, curveGraph.AdjacencyTree);

        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                return Properties.Resources.Deconstruct_Graph;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("66A08FD9-D549-43EA-B3BE-83FBFA1898BF"); }
        }
    }
}