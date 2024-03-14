using System;
using System.Collections.Generic;
using System.Xml.Serialization;
using Ariadne.Utilities;
using Eto.Forms;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Rhino.Collections;
using Rhino.Geometry;
using System.Linq;
using System.Drawing;

namespace Ariadne.FDM
{
    public class DeconstructNetwork : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the DeconstructFDMNetwork class.
        /// </summary>
        public DeconstructNetwork()
          : base("Deconstruct FDM Network", "FDM Network Info",
              "Deconstructed attributes of a FDM Network",
              "Ariadne", "Design")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("FDM Network", "FDM Network", "FDM Network to deconstruct", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Nodes", "Nodes", "Nodes in the FDM Network", GH_ParamAccess.list);
            pManager.AddGenericParameter("Edges", "Edges", "Edges in the FDM Network", GH_ParamAccess.list);
            pManager.AddNumberParameter("Edge Indices", "Edge Indices", "Edges defined by start and end node indices", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Adjacency List", "Adjacency List", "Adjacency List representation of graph", GH_ParamAccess.list);
            pManager.AddNumberParameter("Free Node Indices", "Free Indices", "Indices of free nodes", GH_ParamAccess.list);
            pManager.AddGenericParameter("Free Nodes", "Free Nodes", "Free nodes in the FDM Network", GH_ParamAccess.list);
            pManager.AddNumberParameter("Fixed Node Indices", "Fixed Indices", "Indices of fixed nodes", GH_ParamAccess.list);
            pManager.AddGenericParameter("Fixed Nodes", "Fixed Nodes", "Fixed nodes in the FDM Network", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            FDM_Network network = new FDM_Network();
            if (!DA.GetData(0, ref network)) return;

            network.Graph.EdgeIndicesToTree();
            network.Graph.AdjacencyListToTree();

            IEnumerable<int> freeIndices = Enumerable.Range(0, network.Graph.Nn - network.Anchors.Count);
            IEnumerable<int> fixedIndices = Enumerable.Range(network.Graph.Nn - network.Anchors.Count, network.Graph.Nn);

            DA.SetDataList(0, network.Graph.Nodes.Select(node => node.Position));
            DA.SetDataList(1, network.Graph.Edges.Select(edge => edge.Curve));
            DA.SetDataTree(2, network.Graph.IndicesTree);
            DA.SetDataTree(3, network.Graph.AdjacencyTree);
            DA.SetDataList(4, freeIndices);
            DA.SetDataList(5, network.Free);
            DA.SetDataList(6, fixedIndices);
            DA.SetDataList(7, network.Fixed);
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                return Properties.Resources.Deconstruct_Network;

            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("A94E5464-C2E0-42D7-B6A4-E65BBF37CC09"); }
        }
    }
}