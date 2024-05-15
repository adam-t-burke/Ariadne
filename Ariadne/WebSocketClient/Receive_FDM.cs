using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Grasshopper;
using Rhino;
using Rhino.Geometry;
using System.Drawing;
using Ariadne.FDM;
using System.Text.Json;
using System.Text.Json.Serialization;
using System.Text.Json.Nodes;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using System.Runtime.CompilerServices;
using Ariadne.Graphs;
using System.Linq;
using static System.Runtime.InteropServices.JavaScript.JSType;

namespace Ariadne.WebSocketClient
{
    public class ReceiveFDM : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Receive class.
        /// </summary>
        public ReceiveFDM()
          : base("Receive FDM", "Receive FDM",
              "Receive response from Julia FDM Optimization server.",
              "Ariadne", "Communication")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddTextParameter("Response", "Response", "Response from server.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Network", "Network", "FDM Network", GH_ParamAccess.item);
            //pManager.AddBooleanParameter("Update", "Update", "Update output constantly, not recommended if you are sliding over lots of values, or making adjustments", GH_ParamAccess.item, true);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Nodes", "Nodes", "Nodes in the FDM Network", GH_ParamAccess.list);
            pManager.AddCurveParameter("Edges", "Edges", "Edges in the FDM Network", GH_ParamAccess.list);
            pManager.AddTextParameter("Iterations", "Iter", "Iterations of form finding.", GH_ParamAccess.item);
            pManager.AddNumberParameter("Force Densities", "Q", "Final force density values", GH_ParamAccess.list);
            pManager.AddNumberParameter("Loss", "Loss", "Optimization final loss value", GH_ParamAccess.item);
            pManager.AddNumberParameter("Loss trace", "Loss trace", "Optimization loss over time.", GH_ParamAccess.list);
            pManager.AddPointParameter("Node trace", "Node trace ", "Node Values over time.", GH_ParamAccess.tree);
            pManager.AddGenericParameter("Updated Network", "Network", "New Network with updated values", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            string _response = "";
            FDM_Network oldNetwork = new();
            Graph newGraph = new();           

            if (!DA.GetData(0, ref _response)) return;
            if (!DA.GetData(1, ref oldNetwork)) return;

            newGraph.Tolerance = oldNetwork.Graph.Tolerance;
            newGraph.Nodes = new();
            newGraph.Edges = new();
            newGraph.IndicesTree = oldNetwork.Graph.IndicesTree;
            newGraph.AdjacencyTree = oldNetwork.Graph.AdjacencyTree;
            
            FDM_Response response = JsonSerializer.Deserialize<FDM_Response>(_response);

            if (response.Q.Count != oldNetwork.Graph.Ne)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "Number of edges in network does not match number of force densities in response. It is safe to ignore this warning if you just changed the topology of your network");
                return;
            }

            if (oldNetwork.Graph.Nn != response.X.Count)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "Number of nodes in network does not match number of nodes in response. It is safe to ignore this warning if you just changed the topology of your network");
                return;
            }

            for (int i = 0; i < oldNetwork.Graph.Nn; i++)
            {
                Node node = new()
                {
                    Anchor = oldNetwork.Graph.Nodes[i].Anchor,
                    Value = new Point3d(response.X[i], response.Y[i], response.Z[i])
                }; 
                
                newGraph.Nodes.Add(node);
            }

            for (int i = 0; i < oldNetwork.Graph.Nn; i++)
            {
                foreach(Node neighbor in oldNetwork.Graph.Nodes[i].Neighbors)
                {
                    newGraph.Nodes[i].Neighbors.Add(newGraph.Nodes[oldNetwork.Graph.Nodes.IndexOf(neighbor)]);
                }
            }

            for (int i = 0; i < oldNetwork.Graph.Ne; i++)
            {
                Edge edge = new()
                {
                    Start = newGraph.Nodes[oldNetwork.Graph.Nodes.IndexOf(oldNetwork.Graph.Edges[i].Start)],
                    End = newGraph.Nodes[oldNetwork.Graph.Nodes.IndexOf(oldNetwork.Graph.Edges[i].End)],
                    Q = response.Q[i]
                };
                edge.Value = new LineCurve(edge.Start.Value, edge.End.Value);
                edge.ReferenceID = oldNetwork.Graph.Edges[i].ReferenceID;

                newGraph.Edges.Add(edge);                
            }

            FDM_Network newNetwork = new() {  };
            newNetwork.Graph = newGraph;
            newNetwork.Valid = true;
            newNetwork.FreeNodes = oldNetwork.FreeNodes;
            newNetwork.FixedNodes = oldNetwork.FixedNodes;
            newNetwork.Free = oldNetwork.FreeNodes.Select(x => newGraph.Nodes[x]).ToList();
            newNetwork.Fixed = oldNetwork.FixedNodes.Select(x => newGraph.Nodes[x]).ToList();
            newNetwork.Anchors = newNetwork.Fixed.Select(x => x.Value).ToList();
            newNetwork.ATol = oldNetwork.ATol;
            newNetwork.ETol = oldNetwork.ETol;




            DA.SetDataList(0, newNetwork.Graph.Nodes);
            DA.SetDataList(1, newNetwork.Graph.Edges);
            DA.SetData(2, response.Iter);
            DA.SetDataList(3, response.Q);
            DA.SetData(4, response.Loss);
            DA.SetDataList(5, response.Losstrace);
            if (response.NodeTrace != null)
            {
                DA.SetDataTree(6, response.NodeTraceToTree());
            }
            DA.SetData(7, newNetwork);


        }

        internal class FDM_Response
        {
            public int Iter { get; set; }
            public double Loss { get; set; }
            public List<double> Q { get; set; }
            public List<double> X { get; set; }
            public List<double> Y { get; set; }
            public List<double> Z { get; set; }
            public List<double> Losstrace { get; set; }
            public List<List<List<double>>> NodeTrace { get; set; }

            /// <summary>
            /// Convert node trace to grasshopper tree structure w/points
            /// </summary>
            /// <returns></returns>
            public GH_Structure<GH_Point> NodeTraceToTree()
            {
                int counter = 0;
                GH_Structure<GH_Point> Trace = new GH_Structure<GH_Point>();

                foreach (List<List<double>> iter in NodeTrace)
                {
                    for (int i = 0; i < X.Count; i++)
                    {
                        Point3d pt = new(iter[0][i], iter[1][i], iter[2][i]);
                        GH_Point gh_Pt = new(pt);
                        Trace.Append(gh_Pt, new(counter));
                    }
                    counter++;
                }

                return Trace;

            }
        }


        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                return Properties.Resources.ReceiveFDM;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("A3A719C1-BCA9-4222-9ADC-25EEEE807256"); }
        }
    }
}