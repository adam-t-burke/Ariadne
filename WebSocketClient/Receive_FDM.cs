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
using System.Threading.Tasks;
using System.Collections.Concurrent;
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

            // Create node index mapping for O(1) lookups (eliminates expensive IndexOf calls)
            var nodeIndexMap = new Dictionary<Node, int>(oldNetwork.Graph.Nn);
            for (int i = 0; i < oldNetwork.Graph.Nn; i++)
            {
                nodeIndexMap[oldNetwork.Graph.Nodes[i]] = i;
            }

            // Phase 1: Parallel node creation with pre-allocated capacity
            newGraph.Nodes = new List<Node>(oldNetwork.Graph.Nn);
            for (int i = 0; i < oldNetwork.Graph.Nn; i++)
            {
                newGraph.Nodes.Add(null); // Pre-allocate space
            }

            Parallel.For(0, oldNetwork.Graph.Nn, i =>
            {
                Node node = new()
                {
                    Anchor = oldNetwork.Graph.Nodes[i].Anchor,
                    Value = new Point3d(response.X[i], response.Y[i], response.Z[i])
                };
                newGraph.Nodes[i] = node;
            });

            // Phase 2: Parallel neighbor reconstruction using index mapping
            Parallel.For(0, oldNetwork.Graph.Nn, i =>
            {
                var neighbors = new List<Node>(oldNetwork.Graph.Nodes[i].Neighbors.Count);
                foreach (Node neighbor in oldNetwork.Graph.Nodes[i].Neighbors)
                {
                    int neighborIndex = nodeIndexMap[neighbor];
                    neighbors.Add(newGraph.Nodes[neighborIndex]);
                }
                newGraph.Nodes[i].Neighbors = neighbors;
            });

            // Phase 3: Parallel edge creation using index mapping
            newGraph.Edges = new List<Edge>(oldNetwork.Graph.Ne);
            for (int i = 0; i < oldNetwork.Graph.Ne; i++)
            {
                newGraph.Edges.Add(null); // Pre-allocate space
            }

            Parallel.For(0, oldNetwork.Graph.Ne, i =>
            {
                int startIndex = nodeIndexMap[oldNetwork.Graph.Edges[i].Start];
                int endIndex = nodeIndexMap[oldNetwork.Graph.Edges[i].End];

                Edge edge = new()
                {
                    Start = newGraph.Nodes[startIndex],
                    End = newGraph.Nodes[endIndex],
                    Q = response.Q[i],
                    ReferenceID = oldNetwork.Graph.Edges[i].ReferenceID
                };
                edge.Value = new LineCurve(edge.Start.Value, edge.End.Value);

                newGraph.Edges[i] = edge;
            });

            // Optimized FDM_Network construction with pre-allocated collections
            FDM_Network newNetwork = new() 
            {
                Graph = newGraph,
                Valid = true,
                FreeNodes = oldNetwork.FreeNodes,
                FixedNodes = oldNetwork.FixedNodes,
                ATol = oldNetwork.ATol,
                ETol = oldNetwork.ETol
            };

            // Replace LINQ with direct array allocation for better performance
            newNetwork.Free = new List<Node>(oldNetwork.FreeNodes.Count);
            for (int i = 0; i < oldNetwork.FreeNodes.Count; i++)
            {
                newNetwork.Free.Add(newGraph.Nodes[oldNetwork.FreeNodes[i]]);
            }

            newNetwork.Fixed = new List<Node>(oldNetwork.FixedNodes.Count);
            newNetwork.Anchors = new List<Point3d>(oldNetwork.FixedNodes.Count);
            for (int i = 0; i < oldNetwork.FixedNodes.Count; i++)
            {
                var fixedNode = newGraph.Nodes[oldNetwork.FixedNodes[i]];
                newNetwork.Fixed.Add(fixedNode);
                newNetwork.Anchors.Add(fixedNode.Value);
            }




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
                // Early exit for empty traces
                if (NodeTrace == null || NodeTrace.Count == 0 || X == null || X.Count == 0)
                    return new GH_Structure<GH_Point>();

                GH_Structure<GH_Point> Trace = new GH_Structure<GH_Point>();

                var iterations = NodeTrace.Count;
                var pointsPerIteration = X.Count;
                
                // Pre-allocate arrays to avoid ConcurrentBag overhead
                var allPoints = new GH_Point[iterations][];
                
                // Parallel processing with pre-allocated arrays
                Parallel.For(0, iterations, iter =>
                {
                    var iterData = NodeTrace[iter];
                    var xData = iterData[0];
                    var yData = iterData[1]; 
                    var zData = iterData[2];
                    
                    // Pre-allocate array for this iteration
                    allPoints[iter] = new GH_Point[pointsPerIteration];
                    
                    // Process points for this iteration in parallel
                    Parallel.For(0, pointsPerIteration, i =>
                    {
                        Point3d pt = new(xData[i], yData[i], zData[i]);
                        allPoints[iter][i] = new GH_Point(pt);
                    });
                });

                // Sequential append to maintain tree structure integrity (much faster than grouping)
                for (int iter = 0; iter < iterations; iter++)
                {
                    var path = new GH_Path(iter);
                    for (int i = 0; i < pointsPerIteration; i++)
                    {
                        Trace.Append(allPoints[iter][i], path);
                    }
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