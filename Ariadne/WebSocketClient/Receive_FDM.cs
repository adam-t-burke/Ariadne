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
            pManager.AddPointParameter("Node trace", "Node trace ", "Node positions over time.", GH_ParamAccess.tree);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            string _response = "";
            FDM_Network network = new FDM_Network();
            List<Point3d> nodes = new List<Point3d>();
            List<Curve> edges = new List<Curve>();
           

            if (!DA.GetData(0, ref _response)) return;
            if (!DA.GetData(1, ref network)) return;
            
            FDM_Response response = JsonSerializer.Deserialize<FDM_Response>(_response);

            if (response.Q.Count != network.Graph.Ne)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Number of edges in network does not match number of force densities in response. It is safe to ignore this warning if you just changed the topology of your network");
                return;
            }

            for (int i = 0; i < network.Graph.Nn; i++)
            {
                nodes.Add(new Point3d(response.X[i], response.Y[i], response.Z[i]));
            }
            // get nodes



            foreach (Edge edge in network.Graph.Edges)
            {
                edges.Add(new LineCurve(nodes[network.Graph.Nodes.IndexOf(edge.Start)], nodes[network.Graph.Nodes.IndexOf(edge.End)]));
            }

            DA.SetDataList(0, nodes);
            DA.SetDataList(1, edges);
            DA.SetData(2, response.Iter);
            DA.SetDataList(3, response.Q);
            DA.SetData(4, response.Loss);
            DA.SetDataList(5, response.Losstrace);
            if (response.NodeTrace != null)
            {
                DA.SetDataTree(6, response.NodeTraceToTree());
            }


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
                        Point3d pt = new Point3d(iter[0][i], iter[1][i], iter[2][i]);
                        GH_Point gh_Pt = new GH_Point(pt);
                        Trace.Append(gh_Pt, new GH_Path(counter));
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
                return Properties.Resources.Receive;
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