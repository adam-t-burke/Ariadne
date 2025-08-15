using System;
using System.Collections.Generic;
using Grasshopper.Kernel;
using Rhino.Geometry;
using Ariadne.Optimization;
using System.Text.Json;
using Grasshopper.Kernel.Types;
using System.Text.Json.Serialization;
using System.Diagnostics.Eventing.Reader;
using Grasshopper.Kernel.Parameters;
using System.Drawing;
using System.Security.Cryptography;
using System.Text;

namespace Ariadne.FDM
{
    public class ConstructProblem : GH_Component
    {
        JsonSerializerOptions options = new() { WriteIndented = true, DefaultIgnoreCondition = JsonIgnoreCondition.WhenWritingNull};
        /// <summary>
        /// Initializes a new instance of the FDM_Problem class.
        /// </summary>
        public ConstructProblem()
          : base("FDM Problem", "FDM Problem",
              "Collect all the elements required form find using the Force Density Method",
              "Ariadne", "Design")
        {
        }


        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("FDM Network", "Network", "FDM Network", GH_ParamAccess.item);
            pManager.AddNumberParameter("Force Densities", "q", "Force density of edges", GH_ParamAccess.list, 10.0);
            pManager.AddVectorParameter("Loads", "Loads", "Loads on the free nodes of the network", GH_ParamAccess.list, new Vector3d(0.0,0.0,0.0));
            pManager.AddIntegerParameter("Nodes to apply loads to", "Load Nodes", "Indices of nodes to apply loads to", GH_ParamAccess.list);
            pManager.AddGenericParameter("Optimization Parameters", "Params", "Parameters and objectives for optimization.", GH_ParamAccess.item);
            pManager.AddGenericParameter("Variable Anchors", "VarAnchors", "Variable anchors for optimization", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Cancel Form-Finding", "Cancel", "Cancel the form-finding process", GH_ParamAccess.item, false);
    
            pManager[3].Optional = true;
            pManager[4].Optional = true;
            pManager[5].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("FDM Problem", "FDM", "All of the elements required for form-finding using the Force Density Method", GH_ParamAccess.item);
            pManager.AddTextParameter("Message", "Message", "Message to send to the server", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        private string lastInputHash = null;
        private FDM_Problem lastProblem = null;

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            //Initialize
            FDM_Network Network = new FDM_Network();
            List<double> Q = new List<double>();
            List<Vector3d> P = new List<Vector3d>();
            List<int> LoadNodes = new List<int>();
            OBJParameters Optimization = null;
            FDM_Problem problem;
            List<Anchor> anchors = new List<Anchor>();
            bool cancel = false;

            //assign
            if (!DA.GetData(0, ref Network)) return;
            if (!DA.GetDataList(1, Q)) return;
            if (!DA.GetDataList(2, P)) return;
            DA.GetDataList(3, LoadNodes);
            DA.GetData(4, ref Optimization);
            DA.GetDataList(5, anchors);
            DA.GetData(6, ref cancel);

            if (cancel)
            {
                DA.SetData(1, "cancel");
                return;
            }


            if (P.Count > 1)
            {
                if (!(P.Count == Network.Free.Count || P.Count == LoadNodes.Count))
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Number of loads must match number of load nodes, number of free nodes, or be equal to one.");
                    return;
                }
            }
            if (Optimization == null && anchors.Count > 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Variable Anchors can only be used with optimization. Only a FDM direct solve will be done.");
            }

            string inputHash = ComputeHash(Network, Q, P, LoadNodes, Optimization, anchors);

            if (inputHash == lastInputHash && lastProblem != null)
            {
                DA.SetData(0, lastProblem);
                DA.SetData(1, "cached");
                return;
            }

            if (Optimization != null && LoadNodes.Count > 0 && anchors.Count > 0)
            {
                problem = new FDM_Problem(Network, Q, P, Optimization, LoadNodes, anchors);
            }
            else if (Optimization != null && LoadNodes.Count > 0 && anchors.Count == 0)
            {
                problem = new FDM_Problem(Network, Q, P, Optimization, LoadNodes);
            }
            else if (Optimization != null && LoadNodes.Count == 0 && anchors.Count > 0)
            {
                problem = new FDM_Problem(Network, Q, P, Optimization, anchors);
            }
            else if (Optimization != null && LoadNodes.Count == 0 && anchors.Count == 0)
            {
                problem = new FDM_Problem(Network, Q, P, Optimization);
            }
            else if (Optimization == null && LoadNodes.Count > 0)
            {
                problem = new FDM_Problem(Network, Q, P, LoadNodes);
            }            
            else if (Optimization == null && LoadNodes.Count == 0)
            {
                problem = new FDM_Problem(Network, Q, P);
            }
            else
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Something went wrong with the inputs.");
                return;
            }


            string message = JsonSerializer.Serialize<FDM_Problem>(problem, options);

            lastInputHash = inputHash;
            lastProblem = problem;
            DA.SetData(0, problem);
            DA.SetData(1, message);
        }

        private string ComputeHash(
    FDM_Network network,
    List<double> Q,
    List<Vector3d> P,
    List<int> LoadNodes,
    OBJParameters Optimization,
    List<Anchor> anchors)
{
    var sb = new StringBuilder();

    // Network topology: node and edge count, plus node positions
    if (network?.Graph?.Nodes != null)
    {
        sb.Append("N:").Append(network.Graph.Nodes.Count);
        foreach (var node in network.Graph.Nodes)
        {
            var pt = node.Value;
            sb.AppendFormat("({0:F6},{1:F6},{2:F6})", pt.X, pt.Y, pt.Z);
        }
    }
    if (network?.Graph?.Edges != null)
    {
        sb.Append("E:").Append(network.Graph.Edges.Count);
        foreach (var edge in network.Graph.Edges)
        {
            sb.AppendFormat("({0},{1})", edge.Start.Index, edge.End.Index);
        }
    }

    // Force densities
    sb.Append("Q:").Append(Q.Count);
    foreach (var q in Q)
        sb.AppendFormat("{0:F6}", q);

    // Loads
    sb.Append("P:").Append(P.Count);
    foreach (var v in P)
        sb.AppendFormat("({0:F6},{1:F6},{2:F6})", v.X, v.Y, v.Z);

    // Load nodes
    sb.Append("LN:").Append(LoadNodes.Count);
    foreach (var idx in LoadNodes)
        sb.Append(idx);

    // Optimization parameters (basic fields)
    if (Optimization != null)
    {
        sb.Append("Opt:");
        if (Optimization.LowerBound != null)
            foreach (var lb in Optimization.LowerBound) sb.AppendFormat("{0:F6}", lb);
        if (Optimization.UpperBound != null)
            foreach (var ub in Optimization.UpperBound) sb.AppendFormat("{0:F6}", ub);
        sb.AppendFormat("{0:F6}{1:F6}", Optimization.AbsTol, Optimization.RelTol);
        sb.Append(Optimization.MaxIterations);
        sb.Append(Optimization.UpdateFrequency);
        sb.Append(Optimization.ShowIterations);
        sb.Append(Optimization.NodeTrace);
        if (Optimization.Objectives != null)
            sb.Append("ObjCount:").Append(Optimization.Objectives.Count);
    }

    // Anchors
    sb.Append("A:").Append(anchors.Count);
    foreach (var a in anchors)
        sb.AppendFormat("({0},{1:F6},{2:F6},{3:F6})", a.NodeIndex, a.InitialX, a.InitialY, a.InitialZ);

    // Compute SHA256 hash for stability
    using var sha = SHA256.Create();
    var hashBytes = sha.ComputeHash(Encoding.UTF8.GetBytes(sb.ToString()));
    return Convert.ToHexString(hashBytes);
}



        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                return Properties.Resources.Create;

            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("954DE917-FA4B-442F-9B83-F0C83AE7E12E"); }
        }
    }
}