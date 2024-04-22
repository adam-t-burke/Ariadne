using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Linq;
using Ariadne.FDM;
using System.Drawing;

namespace Ariadne.Utilities
{
    public class NetworkInfo : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the VisualizeNetwork class.
        /// </summary>
        public NetworkInfo()
          : base("NetworkInformation", "Info",
              "Information about a FDM network",
              "Ariadne", "Utilities")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("FDM Network", "Network", "FDM network class", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Nodes", "N", "Nodes", GH_ParamAccess.list);
            pManager.AddCurveParameter("Edges", "E", "Edges", GH_ParamAccess.list);
            pManager.AddNumberParameter("Edge Lengths", "L", "Length of edges", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Start Index", "iStart", "Index of start point in N", GH_ParamAccess.list);
            pManager.AddIntegerParameter("End Index", "iStart", "Index of end point in N", GH_ParamAccess.list);
            pManager.AddPointParameter("Anchors", "A", "Anchor nodes", GH_ParamAccess.list);
            pManager.AddNumberParameter("ForceDensities", "q", "Force densities of edges", GH_ParamAccess.list);
            pManager.AddIntegerParameter("NumberEdges", "Ne", "Number of edges", GH_ParamAccess.item);
            pManager.AddIntegerParameter("NumberNodes", "Nn", "Number of nodes", GH_ParamAccess.item);
            pManager.AddIntegerParameter("FreeIndices", "iN", "Indices of free points", GH_ParamAccess.list);
            pManager.AddIntegerParameter("FixedIndices", "iF", "Indices of fixed points", GH_ParamAccess.list);
            //pManager.AddNumberParameter("MemberForces", "Force", "Internal forces assuming L = stressed length", GH_ParamAccess.list);
            //pManager.AddVectorParameter("Reaction Forces", "Reactions", "Force vectors acting at anchor points", GH_ParamAccess.list);

            pManager.HideParameter(0);
            pManager.HideParameter(1);
            pManager.HideParameter(5);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            FDM_Problem fdmNetwork = new FDM_Problem();
            if (!DA.GetData(0, ref fdmNetwork)) return;

            //List<double> forces = Solver.Forces(fdmNetwork);
            //List<Vector3d> reactions = Solver.Reactions(fdmNetwork);
            List<double> lengths = UtilityFunctions.GetLengths(fdmNetwork.Network.Graph.Edges);

            List<Point3d> startPt = fdmNetwork.Network.Graph.Edges.Select(st => st.Start.Value).ToList();
            List<Point3d> endPt = fdmNetwork.Network.Graph.Edges.Select(end => end.End.Value).ToList();

            DA.SetDataList(0, fdmNetwork.Network.Graph.Nodes);
            DA.SetDataList(1, fdmNetwork.Network.Graph.Edges);
            DA.SetDataList(2, lengths);
            DA.SetDataList(3, startPt);
            DA.SetDataList(4, endPt);
            DA.SetDataList(5, fdmNetwork.Network.Anchors);
            DA.SetDataList(6, fdmNetwork.Q);
            DA.SetData(7, fdmNetwork.Network.Graph.Ne);
            DA.SetData(8, fdmNetwork.Network.Graph.Nn);
            DA.SetDataList(9, fdmNetwork.Network.Free);
            DA.SetDataList(10, fdmNetwork.Network.Fixed);
            //DA.SetDataList(11, forces);
            //DA.SetDataList(12, reactions);

        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override Bitmap Icon => Properties.Resources.Information;

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("4E02AFBA-6E24-4855-9A0F-91E2042ED22B"); }
        }
    }
}