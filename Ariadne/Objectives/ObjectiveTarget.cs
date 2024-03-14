using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Ariadne.Objectives;
using Ariadne.FDM;
using Rhino.DocObjects;
using System.Linq;
using Ariadne.Graphs;
using System.Drawing;

namespace Ariadne.Objectives
{
    public class ObjectiveTarget : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ObjectiveTarget class.
        /// </summary>
        public ObjectiveTarget()
          : base("Target Geometry", "TargetGeo",
              "Minimize the deviation from starting nodal positions",
              "Ariadne", "ObjectiveFunctions")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Network", "Network","Network to select target nodes from", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Nodes", "Nodes", "Nodes to target", GH_ParamAccess.list);
            pManager.AddNumberParameter("Weight", "Weight", "Weight of this objective function", GH_ParamAccess.item, 1.0);

            pManager[1].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Target Geometry Objective Function", "OBJ", "Target Geometry Objective Function", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            FDM_Network network = new FDM_Network();
            List<Node> nodes = new List<Node>();
            double weight = 1.0;

            if (!DA.GetData(0, ref network))
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Network input required for specifying target nodes.");
                return;
            }
                
            DA.GetDataList(1, nodes);
            if (!DA.GetData(2, ref weight)) return;


            if (network.Valid && nodes.Count > 0)
            {
                List<Point3d> points = nodes.Select(node => node.Position).ToList();

                List<int> indices = new List<int>(); 
                foreach(Node node in nodes)
                {
                    int index = network.Graph.Nodes.IndexOf(node);
                    indices.Add(index);
                }
                    
                OBJTarget obj = new OBJTarget(weight, indices, points);
                DA.SetData(0, obj);
            }
            else if (network.Valid)
            {
                OBJTarget obj = new OBJTarget(weight, network.Free.Select(node => node.Position).ToList());
                DA.SetData(0, obj);
            }
            else
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Network invalid.");
            }
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override Bitmap Icon => Properties.Resources.Target;

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("6AA55256-288F-49A6-9648-5314691BB659"); }
        }
    }
}