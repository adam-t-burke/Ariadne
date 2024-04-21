using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Rhino;
using Ariadne.Optimization;
using Ariadne.FDM;
using Ariadne.Graphs;
using System.Drawing;
using System.Linq;
using System.Threading.Tasks;
using Grasshopper.Plugin;
using System.Threading;

namespace Ariadne.Utilities
{

    public class UpdateRefs : GH_Component
    {
        FDM_Network oldNet;
        FDM_Network newNet;
        RhinoDoc doc;
        GH_RhinoScriptInterface rs = new();
        GH_Document gh_doc;


        /// <summary>
        /// Initializes a new instance of the UpdateCurves class.
        /// </summary>
        public UpdateRefs()
          : base("Update Reference Geometry", "Update",
              "Update ",
              "Ariadne", "Utilities")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddBooleanParameter("Update Curves", "Update", "Button to trigger the update the drawn network to its new form found position.", GH_ParamAccess.item, false);
            pManager.AddGenericParameter("Old Network", "Old Net", "Old form found network.", GH_ParamAccess.item);
            pManager.AddGenericParameter("New Network", "New Net", "New form found network.", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Difference", "Diff", "Difference between old and new network", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {

            bool update = false;
            doc = RhinoDoc.ActiveDoc;
            gh_doc = OnPingDocument();            
            
            double diff = 0;           

            DA.GetData(0, ref update);
            if (!DA.GetData(1, ref oldNet)) return;
            if (!DA.GetData(2, ref newNet)) return;

            if (update)
            {
                try
                {                    
                    if (oldNet.Graph.Nn == newNet.Graph.Nn && oldNet.Graph.Ne == newNet.Graph.Ne)
                    {
                        diff = NetworkDifference(oldNet, newNet);
                        if (diff > 1e-3)
                        {
                            var t = Task.Run(() =>
                            {
                                Thread.Sleep(10);
                                Update();
                                Thread.Sleep(10);
                                RhinoApp.InvokeOnUiThread((Action)delegate { rs.EnableSolver(); });
                            });

                            rs.DisableSolver();
                                                          
                        }
                    }
                    else
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Remark, "Network has not changed.");
                    }
                    
                }
                catch (Exception e)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, e.Message);
                    return; 
                }
            }
            DA.SetData(0, diff);
        }

        private static double NetworkDifference(FDM_Network oldNet, FDM_Network newNet)
        {
            double difference = 0;

            for(int i = 0; i < oldNet.Graph.Nn; i++)
            {
                difference += oldNet.Graph.Nodes[i].Value.DistanceTo(newNet.Graph.Nodes[i].Value);
            }

            return difference / oldNet.Graph.Nn;
        }   


        private void Update()
        {
            foreach(Edge edge in newNet.Graph.Edges)
            {
                Curve curve = edge.Value;
                Guid guid = edge.ReferenceID;
                doc.Objects.Replace(guid, curve);
            };            
        }


        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override Bitmap Icon => Properties.Resources.update;

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("11D6D7E9-EE9F-4B5E-A35F-087F852403CC"); }
        }
    }
}