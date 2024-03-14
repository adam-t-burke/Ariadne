using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Ariadne.Objectives
{
    public class ObjectiveMaxForce : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ObjectiveMaxForce class.
        /// </summary>
        public ObjectiveMaxForce()
          : base("Maximum Force", "MaxForce",
              "Penalizes force values above threshold",
              "Ariadne", "ObjectiveFunctions")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddIntegerParameter("Edge Indices", "Edge Indices", "Edge Indices", GH_ParamAccess.list);
            pManager.AddNumberParameter("Maximum Force", "Force", "Target maximum force", GH_ParamAccess.list, 10000.0);
            pManager.AddNumberParameter("Weight", "W", "Weight of objective", GH_ParamAccess.item, 1.0);

            pManager[0].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Maximum Force Objective", "OBJ", "Maximum force function", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<int> edges = new List<int>();
            List<double> force = new List<double> { 1.0 };
            double weight = 1.0;

            DA.GetDataList(0, edges);
            if (!DA.GetData(1, ref force)) { return; }
            if (!DA.GetData(2, ref weight)) { return; }


            if (edges.Count > 0 && force.Count == 1)
            {
                OBJMaxforce obj = new OBJMaxforce(weight, force[0], edges);
                if (obj.IsValid)
                {
                    DA.SetData(0, obj);
                }
            }
            else if (edges.Count > 0 && force.Count > 1)
            {
                OBJMaxforce obj = new OBJMaxforce(weight, force, edges);
                if (obj.IsValid)
                {
                    DA.SetData(0, obj);
                }
                else
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Objective creation failed. Check to ensure number of forces supplied matches the number of edge indices.");
                    return;
                }
            }
            else
            {
                OBJMaxforce obj = new OBJMaxforce(weight, force);
                if (obj.IsValid)
                {
                    DA.SetData(0, obj);
                }
                else
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Objective creation failed.");
                    return;
                }
            }


        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override Bitmap Icon => Properties.Resources.maxforce;

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("EA17A78F-BB61-4B5C-861F-0F76247DCC3B"); }
        }
    }
}