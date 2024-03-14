using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Ariadne.Objectives;
using Ariadne.FDM;
using System.Drawing;

namespace Ariadne.Objectives
{
    public class ObjectiveForcevariation : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ObjectiveForcevariation class.
        /// </summary>
        public ObjectiveForcevariation()
          : base("Force Variation", "Force",
              "Minimize variance of forces in the network",
              "Ariadne", "ObjectiveFunctions")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddIntegerParameter("Edge Indices", "Edges", "Indices of edges to equalized forces in.", GH_ParamAccess.list);
            pManager.AddNumberParameter("Weight", "W", "Weight of objective", GH_ParamAccess.item, 1.0);

            pManager[0].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Force Variation Objective", "OBJ", "Force Variation Objective Function", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<int> indices = new List<int>();
            double weight = 1.0;

            DA.GetDataList(0, indices);
            if (!DA.GetData(1, ref weight)) return;

            if (indices.Count > 0)
            {
                OBJforcevariation obj = new OBJforcevariation(weight, indices);
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
            else
            {
                OBJforcevariation obj = new OBJforcevariation(weight);
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
        protected override Bitmap Icon
        {
            get
            {
                return Properties.Resources.Forcevar;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("DDF4BBFC-685F-4D1F-A4BB-0E91B539D149"); }
        }
    }
}