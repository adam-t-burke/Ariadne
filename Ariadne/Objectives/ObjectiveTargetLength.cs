using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Ariadne.FDM;
using System.Drawing;

namespace Ariadne.Objectives
{
    public class ObjectiveTargetLength : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ObjectiveTargetLength class.
        /// </summary>
        public ObjectiveTargetLength()
          : base("Target Length", "Target Length",
              "Match a specified target length for each edge.",
              "Ariadne", "ObjectiveFunctions")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddIntegerParameter("Edge Indices", "Edge Indices", "Indices of edges to equalized forces in.", GH_ParamAccess.list);
            pManager.AddNumberParameter("Target Length", "Length(s)", "Target length for each edge", GH_ParamAccess.list, 1.0);
            pManager.AddNumberParameter("Weight", "Weight", "Weight of objective", GH_ParamAccess.item, 1.0);
            

            pManager[0].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Target Length Objective", "OBJ", "Target length function", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double weight = 1.0;
            List<double> targetLengths = new List<double>();
            List<int> edgeIndices = new List<int>();

            DA.GetDataList(0, edgeIndices);
            DA.GetDataList(1, targetLengths);
            DA.GetData(2, ref weight);
            
            if (edgeIndices.Count == 0 && targetLengths.Count > 1)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Number of target lengths must be 1 when no edge indices are provided.");
                return;
            }
            else if(edgeIndices.Count > 0)
            {
                OBJTargetLength obj = new OBJTargetLength(weight, targetLengths, edgeIndices);
                if (obj.IsValid)
                {
                    DA.SetData(0, obj);
                }
                else
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Objective creation failed. Check to make sure the number of edge indices and target lengths supplied are equal.");
                    return;
                }
            }
            else
            {
                OBJTargetLength obj = new OBJTargetLength(weight, targetLengths);
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
                //You can add image files to your project resources and access them like this:
                return Properties.Resources.Target_Length;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("BB794D01-517E-44F8-8D3C-81801EF74892"); }
        }
    }
}