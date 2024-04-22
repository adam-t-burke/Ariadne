using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Ariadne.Objectives;
using System.Drawing;
using Ariadne.Graphs;

namespace Ariadne.Objectives
{
    public class ObjectiveMaxLength : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ObjectiveMaxLength class.
        /// </summary>
        public ObjectiveMaxLength()
          : base("Maximum Length", "MaxLength",
              "Penalizes edge lengths above threshold",
              "Ariadne", "ObjectiveFunctions")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Edges", "Edges", "Edges with maximum length objective", GH_ParamAccess.list);
            pManager.AddNumberParameter("Maximum Length", "Length", "Target maximum length", GH_ParamAccess.list, 1000.0);
            pManager.AddNumberParameter("Weight", "W", "Weight of objective", GH_ParamAccess.item, 1.0);

            pManager[0].Optional = true;
            
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Max Length Objective", "OBJ", "Maximum length function", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Edge> edges = new();
            List<double> length = new List<double> { 1.0 };
            double weight = 1.0;

            DA.GetDataList(0, edges);
            if (!DA.GetDataList(1, length)) { return; }
            if (!DA.GetData(2, ref weight)) { return; }


            if (edges.Count > 0 && length.Count == 1)
            {
                OBJMaxlength obj = new OBJMaxlength(weight, length, edges);
                if (obj.IsValid)
                {
                    DA.SetData(0, obj);
                }
            }
            else if (edges.Count > 0 && length.Count > 1)
            {
                OBJMaxlength obj = new OBJMaxlength(weight, length, edges);
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
                OBJMaxlength obj = new OBJMaxlength(weight, length);
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
        protected override Bitmap Icon => Properties.Resources.maxlength;
        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("A4EF19AE-A001-4902-B8B7-65117EBC4F45"); }
        }
    }
}