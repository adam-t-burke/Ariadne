using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Ariadne.Objectives;
using System.Drawing;
using Ariadne.Graphs;

namespace Ariadne.Objectives
{
    public class ObjectiveMinLength : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ObjectiveMinLength class.
        /// </summary>
        public ObjectiveMinLength()
          : base("Minimum Length", "MinLength",
              "Penalizes edge lengths below threshold",
              "Ariadne", "ObjectiveFunctions")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Edges", "Edges", "Edges with minimum length objective", GH_ParamAccess.list);
            pManager.AddNumberParameter("Minimum Length", "Length", "Target minimum length", GH_ParamAccess.list, 1.0);
            pManager.AddNumberParameter("Weight", "W", "Weight of objective", GH_ParamAccess.item, 1.0);

        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Minimum Length Objective", "OBJ", "Minimum length function", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Edge> edges = new();
            List<double> length = new();
            double weight = 1.0;

            DA.GetDataList(0, edges);
            if (!DA.GetDataList(1, length)) { return; }
            if (!DA.GetData(2, ref weight)) { return; }


            if (edges.Count > 1)
            {
                OBJMinlength obj = new OBJMinlength(weight, length, edges);
                if (obj.IsValid)
                {
                    DA.SetData(0, obj);
                }
            }
            else
            {
                OBJMinlength obj = new OBJMinlength(weight, length);
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
        protected override Bitmap Icon => Properties.Resources.minlength;
        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("C3805761-EB4C-4502-A4C3-495850141928"); }
        }
    }
}