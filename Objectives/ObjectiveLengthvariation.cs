using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Ariadne.Objectives;
using Ariadne.FDM;
using System.Drawing;
using Ariadne.Graphs;

namespace Ariadne.Objectives
{
    public class ObjectiveLengthvariation : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ObjectiveLengthvariation class.
        /// </summary>
        public ObjectiveLengthvariation()
          : base("Length Variation", "LengthVar",
              "Minimize variance of lengths in the network",
              "Ariadne", "ObjectiveFunctions")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Edges", "Edges", "Edges to equalize lengths.", GH_ParamAccess.list);
            pManager.AddNumberParameter("Weight", "W", "Weight of objective", GH_ParamAccess.item, 1.0);

            pManager[0].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Length Variation Objective", "OBJ", "Length Objective Function", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Edge> indices = new();
            double weight = 1.0;

            DA.GetDataList(0, indices);
            if (!DA.GetData(1, ref weight)) return;

            if (indices.Count > 0)
            {
                OBJlengthvariation obj = new OBJlengthvariation(weight, indices);
                if (obj.IsValid)
                {
                    DA.SetData(0, obj);
                }
                else
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid Objective");
                    return;
                }
            }
            else
            {
                OBJlengthvariation obj = new OBJlengthvariation(weight);
                if (obj.IsValid)
                {
                    DA.SetData(0, obj);
                }
                else
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Invalid Objective");
                    return;
                }
            }

        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override Bitmap Icon => Properties.Resources.lengthvar;

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("2F4C6C0C-8D39-4764-979C-FEF17F1DB38D"); }
        }
    }
}