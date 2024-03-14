using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Ariadne.Optimization
{
    public class ConstructAnchor : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ConstructAnchor class.
        /// </summary>
        public ConstructAnchor()
          : base("Construct Anchor", "Anchor",
              "Construct anchor object for variable anchors in optimization",
              "Ariadne", "Optimization")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddIntegerParameter("Node Index", "Node", "Index of node to anchor", GH_ParamAccess.item);
            pManager.AddNumberParameter("X", "X", "X coordinate of anchor", GH_ParamAccess.item);
            pManager.AddNumberParameter("Y", "Y", "Y coordinate of anchor", GH_ParamAccess.item);
            pManager.AddNumberParameter("Z", "Z", "Z coordinate of anchor", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Variable Anchor", "VarAnchor", "Variable Anchor Object", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            int nodeIndex = 0;
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;

            if(!DA.GetData(0, ref nodeIndex)) { return; }
            if(!DA.GetData(1, ref x)) { return; }
            if(!DA.GetData(2, ref y)) { return; }
            if(!DA.GetData(3, ref z)) { return; }

            Anchor anchor = new Anchor(nodeIndex, x, y, z);

            DA.SetData(0, anchor);
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("5A04038F-2EC4-488C-BB18-77987459D39A"); }
        }
    }
}