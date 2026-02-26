using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Text.Json;
using System.Text.Json.Serialization;
using Ariadne.Properties;
using System.Drawing;

namespace Ariadne.Utilities
{
    /// <summary>
    /// Serializes an object to JSON text for inspection or debugging.
    /// </summary>
    public class SeralizeObjects : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the SeralizeObjects class.
        /// </summary>
        public SeralizeObjects()
          : base("SeralizeObjects", "Serialize",
              "Serialize any object to inspect its contents.",
              "Ariadne", "Utilities")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Object", "Object", "Object to seralize", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("Seralized Object", "Obj", "Seralized object", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            object obj = null;

            if (!DA.GetData(0, ref obj)) { return; }

            var options = new JsonSerializerOptions { WriteIndented = true };
            string json = JsonSerializer.Serialize(obj, options);

            DA.SetData(0, json);
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
                return Resources.Serialize;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("958DC5D7-C4B4-473D-BF4E-BC525AF54D96"); }
        }
    }
}