using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Ariadne.Utilities;
using Grasshopper.Kernel.Parameters;
using Grasshopper.Kernel.Types;
using Ariadne.Objectives;
using System.Drawing;

namespace Ariadne.Optimization
{
    /// <summary>
    /// Parameters for remote optimization
    /// </summary>
    public class ConstructOptimizationParameters : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the OptimizationParameters class.
        /// </summary>
        public ConstructOptimizationParameters()
          : base("FDM Optimization", "Optimization",
              "Optimization Parameters",
              "Ariadne", "Optimization")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Objective Functions", "Obj", "Objective function(s) to minimize.", GH_ParamAccess.list);
            pManager.AddNumberParameter("Lower Bound", "LB", "Lower bound of force densities", GH_ParamAccess.item, 1.0);
            pManager.AddNumberParameter("Upper Bound", "UB", "Upper bound of force densities", GH_ParamAccess.item, 100.0);
            pManager.AddNumberParameter("Absolute Tolerance", "AbsTol", "Absolute stopping tolerance", GH_ParamAccess.item, 1e-6);
            pManager.AddNumberParameter("Relative Tolerance", "RelTol", "Relative stopping tolerance", GH_ParamAccess.item, 1e-6);
            pManager.AddIntegerParameter("Maximum Iterations", "MaxIter", "Maximum number of iterations", GH_ParamAccess.item, 400);
            pManager.AddIntegerParameter("Update Frequency", "Frequency", "Frequency of return reports", GH_ParamAccess.item, 20);
            pManager.AddBooleanParameter("Show Iterations", "ShowIter", "Show intermittent solutions", GH_ParamAccess.item, false);
            pManager.AddBooleanParameter("Trace Optimization Geometry", "Node Trace", "Trace node positions through optimization", GH_ParamAccess.item, false);

            pManager[0].Optional = true;
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Optimization Parameters", "Params", "Parameters for optimization", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            double lb = 1e-2;
            double ub = 1e-0;
            double abstol = 1e-4;
            double reltol = 1e-4;
            int maxiter = 400;
            int freq = 20;
            bool show = true;
            bool trace = false;

            List<OBJ> objs = new List<OBJ>();

            DA.GetDataList(0, objs);
            DA.GetData(1, ref lb);
            DA.GetData(2, ref ub);
            DA.GetData(3, ref abstol);
            DA.GetData(4, ref reltol);
            DA.GetData(5, ref maxiter);
            DA.GetData(6, ref freq);
            DA.GetData(7, ref show);
            DA.GetData(8, ref trace);


            if (show && freq < 10)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Update frequency is small--may cause visualization issues for large networks");
            }
            if (freq < 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Update frequency must be positive");
                return;
            }

            OBJParameters objparams = new OBJParameters(lb, ub, abstol, reltol, objs, show, freq, maxiter, trace);
            DA.SetData(0, objparams);
         


            
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override Bitmap Icon => Properties.Resources.parameters;

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("C94AAD07-60DA-4A0A-8E30-1CE1B46B685A"); }
        }
    }
}