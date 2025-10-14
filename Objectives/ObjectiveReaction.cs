using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using Ariadne.Objectives;
using Ariadne.FDM;
using Rhino.DocObjects;
using System.Linq;
using Ariadne.Graphs;
using System.Drawing;

namespace Ariadne.Objectives
{
    public class ObectiveReaction : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the ObjectiveTarget class.
        /// </summary>
        public ObectiveReaction()
          : base("Reaction Objective", "Reaction",
              "Minimize the angular deviation of a reaction at an anchor from a target vector (magnitude upper bound optional)",
              "Ariadne", "ObjectiveFunctions")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Nodes", "Nodes", "Anchor nodes to evaluate reaction at", GH_ParamAccess.list);
            pManager.AddVectorParameter("Target Vector", "Target", "Target vector(s) for reaction direction. If more than one vector is provided then the count must match the number of anchor points.", GH_ParamAccess.list);
            pManager.AddBooleanParameter("Use Magnitude", "UseMag", "If true then the magnitude of the reaction force is also considered in the objective function (i.e. force must be less than target magnitude)", GH_ParamAccess.item, false);
            pManager.AddNumberParameter("Weight", "Weight", "Weight of this objective function", GH_ParamAccess.item, 1.0);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Target Geometry Objective Function", "OBJ", "Target Geometry Objective Function", GH_ParamAccess.item);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Node> nodes = new();
            List<Vector3d> targets = new();
            bool useMag = false;
            double weight = 1.0;

                
            if(!DA.GetDataList(0, nodes)) return;
            if(!DA.GetDataList(1, targets)) return;
            DA.GetData(2, ref useMag);
            if (!DA.GetData(3, ref weight)) return;

            OBJReaction obj = new(weight, nodes, targets, useMag);


            if (obj.IsValid)
            {                    
                DA.SetData(0, obj);
            }
            else 
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Objective function not valid. Check inputs. This objective is only valid on anchor nodes.");
                return;
            }
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override Bitmap Icon => Properties.Resources.Reactions;

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("81BADEDE-04C7-4BEB-BB3A-4D562C71320F"); }
        }
    }
}