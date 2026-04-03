using System;
using System.Collections.Generic;
using System.Drawing;
using Ariadne.Graphs;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    public class SupportComponent : GH_Component
    {
        public SupportComponent()
            : base("Point Support", "PtSup",
                "Define a concentrated nodal support from a Node object or Point3d location.",
                "Theseus-FEA", "Setup")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Node", "N", "Node object (from mesh/network deconstruct) or Point3d", GH_ParamAccess.item);
            pManager.AddBooleanParameter("Fix X", "X", "Constrain X translation", GH_ParamAccess.item, true);
            pManager.AddBooleanParameter("Fix Y", "Y", "Constrain Y translation", GH_ParamAccess.item, true);
            pManager.AddBooleanParameter("Fix Z", "Z", "Constrain Z translation", GH_ParamAccess.item, true);
            pManager.AddBooleanParameter("Fix RX", "RX", "Constrain X rotation (shell/beam only)", GH_ParamAccess.item, false);
            pManager.AddBooleanParameter("Fix RY", "RY", "Constrain Y rotation (shell/beam only)", GH_ParamAccess.item, false);
            pManager.AddBooleanParameter("Fix RZ", "RZ", "Constrain Z rotation (shell/beam only)", GH_ParamAccess.item, false);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Support", "S", "FEA support", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            object? nodeInput = null;
            bool fixX = true, fixY = true, fixZ = true;
            bool fixRX = false, fixRY = false, fixRZ = false;

            if (!DA.GetData(0, ref nodeInput)) return;
            DA.GetData(1, ref fixX);
            DA.GetData(2, ref fixY);
            DA.GetData(3, ref fixZ);
            DA.GetData(4, ref fixRX);
            DA.GetData(5, ref fixRY);
            DA.GetData(6, ref fixRZ);

            if (!fixX && !fixY && !fixZ && !fixRX && !fixRY && !fixRZ)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Support has no constrained DOFs.");
                return;
            }

            if (nodeInput is GH_ObjectWrapper wrapper)
                nodeInput = wrapper.Value;

            Point3d location;
            int? nodeIndex = null;

            switch (nodeInput)
            {
                case Node node:
                    location = node.Value;
                    nodeIndex = node.Index;
                    break;
                case GH_Point ghPt:
                    location = ghPt.Value;
                    break;
                case Point3d pt:
                    location = pt;
                    break;
                default:
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                        "Input must be a Node object (from mesh/network deconstruct) or a Point3d.");
                    return;
            }

            DA.SetData(0, new FeaSupport
            {
                Location = location,
                NodeIndex = nodeIndex,
                FixX = fixX,
                FixY = fixY,
                FixZ = fixZ,
                FixRX = fixRX,
                FixRY = fixRY,
                FixRZ = fixRZ,
            });
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-000000000001");
    }
}
