using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    public class DeconstructShellResultsComponent : GH_Component
    {
        public DeconstructShellResultsComponent()
            : base("Deconstruct Shell Results", "ShellRes",
                "Break down shell-specific fields from a unified FEA result",
                "Theseus-FEA", "Solver")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Results", "Res", "Unified FEA result (must contain shell data)", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddVectorParameter("Displacements", "U", "Nodal displacement vectors", GH_ParamAccess.list);
            pManager.AddVectorParameter("Rotations", "R", "Nodal rotation vectors", GH_ParamAccess.list);
            pManager.AddVectorParameter("Reaction Forces", "RF", "Reaction force vectors at supports", GH_ParamAccess.list);
            pManager.AddVectorParameter("Reaction Moments", "RM", "Reaction moment vectors at supports", GH_ParamAccess.list);
            pManager.AddPointParameter("Deformed Nodes", "DN", "Deformed node positions", GH_ParamAccess.list);
            pManager.AddNumberParameter("S1", "S1", "Major principal membrane stress per element", GH_ParamAccess.list);
            pManager.AddNumberParameter("S2", "S2", "Minor principal membrane stress per element", GH_ParamAccess.list);
            pManager.AddNumberParameter("Von Mises", "VM", "Von Mises stress per element (max of top/bottom)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Utilization", "Util", "Stress utilization per element", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            object? input = null;
            if (!DA.GetData(0, ref input)) return;
            if (input is GH_ObjectWrapper wrapper)
                input = wrapper.Value;
            if (input is not FeaResult result)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Input must be a unified FeaResult.");
                return;
            }
            if (!result.Model.HasShellElements)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Result does not contain shell data.");
                return;
            }

            DA.SetDataList(0, result.Displacements);
            DA.SetDataList(1, result.ShellRotations ?? Array.Empty<Vector3d>());
            DA.SetDataList(2, result.Reactions);
            DA.SetDataList(3, result.ShellReactionMoments ?? Array.Empty<Vector3d>());
            DA.SetDataList(4, result.DeformedNodes);
            DA.SetDataList(5, result.ShellS1 ?? Array.Empty<double>());
            DA.SetDataList(6, result.ShellS2 ?? Array.Empty<double>());
            DA.SetDataList(7, result.ShellVonMises ?? Array.Empty<double>());
            DA.SetDataList(8, result.Utilization);
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-300000000003");
    }
}
