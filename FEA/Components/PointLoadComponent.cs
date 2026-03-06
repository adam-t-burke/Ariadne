using System;
using System.Collections.Generic;
using System.Drawing;
using Ariadne.Graphs;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    public class PointLoadComponent : GH_Component
    {
        public PointLoadComponent()
            : base("Point Load", "PtLoad",
                "Create concentrated nodal loads. Useful for bars or advanced solid use; surface loads are preferred for solids.",
                "Theseus-FEA", "Setup")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddVectorParameter("Forces", "F", "Force vectors", GH_ParamAccess.list);
            pManager.AddGenericParameter("Nodes", "N", "Target nodes or node indices", GH_ParamAccess.list);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Node Loads", "NL", "Canonical nodal load objects", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            var forces = new List<Vector3d>();
            var nodeInputs = new List<object>();

            if (!DA.GetDataList(0, forces) || forces.Count == 0) return;
            if (!DA.GetDataList(1, nodeInputs) || nodeInputs.Count == 0) return;

            var loads = new List<FeaLoad>(Math.Max(forces.Count, nodeInputs.Count));
            int count = Math.Max(forces.Count, nodeInputs.Count);

            for (int i = 0; i < count; i++)
            {
                var force = forces[Math.Min(i, forces.Count - 1)];
                int nodeIndex = ParseNodeIndex(nodeInputs[Math.Min(i, nodeInputs.Count - 1)]);
                if (nodeIndex < 0)
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                        "Nodes input must contain Node objects or integer node indices.");
                    return;
                }

                loads.Add(new FeaLoad
                {
                    NodeIndex = nodeIndex,
                    Force = force,
                });
            }

            DA.SetDataList(0, loads);
        }

        private static int ParseNodeIndex(object input)
        {
            if (input is GH_ObjectWrapper wrapper)
                input = wrapper.Value!;

            return input switch
            {
                Node node => node.Index,
                GH_Integer ghInt => ghInt.Value,
                int value => value,
                _ => -1,
            };
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-200000000011");
    }
}
