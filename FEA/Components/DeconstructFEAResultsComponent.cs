using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Ariadne.FEA;
using Ariadne.Graphs;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    public class DeconstructFEAResultsComponent : GH_Component
    {
        public DeconstructFEAResultsComponent()
            : base("Deconstruct FEA Results", "FEA Info",
                "Extract detailed results from an FEA solve",
                "Theseus-FEA", "Utilities")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Results", "R", "FEA solve result", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddPointParameter("Nodes", "N", "Original node positions", GH_ParamAccess.list);
            pManager.AddCurveParameter("Edges", "E", "Original edge curves (bar only; empty for solid)", GH_ParamAccess.list);
            pManager.AddPointParameter("Deformed Nodes", "DN", "Deformed node positions", GH_ParamAccess.list);
            pManager.AddVectorParameter("Displacements", "U", "Displacement vector per node", GH_ParamAccess.list);
            pManager.AddVectorParameter("Reactions", "R", "Reaction force vectors at supports", GH_ParamAccess.list);
            pManager.AddGenericParameter("Reaction Nodes", "RN", "Node objects corresponding to each reaction", GH_ParamAccess.list);
            pManager.AddNumberParameter("Axial Forces", "F", "Axial force per element (bar only; empty for solid)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Stresses", "σ", "Axial stress (bar) or von Mises (solid) per element", GH_ParamAccess.list);
            pManager.AddNumberParameter("Strains", "ε", "Axial strain per element (bar only; empty for solid)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Utilization", "Util", "Stress utilization per element", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Num Nodes", "Nn", "Number of nodes", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Num Elements", "Ne", "Number of elements", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            object? input = null;
            if (!DA.GetData(0, ref input)) return;
            if (input == null) return;

            // Unwrap GH_ObjectWrapper when Grasshopper wraps custom types for generic parameters
            if (input is GH_ObjectWrapper wrapper)
                input = wrapper.Value;

            if (input is FeaResult result)
            {
                SolveResult(DA, result);
                return;
            }

            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Input must be a FeaResult.");
        }

        private void SolveResult(IGH_DataAccess DA, FeaResult result)
        {
            var model = result.Model;
            bool isBar = model.HasBarElements;

            if (isBar)
            {
                var graph = model.BarGraph!;
                DA.SetDataList(0, graph.Nodes.ConvertAll(n => n.Value));
                DA.SetDataList(1, graph.Edges.ConvertAll(e => (LineCurve)e.Value));
                DA.SetDataList(2, result.DeformedNodes);
                DA.SetDataList(3, result.Displacements);
                DA.SetDataList(4, result.Reactions);

                var reactionNodes = new List<Node>();
                foreach (int ni in result.ReactionNodeIndices)
                    reactionNodes.Add(model.CreateNode(ni));
                DA.SetDataList(5, reactionNodes);

                DA.SetDataList(6, result.AxialForces ?? []);
                DA.SetDataList(7, result.BarStresses ?? []);
                DA.SetDataList(8, result.BarStrains ?? []);
                DA.SetDataList(9, result.Utilization);
                DA.SetData(10, graph.Nn);
                DA.SetData(11, graph.Ne);
            }
            else
            {
                var mesh = model.SolidMesh!;
                DA.SetDataList(0, mesh.TetNodes);
                DA.SetDataList(1, new List<LineCurve>());
                DA.SetDataList(2, result.DeformedNodes);
                DA.SetDataList(3, result.Displacements);
                DA.SetDataList(4, result.Reactions);

                var reactionNodes = new List<Node>();
                foreach (int ni in result.ReactionNodeIndices)
                    reactionNodes.Add(model.CreateNode(ni));
                DA.SetDataList(5, reactionNodes);

                DA.SetDataList(6, new List<double>());
                DA.SetDataList(7, result.VonMises ?? []);
                DA.SetDataList(8, new List<double>());
                DA.SetDataList(9, result.Utilization);
                DA.SetData(10, mesh.Nn);
                DA.SetData(11, mesh.Ne);
            }
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-000000000020");
    }
}
