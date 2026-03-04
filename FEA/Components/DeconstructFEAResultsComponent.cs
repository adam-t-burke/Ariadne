using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Grasshopper.Kernel;
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
            pManager.AddCurveParameter("Edges", "E", "Original edge curves", GH_ParamAccess.list);
            pManager.AddPointParameter("Deformed Nodes", "DN", "Deformed node positions", GH_ParamAccess.list);
            pManager.AddCurveParameter("Deformed Edges", "DE", "Deformed edge curves", GH_ParamAccess.list);
            pManager.AddNumberParameter("Displacements", "U", "Nodal displacements (flat)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Reactions", "R", "Reaction forces (flat)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Axial Forces", "F", "Axial force per element", GH_ParamAccess.list);
            pManager.AddNumberParameter("Stresses", "σ", "Axial stress per element", GH_ParamAccess.list);
            pManager.AddNumberParameter("Strains", "ε", "Axial strain per element", GH_ParamAccess.list);
            pManager.AddNumberParameter("Utilization", "Util", "Stress utilization per element", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Num Nodes", "Nn", "Number of nodes", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Num Elements", "Ne", "Number of elements", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            FeaSolveResult result = null;
            if (!DA.GetData(0, ref result)) return;
            if (result == null) return;

            var graph = result.Network.Graph;
            DA.SetDataList(0, graph.Nodes.ConvertAll(n => n.Value));
            DA.SetDataList(1, graph.Edges.ConvertAll(e => (LineCurve)e.Value));
            DA.SetDataList(2, result.DeformedNodes);
            DA.SetDataList(3, result.DeformedEdges);
            DA.SetDataList(4, result.Displacements);
            DA.SetDataList(5, result.Reactions);
            DA.SetDataList(6, result.AxialForces);
            DA.SetDataList(7, result.Stresses);
            DA.SetDataList(8, result.Strains);
            DA.SetDataList(9, result.Utilization);
            DA.SetData(10, graph.Nn);
            DA.SetData(11, graph.Ne);
        }

        protected override Bitmap Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-000000000020");
    }
}
