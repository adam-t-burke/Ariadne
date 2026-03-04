using System;
using System.Collections.Generic;
using System.Drawing;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Ariadne.Graphs;

namespace Ariadne.FEA.Components
{
    public class ConstructFEANetworkComponent : GH_Component
    {
        public ConstructFEANetworkComponent()
            : base("Construct FEA Network", "FEA Network",
                "Create an FEA network from edges, supports, materials, and sections",
                "Theseus-FEA", "Design")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Edges", "E", "Edge curves defining the truss", GH_ParamAccess.tree);
            pManager.AddNumberParameter("Tolerance", "Tol", "Edge merging tolerance", GH_ParamAccess.item, 0.001);
            pManager.AddGenericParameter("Supports", "Sup", "Support conditions", GH_ParamAccess.list);
            pManager.AddGenericParameter("Material", "Mat", "Material properties", GH_ParamAccess.list);
            pManager.AddGenericParameter("Section", "Sec", "Section properties", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Mat Assign", "MA", "Material index per edge (default 0)", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Sec Assign", "SA", "Section index per edge (default 0)", GH_ParamAccess.list);

            pManager[3].Optional = true;
            pManager[4].Optional = true;
            pManager[5].Optional = true;
            pManager[6].Optional = true;
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Network", "N", "FEA network", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            DA.DisableGapLogic();

            GH_Structure<GH_Curve> edges = new();
            double tol = 0.001;
            var supports = new List<FeaSupport>();
            var materials = new List<FeaMaterial>();
            var sections = new List<FeaSection>();
            var matAssign = new List<int>();
            var secAssign = new List<int>();

            if (!DA.GetDataTree(0, out edges)) return;
            DA.GetData(1, ref tol);
            if (!DA.GetDataList(2, supports)) return;

            if (!DA.GetDataList(3, materials) || materials.Count == 0)
                materials = [FeaMaterial.Steel()];
            if (!DA.GetDataList(4, sections) || sections.Count == 0)
                sections = [new FeaSection()];

            DA.GetDataList(5, matAssign);
            DA.GetDataList(6, secAssign);

            var graph = new Graph(edges, tol);

            if (supports.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "At least one support is required.");
                return;
            }

            var network = new FEA_Network(
                graph, supports, materials, sections,
                matAssign.ToArray(), secAssign.ToArray(), tol);

            if (!network.Valid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, network.ValidationMessage);
                return;
            }

            DA.SetData(0, network);
        }

        protected override Bitmap Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-000000000004");
    }
}
