using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;
using Ariadne.Graphs;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    public class ConstructFEAModelComponent : GH_Component
    {
        public ConstructFEAModelComponent()
            : base("Construct FEA Model", "FEA Model",
                "Create a unified FEA model from bar edges and/or solid mesh",
                "Theseus-FEA", "Setup")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Edges", "E", "Edge curves defining bar elements", GH_ParamAccess.tree);
            pManager.AddGenericParameter("Solid Mesh", "SM", "Tetrahedral solid mesh", GH_ParamAccess.item);
            pManager.AddNumberParameter("Tolerance", "Tol", "Node merging tolerance", GH_ParamAccess.item, 0.001);
            pManager.AddGenericParameter("Supports", "Sup", "Support conditions", GH_ParamAccess.list);
            pManager.AddGenericParameter("Node Loads", "NL", "Nodal loads applied to the model", GH_ParamAccess.list);
            pManager.AddGenericParameter("Material", "Mat", "Material properties", GH_ParamAccess.list);
            pManager.AddGenericParameter("Section", "Sec", "Section properties (bar elements only)", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Mat Assign", "MA", "Material index per element (default 0)", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Sec Assign", "SA", "Section index per bar element (default 0)", GH_ParamAccess.list);

            pManager[0].Optional = true;
            pManager[1].Optional = true;
            pManager[4].Optional = true;
            pManager[5].Optional = true;
            pManager[6].Optional = true;
            pManager[7].Optional = true;
            pManager[8].Optional = true;
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("FEA Model", "M", "Unified FEA model", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            DA.DisableGapLogic();

            var edgeTree = new GH_Structure<GH_Curve>();
            SolidMesh? solidMesh = null;
            double tol = 0.001;
            var supportInputs = new List<object>();
            var pointSupports = new List<FeaSupport>();
            var solidSupportRegions = new List<FeaSolidSupportRegion>();
            var loads = new List<FeaLoad>();
            var materials = new List<FeaMaterial>();
            var sections = new List<FeaSection>();
            var matAssign = new List<int>();
            var secAssign = new List<int>();

            bool hasEdges = DA.GetDataTree(0, out edgeTree) && edgeTree.DataCount > 0;
            DA.GetData(1, ref solidMesh);
            DA.GetData(2, ref tol);

            bool hasSolid = solidMesh != null && solidMesh.Ne > 0;

            if (!hasEdges && !hasSolid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                    "Provide at least one of Edges (bar) or Solid Mesh.");
                return;
            }

            if (!DA.GetDataList(3, supportInputs) || supportInputs.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                    "At least one support is required.");
                return;
            }

            foreach (var input in supportInputs)
            {
                object value = input;
                if (value is GH_ObjectWrapper wrapper)
                    value = wrapper.Value!;

                switch (value)
                {
                    case FeaSupport support:
                        pointSupports.Add(support);
                        break;
                    case FeaSolidSupportRegion region:
                        solidSupportRegions.Add(region);
                        break;
                    default:
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                            "Supports input must contain FeaSupport or FeaSolidSupportRegion objects.");
                        return;
                }
            }

            if (!hasSolid && solidSupportRegions.Count > 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                    "Solid support regions require a Solid Mesh input.");
                return;
            }

            if (hasSolid && pointSupports.Count > 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                    "Point supports on solids are concentrated nodal restraints. Prefer face support components for solids.");
            }

            DA.GetDataList(4, loads);

            if (!DA.GetDataList(5, materials) || materials.Count == 0)
                materials = [FeaMaterial.Steel()];

            if (!DA.GetDataList(6, sections) || sections.Count == 0)
                sections = [new FeaSection()];

            DA.GetDataList(7, matAssign);
            DA.GetDataList(8, secAssign);

            FEA_Model model;

            if (hasEdges && !hasSolid)
            {
                var graph = new Graph(edgeTree, tol);
                model = new FEA_Model(graph, pointSupports, materials, sections,
                    matAssign.ToArray(), secAssign.ToArray(), tol, loads, solidSupportRegions);
            }
            else if (hasSolid && !hasEdges)
            {
                model = new FEA_Model(solidMesh!, pointSupports, materials,
                    matAssign.ToArray(), tol, loads, solidSupportRegions);
            }
            else
            {
                var graph = new Graph(edgeTree, tol);
                model = new FEA_Model(graph, pointSupports, materials, sections,
                    matAssign.ToArray(), secAssign.ToArray(), tol, loads, solidSupportRegions)
                {
                    SolidMesh = solidMesh
                };
            }

            if (!model.Valid)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, model.ValidationMessage);
                return;
            }

            DA.SetData(0, model);
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-200000000001");
    }
}
