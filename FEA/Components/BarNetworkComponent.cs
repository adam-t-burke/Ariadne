using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Ariadne.Graphs;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Special;
using Grasshopper.Kernel.Types;

namespace Ariadne.FEA.Components
{
    public class BarNetworkComponent : GH_Component
    {
        public BarNetworkComponent()
            : base("Linear Element Network", "LinNet",
                "Build a linear element network from curves with material and section assignments. Set Bending=true for frame (beam) elements.",
                "Theseus-FEA", "Meshing")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Edges", "E", "Edge curves defining elements", GH_ParamAccess.tree);                     // 0
            pManager.AddNumberParameter("Tolerance", "Tol", "Node merging tolerance", GH_ParamAccess.item, 0.001);               // 1
            pManager.AddBooleanParameter("Bending", "B", "Enable bending (true=frame, false=truss)", GH_ParamAccess.item, false); // 2
            pManager.AddIntegerParameter("Formulation", "BF", "Beam formulation: 0=Euler-Bernoulli, 1=Timoshenko", GH_ParamAccess.item, 0); // 3
            pManager.AddGenericParameter("Material", "Mat", "Material(s) — single or list", GH_ParamAccess.list);                 // 4
            pManager.AddGenericParameter("Section", "Sec", "Section(s) — single or list", GH_ParamAccess.list);                   // 5
            pManager.AddIntegerParameter("Mat Assign", "MA", "Material index per edge (default 0)", GH_ParamAccess.list);         // 6
            pManager.AddIntegerParameter("Sec Assign", "SA", "Section index per edge (default 0)", GH_ParamAccess.list);          // 7

            pManager[3].Optional = true;
            pManager[4].Optional = true;
            pManager[5].Optional = true;
            pManager[6].Optional = true;
            pManager[7].Optional = true;
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Network", "Net", "Linear element network for FEA model construction", GH_ParamAccess.item);
            pManager.AddGenericParameter("Nodes", "N", "Node objects with Index (wire into Point Load or Support)", GH_ParamAccess.list);
            pManager.AddCurveParameter("Edges", "E", "Edge curves (tree, preserving input structure)", GH_ParamAccess.tree);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            if (!DA.GetDataTree(0, out GH_Structure<GH_Curve> edgeTree)) return;
            if (edgeTree.DataCount == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "No edge curves provided.");
                return;
            }

            double tol = 0.001;
            bool bending = false;
            int formulationInt = 0;

            DA.GetData(1, ref tol);
            DA.GetData(2, ref bending);
            DA.GetData(3, ref formulationInt);

            var formulation = formulationInt switch
            {
                1 => BeamFormulation.Timoshenko,
                _ => BeamFormulation.EulerBernoulli,
            };

            var materials = new List<FeaMaterial>();
            var sections = new List<FeaSection>();
            var matAssign = new List<int>();
            var secAssign = new List<int>();

            if (!DA.GetDataList(4, materials) || materials.Count == 0)
                materials = [FeaMaterial.Steel()];
            if (!DA.GetDataList(5, sections) || sections.Count == 0)
                sections = [new FeaSection()];
            DA.GetDataList(6, matAssign);
            DA.GetDataList(7, secAssign);

            var graph = new Graph(edgeTree, tol);

            if (graph.Ne == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Graph construction produced no edges. Check input curves and tolerance.");
                return;
            }

            if (bending)
            {
                foreach (var sec in sections)
                {
                    if (sec.Iy <= 0 || sec.Iz <= 0 || sec.J <= 0)
                    {
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Warning,
                            "Bending is enabled but section has Iy/Iz/J = 0. Beam elements require nonzero second moments and torsion constant.");
                        break;
                    }
                }
            }

            var network = new BarNetwork(graph, bending, formulation,
                materials, sections, matAssign.ToArray(), secAssign.ToArray());

            string type = bending ? $"Frame ({formulation})" : "Truss";
            AddRuntimeMessage(GH_RuntimeMessageLevel.Remark,
                $"{type} network: {graph.Nn} nodes, {graph.Ne} edges.");

            DA.SetData(0, network);
            DA.SetDataList(1, graph.Nodes);
            DA.SetDataTree(2, graph.OutputEdgeTree);
        }

        public override void AddedToDocument(GH_Document document)
        {
            base.AddedToDocument(document);
            AutoAttachValueList(document, 3, "Beam Formulation",
                ("Euler-Bernoulli", 0), ("Timoshenko", 1));
        }

        private void AutoAttachValueList(GH_Document doc, int paramIndex, string name,
            params (string label, int value)[] entries)
        {
            var param = Params.Input[paramIndex];
            if (param.SourceCount > 0) return;

            var valueList = new GH_ValueList();
            valueList.CreateAttributes();
            valueList.NickName = name;
            valueList.ListItems.Clear();
            foreach (var (label, value) in entries)
                valueList.ListItems.Add(new GH_ValueListItem(label, value.ToString()));

            valueList.Attributes.Pivot = new System.Drawing.PointF(
                param.Attributes.InputGrip.X - valueList.Attributes.Bounds.Width - 30,
                param.Attributes.InputGrip.Y - valueList.Attributes.Bounds.Height / 2);

            doc.AddObject(valueList, false);
            param.AddSource(valueList);
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-100000000052");
    }
}
