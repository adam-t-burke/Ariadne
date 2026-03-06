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
    public class DeconstructFEAModelComponent : GH_Component
    {
        public DeconstructFEAModelComponent()
            : base("Deconstruct FEA Model", "Decon Model",
                "Extract nodes, elements, supports, and properties from a unified FEA model",
                "Theseus-FEA", "Utilities")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("FEA Model", "M", "Unified FEA model to deconstruct", GH_ParamAccess.item);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Nodes", "Nodes", "All node objects", GH_ParamAccess.list);
            pManager.AddGenericParameter("Bar Elements", "Bar", "Bar elements with material/section info", GH_ParamAccess.list);
            pManager.AddGenericParameter("Solid Mesh", "SM", "Solid mesh (null if no solid elements)", GH_ParamAccess.item);
            pManager.AddGenericParameter("Supports", "Sup", "Support boundary condition objects", GH_ParamAccess.list);
            pManager.AddGenericParameter("Node Loads", "NL", "Canonical nodal load objects", GH_ParamAccess.list);
            pManager.AddVectorParameter("Loads", "L", "Load vectors", GH_ParamAccess.list);
            pManager.AddGenericParameter("Load Nodes", "LN", "Node objects corresponding to each load", GH_ParamAccess.list);
            pManager.AddGenericParameter("Support Nodes", "SN", "Node objects at support locations", GH_ParamAccess.list);
            pManager.AddGenericParameter("Free Nodes", "FN", "Unsupported node objects", GH_ParamAccess.list);
            pManager.AddGenericParameter("Materials", "Mat", "Material definitions", GH_ParamAccess.list);
            pManager.AddGenericParameter("Sections", "Sec", "Section definitions", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Num Nodes", "Nn", "Total number of nodes", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Num Elements", "Ne", "Total number of elements", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            FEA_Model? model = null;
            if (!DA.GetData(0, ref model)) return;
            if (model == null)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Null FEA model.");
                return;
            }

            var allNodes = new List<Node>();
            var barElements = new List<FeaBarElement>();

            if (model.HasBarElements)
            {
                var graph = model.BarGraph!;
                allNodes.AddRange(graph.Nodes);

                for (int e = 0; e < graph.Ne; e++)
                {
                    int matIdx = model.MaterialAssignment[e];
                    int secIdx = model.SectionAssignment[e];
                    barElements.Add(new FeaBarElement(
                        graph.Edges[e], matIdx, secIdx,
                        model.Materials[matIdx], model.Sections[secIdx]));
                }
            }

            SolidMesh? solidMesh = null;
            if (model.HasSolidElements)
            {
                solidMesh = model.SolidMesh!;
                for (int i = 0; i < solidMesh.TetNodes.Count; i++)
                {
                    allNodes.Add(new Node
                    {
                        Value = solidMesh.TetNodes[i],
                        Index = model.HasBarElements
                            ? model.BarGraph!.Nn + i
                            : i
                    });
                }
            }

            var supportSet = new HashSet<int>(model.SupportNodeIndices);
            var supportNodes = new List<Node>(model.SupportNodeIndices.Count);
            var freeNodes = new List<Node>();
            var loadVectors = model.Loads.Select(load => load.Force).ToList();
            var loadNodes = model.Loads.Select(load => model.CreateNode(load.NodeIndex)).ToList();

            foreach (var node in allNodes)
            {
                if (supportSet.Contains(node.Index))
                    supportNodes.Add(node);
                else
                    freeNodes.Add(node);
            }

            var sections = model.HasSolidElements && !model.HasBarElements
                ? new List<FeaSection>()
                : model.Sections;

            DA.SetDataList(0, allNodes);
            DA.SetDataList(1, barElements);
            DA.SetData(2, solidMesh);
            DA.SetDataList(3, model.GetSupportDefinitions());
            DA.SetDataList(4, model.Loads);
            DA.SetDataList(5, loadVectors);
            DA.SetDataList(6, loadNodes);
            DA.SetDataList(7, supportNodes);
            DA.SetDataList(8, freeNodes);
            DA.SetDataList(9, model.Materials);
            DA.SetDataList(10, sections);
            DA.SetData(11, model.TotalNodes);
            DA.SetData(12, model.TotalElements);
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-200000000002");
    }
}
