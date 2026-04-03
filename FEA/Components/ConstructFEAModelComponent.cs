using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Ariadne.Graphs;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    public class ConstructFEAModelComponent : GH_Component
    {
        public ConstructFEAModelComponent()
            : base("Construct FEA Model", "FEA Model",
                "Assemble an FEA model from element objects, supports, and loads",
                "Theseus-FEA", "Setup")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Elements", "E", "Element objects: BarNetwork, SolidMesh, and/or ShellMesh", GH_ParamAccess.list);
            pManager.AddGenericParameter("Supports", "Sup", "Support conditions (from Point Support or face support components)", GH_ParamAccess.list);
            pManager.AddGenericParameter("Node Loads", "NL", "Nodal loads (from Point Load component)", GH_ParamAccess.list);

            pManager[2].Optional = true;
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("FEA Model", "M", "Unified FEA model", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            var elementInputs = new List<object>();
            var supportInputs = new List<object>();
            var loads = new List<FeaLoad>();

            if (!DA.GetDataList(0, elementInputs) || elementInputs.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "At least one element object is required.");
                return;
            }

            if (!DA.GetDataList(1, supportInputs) || supportInputs.Count == 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "At least one support is required.");
                return;
            }

            DA.GetDataList(2, loads);

            // ── Classify element objects ──
            BarNetwork? barNetwork = null;
            SolidMesh? solidMesh = null;
            ShellMesh? shellMeshData = null;
            var pointSupports = new List<FeaSupport>();
            var solidSupportRegions = new List<FeaSolidSupportRegion>();

            foreach (var raw in elementInputs)
            {
                object item = raw;
                if (item is GH_ObjectWrapper w) item = w.Value!;

                switch (item)
                {
                    case BarNetwork bn:
                        if (barNetwork != null)
                        {
                            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Only one BarNetwork is supported per model.");
                            return;
                        }
                        barNetwork = bn;
                        break;

                    case SolidMesh sm:
                        if (solidMesh != null)
                        {
                            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Only one SolidMesh is supported per model.");
                            return;
                        }
                        solidMesh = sm;
                        break;

                    case ShellMesh shm:
                        if (shellMeshData != null)
                        {
                            AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Only one ShellMesh is supported per model.");
                            return;
                        }
                        shellMeshData = shm;
                        break;

                    default:
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                            $"Unsupported element type: {item.GetType().Name}. Expected BarNetwork, SolidMesh, or ShellMesh.");
                        return;
                }
            }

            // ── Classify supports ──
            foreach (var raw in supportInputs)
            {
                object item = raw;
                if (item is GH_ObjectWrapper w) item = w.Value!;

                switch (item)
                {
                    case FeaSupport support:
                        pointSupports.Add(support);
                        break;
                    case FeaSolidSupportRegion region:
                        solidSupportRegions.Add(region);
                        break;
                    default:
                        AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                            "Supports must be FeaSupport or FeaSolidSupportRegion objects.");
                        return;
                }
            }

            if (solidMesh == null && solidSupportRegions.Count > 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error, "Solid support regions require a SolidMesh element.");
                return;
            }

            // ── Build the model ──
            var model = new FEA_Model
            {
                Supports = pointSupports,
                SolidSupportRegions = solidSupportRegions,
                Loads = loads,
            };

            if (barNetwork != null)
            {
                model.BarGraph = barNetwork.Graph;
                model.BarNetworkData = barNetwork;
                model.Materials = barNetwork.Materials;
                model.Sections = barNetwork.Sections;
                model.MaterialAssignment = barNetwork.MaterialAssignment;
                model.SectionAssignment = barNetwork.SectionAssignment;
            }

            if (solidMesh != null)
            {
                model.SolidMesh = solidMesh;
                var mat = solidMesh.Material ?? FeaMaterial.Steel();
                if (barNetwork == null)
                {
                    model.Materials = [mat];
                    model.MaterialAssignment = Enumerable.Repeat(0, solidMesh.Ne).ToArray();
                }
            }

            if (shellMeshData != null)
            {
                model.ShellMesh = shellMeshData.OriginalMesh;
                model.ShellMeshData = shellMeshData;
                model.ShellThickness = shellMeshData.NodeThicknesses;
                if (barNetwork == null && solidMesh == null)
                {
                    model.Materials = [shellMeshData.Material];
                    model.MaterialAssignment = Enumerable.Repeat(0, shellMeshData.OriginalFaceCount).ToArray();
                }
            }

            model.Initialize();

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
