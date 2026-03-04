using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

namespace Ariadne.FEA.Components
{
    public class FeaOptConfigComponent : GH_Component
    {
        public FeaOptConfigComponent()
            : base("FEA Optimization Config", "FEA Opt",
                "Configure FEA optimization parameters and objectives",
                "Theseus-FEA", "Design")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Objectives", "Obj", "FEA objectives", GH_ParamAccess.tree);
            pManager.AddBooleanParameter("Opt Positions", "Pos", "Optimize node positions", GH_ParamAccess.item, true);
            pManager.AddBooleanParameter("Opt Areas", "Area", "Optimize cross-section areas", GH_ParamAccess.item, false);
            pManager.AddBooleanParameter("Opt Supports", "Sup", "Optimize support positions", GH_ParamAccess.item, false);
            pManager.AddIntegerParameter("Max Iterations", "Iter", "Maximum iterations", GH_ParamAccess.item, 500);
            pManager.AddNumberParameter("Abs Tol", "ATol", "Absolute tolerance", GH_ParamAccess.item, 1e-6);
            pManager.AddNumberParameter("Rel Tol", "RTol", "Relative tolerance", GH_ParamAccess.item, 1e-6);
            pManager.AddNumberParameter("Barrier Weight", "BW", "Barrier weight", GH_ParamAccess.item, 10.0);
            pManager.AddNumberParameter("Barrier Sharpness", "BS", "Barrier sharpness", GH_ParamAccess.item, 10.0);
            pManager.AddIntegerParameter("Report Frequency", "RF", "Progress report frequency", GH_ParamAccess.item, 10);
            pManager.AddBooleanParameter("Run", "Run", "Run optimization", GH_ParamAccess.item, false);
            pManager.AddBooleanParameter("Stream Preview", "SP", "Stream intermediate results", GH_ParamAccess.item, true);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddGenericParameter("Config", "C", "FEA optimization config", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            GH_Structure<IGH_Goo> objTree = new();
            bool optPos = true, optAreas = false, optSup = false;
            int maxIter = 500, reportFreq = 10;
            double absTol = 1e-6, relTol = 1e-6, bw = 10.0, bs = 10.0;
            bool run = false, stream = true;

            if (!DA.GetDataTree(0, out objTree)) return;
            DA.GetData(1, ref optPos);
            DA.GetData(2, ref optAreas);
            DA.GetData(3, ref optSup);
            DA.GetData(4, ref maxIter);
            DA.GetData(5, ref absTol);
            DA.GetData(6, ref relTol);
            DA.GetData(7, ref bw);
            DA.GetData(8, ref bs);
            DA.GetData(9, ref reportFreq);
            DA.GetData(10, ref run);
            DA.GetData(11, ref stream);

            var objectives = new List<FeaObjective>();
            foreach (var branch in objTree.Branches)
            {
                foreach (var goo in branch)
                {
                    if (goo is GH_ObjectWrapper wrapper && wrapper.Value is FeaObjective obj)
                        objectives.Add(obj);
                    else if (goo is FeaObjective directObj)
                        objectives.Add(directObj);
                }
            }

            if (objectives.Count == 0)
                AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "No objectives provided.");

            DA.SetData(0, new FeaOptimizationConfig
            {
                Objectives = objectives,
                OptimizeNodePositions = optPos,
                OptimizeAreas = optAreas,
                OptimizeSupportPositions = optSup,
                MaxIterations = maxIter,
                AbsTol = absTol,
                RelTol = relTol,
                BarrierWeight = bw,
                BarrierSharpness = bs,
                ReportFrequency = reportFreq,
                Run = run,
                StreamPreview = stream,
            });
        }

        protected override Bitmap Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-000000000006");
    }
}
