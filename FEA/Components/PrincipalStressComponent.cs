using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Runtime.InteropServices;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Types;
using Rhino.Display;
using Rhino.Geometry;

namespace Ariadne.FEA.Components
{
    /// <summary>
    /// P/Invoke declarations for the SPR recovery FFI in theseus.dll.
    /// </summary>
    internal static class SolidSprInterop
    {
        const string DLL = "theseus";

        [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
        public static extern int theseus_solid_spr_recover(
            nuint num_nodes,
            nuint num_elements,
            double[] node_positions,
            int[] elements,
            int[] num_nodes_per_element,
            double[] element_stresses,
            double[] nodal_stresses,
            double[] element_errors,
            double[] principal_values,
            double[] principal_vectors,
            double[] nodal_von_mises);
    }

    public class PrincipalStressComponent : GH_Component
    {
        private FeaResult? _cachedResult;
        private double[]? _cachedS1, _cachedS2, _cachedS3;
        private Vector3d[]? _cachedV1, _cachedV2, _cachedV3;
        private double[]? _cachedVm;
        private Point3d[]? _cachedNodes;
        private double _arrowScale = 1.0;
        private bool _showArrows = true;

        public PrincipalStressComponent()
            : base("Principal Stresses", "σ-Princ",
                "SPR-recovered principal stresses and directions for solid or shell FEA results",
                "Theseus-FEA", "Utilities")
        { }

        protected override void RegisterInputParams(GH_InputParamManager pManager)
        {
            pManager.AddGenericParameter("Results", "R", "Solid or shell FEA result", GH_ParamAccess.item);
            pManager.AddNumberParameter("Arrow Scale", "AS", "Scale factor for principal direction arrows", GH_ParamAccess.item, 1.0);
            pManager.AddBooleanParameter("Show Arrows", "Show", "Display principal direction arrows in viewport", GH_ParamAccess.item, true);
        }

        protected override void RegisterOutputParams(GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("S1", "S1", "First principal stress per node (most tensile)", GH_ParamAccess.list);
            pManager.AddNumberParameter("S2", "S2", "Second principal stress per node", GH_ParamAccess.list);
            pManager.AddNumberParameter("S3", "S3", "Third principal stress per node (most compressive)", GH_ParamAccess.list);
            pManager.AddVectorParameter("V1", "V1", "First principal direction per node", GH_ParamAccess.list);
            pManager.AddVectorParameter("V2", "V2", "Second principal direction per node", GH_ParamAccess.list);
            pManager.AddVectorParameter("V3", "V3", "Third principal direction per node", GH_ParamAccess.list);
            pManager.AddNumberParameter("VM", "VM", "Von Mises stress per node (SPR-recovered)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Err", "Err", "SPR error estimator per element", GH_ParamAccess.list);
            pManager.AddNumberParameter("Tensor", "T", "Flat 6-component recovered stress tensor per node [xx,yy,zz,xy,yz,xz]", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            object? input = null;
            if (!DA.GetData(0, ref input)) return;
            if (input is GH_ObjectWrapper wrapper)
                input = wrapper.Value;
            if (input is not FeaResult result) return;

            DA.GetData(1, ref _arrowScale);
            DA.GetData(2, ref _showArrows);

            if (result.Model.HasSolidElements && result.SolidStresses != null)
                SolveSolid(DA, result);
            else if (result.Model.HasShellElements && result.ShellTopStresses != null)
                SolveShell(DA, result);
            else
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                    "Results must be a solid or shell FeaResult with stress data.");
            }
        }

        private void SolveSolid(IGH_DataAccess DA, FeaResult result)
        {
            var mesh = result.Model.SolidMesh!;
            int nn = mesh.Nn;
            int ne = mesh.Ne;

            var doc = Rhino.RhinoDoc.ActiveDoc;
            double toMeters = doc != null
                ? Rhino.RhinoMath.UnitScale(doc.ModelUnitSystem, Rhino.UnitSystem.Meters)
                : 1.0;

            var nodePos = new double[nn * 3];
            for (int i = 0; i < nn; i++)
            {
                var pt = mesh.TetNodes[i];
                nodePos[i * 3] = pt.X * toMeters;
                nodePos[i * 3 + 1] = pt.Y * toMeters;
                nodePos[i * 3 + 2] = pt.Z * toMeters;
            }

            var elements = new int[mesh.Elements.Sum(e => e.Length)];
            var numNodesPerElement = new int[ne];
            int elemOffset = 0;
            for (int e = 0; e < ne; e++)
            {
                var elem = mesh.Elements[e];
                numNodesPerElement[e] = elem.Length;
                for (int j = 0; j < elem.Length; j++)
                    elements[elemOffset++] = elem[j];
            }

            int ngp = Theseus.Interop.SolidConstants.NumGP;
            var elemStresses = new double[ne * ngp * 6];
            for (int e = 0; e < ne; e++)
                for (int gp = 0; gp < ngp; gp++)
                    for (int c = 0; c < 6; c++)
                        elemStresses[e * ngp * 6 + gp * 6 + c] = result.SolidStresses![e, c];

            var nodalStresses = new double[nn * 6];
            var elementErrors = new double[ne];
            var principalValues = new double[nn * 3];
            var principalVectors = new double[nn * 9];
            var nodalVonMises = new double[nn];

            int rc = SolidSprInterop.theseus_solid_spr_recover(
                (nuint)nn, (nuint)ne,
                nodePos, elements, numNodesPerElement, elemStresses,
                nodalStresses, elementErrors,
                principalValues, principalVectors, nodalVonMises);

            if (rc != 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                    $"Solid SPR recovery failed (code {rc}).");
                return;
            }

            UnpackAndOutput(DA, result, nn, principalValues, principalVectors,
                nodalVonMises, elementErrors, nodalStresses, mesh.TetNodes.ToArray());
        }

        private void SolveShell(IGH_DataAccess DA, FeaResult result)
        {
            var model = result.Model;
            var shellData = model.ToShellData();
            int nn = shellData.NumNodes;
            int ne = shellData.NumElements;

            var memStresses = new double[ne * 6];
            var topStresses = result.ShellTopStresses!;
            var botStresses = result.ShellBottomStresses!;

            // Build membrane stresses from average of top and bottom
            for (int e = 0; e < ne; e++)
                for (int c = 0; c < 6; c++)
                    memStresses[e * 6 + c] = (topStresses[e * 6 + c] + botStresses[e * 6 + c]) / 2.0;

            var nodalMembrane = new double[nn * 6];
            var nodalTop = new double[nn * 6];
            var nodalBottom = new double[nn * 6];
            var elementErrors = new double[ne];
            var principalValuesTop = new double[nn * 3];
            var principalValuesBot = new double[nn * 3];
            var principalVectorsTop = new double[nn * 9];
            var principalVectorsBot = new double[nn * 9];
            var vonMisesTop = new double[nn];
            var vonMisesBot = new double[nn];

            int rc = Theseus.Interop.ShellInterop.theseus_shell_spr_recover(
                (nuint)nn, (nuint)ne,
                shellData.NodePositions, shellData.Elements, shellData.NumNodesPerElement,
                memStresses, topStresses, botStresses,
                nodalMembrane, nodalTop, nodalBottom,
                elementErrors,
                principalValuesTop, principalValuesBot,
                principalVectorsTop, principalVectorsBot,
                vonMisesTop, vonMisesBot);

            if (rc != 0)
            {
                AddRuntimeMessage(GH_RuntimeMessageLevel.Error,
                    $"Shell SPR recovery failed (code {rc}).");
                return;
            }

            // Use top-fiber principals for visualization (typically the critical fiber)
            // VM envelope = max(top, bottom) per node
            var vmEnvelope = new double[nn];
            for (int i = 0; i < nn; i++)
                vmEnvelope[i] = Math.Max(vonMisesTop[i], vonMisesBot[i]);

            var shellNodes = new Point3d[nn];
            for (int i = 0; i < nn; i++)
            {
                var v = model.ShellMesh!.Vertices[i];
                shellNodes[i] = new Point3d(v.X, v.Y, v.Z);
            }

            UnpackAndOutput(DA, result, nn, principalValuesTop, principalVectorsTop,
                vmEnvelope, elementErrors, nodalTop, shellNodes);
        }

        private void UnpackAndOutput(
            IGH_DataAccess DA, FeaResult result, int nn,
            double[] principalValues, double[] principalVectors,
            double[] nodalVonMises, double[] elementErrors,
            double[] nodalStresses, Point3d[] nodePositions)
        {
            var s1 = new double[nn];
            var s2 = new double[nn];
            var s3 = new double[nn];
            var v1 = new Vector3d[nn];
            var v2 = new Vector3d[nn];
            var v3 = new Vector3d[nn];

            for (int i = 0; i < nn; i++)
            {
                s1[i] = principalValues[i * 3];
                s2[i] = principalValues[i * 3 + 1];
                s3[i] = principalValues[i * 3 + 2];

                v1[i] = new Vector3d(
                    principalVectors[i * 9], principalVectors[i * 9 + 1], principalVectors[i * 9 + 2]);
                v2[i] = new Vector3d(
                    principalVectors[i * 9 + 3], principalVectors[i * 9 + 4], principalVectors[i * 9 + 5]);
                v3[i] = new Vector3d(
                    principalVectors[i * 9 + 6], principalVectors[i * 9 + 7], principalVectors[i * 9 + 8]);
            }

            _cachedResult = result;
            _cachedS1 = s1;
            _cachedS2 = s2;
            _cachedS3 = s3;
            _cachedV1 = v1;
            _cachedV2 = v2;
            _cachedV3 = v3;
            _cachedVm = nodalVonMises;
            _cachedNodes = nodePositions;

            DA.SetDataList(0, s1);
            DA.SetDataList(1, s2);
            DA.SetDataList(2, s3);
            DA.SetDataList(3, v1);
            DA.SetDataList(4, v2);
            DA.SetDataList(5, v3);
            DA.SetDataList(6, nodalVonMises);
            DA.SetDataList(7, elementErrors);
            DA.SetDataList(8, nodalStresses);
        }

        public override void DrawViewportWires(IGH_PreviewArgs args)
        {
            if (!_showArrows || _cachedResult == null || _cachedNodes == null) return;

            int nn = _cachedNodes.Length;
            if (_cachedS1 == null || _cachedS1.Length != nn) return;

            if (_cachedS2 == null || _cachedS3 == null ||
                _cachedV1 == null || _cachedV2 == null || _cachedV3 == null ||
                _cachedNodes == null) return;

            double maxStress = 1e-30;
            for (int i = 0; i < nn; i++)
            {
                double m = Math.Max(Math.Abs(_cachedS1![i]), Math.Max(Math.Abs(_cachedS2[i]), Math.Abs(_cachedS3[i])));
                if (m > maxStress) maxStress = m;
            }

            double invMax = _arrowScale / maxStress;

            for (int i = 0; i < nn; i++)
            {
                var pt = _cachedNodes[i];
                DrawPrincipalArrow(args, pt, _cachedV1[i], _cachedS1![i], invMax);
                DrawPrincipalArrow(args, pt, _cachedV2[i], _cachedS2[i], invMax);
                DrawPrincipalArrow(args, pt, _cachedV3[i], _cachedS3[i], invMax);
            }
        }

        private static void DrawPrincipalArrow(IGH_PreviewArgs args, Point3d pt, Vector3d dir, double stress, double scale)
        {
            double len = Math.Abs(stress) * scale;
            if (len < 1e-12) return;

            var color = stress >= 0 ? Color.Red : Color.Blue;
            var tip = pt + dir * len;
            args.Display.DrawArrow(new Line(pt, tip), color);
        }

        public override BoundingBox ClippingBox
        {
            get
            {
                if (_cachedNodes == null || _cachedNodes.Length == 0) return BoundingBox.Empty;
                var bb = new BoundingBox(_cachedNodes);
                bb.Inflate(_arrowScale);
                return bb;
            }
        }

        protected override Bitmap? Icon => null;
        public override Guid ComponentGuid => new Guid("A1B2C3D4-FEA1-4000-8001-100000000040");
    }
}
