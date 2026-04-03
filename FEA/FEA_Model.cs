namespace Ariadne.FEA;

using System;
using System.Collections.Generic;
using System.Linq;
using Ariadne.Graphs;
using Rhino.Geometry;

public class FEA_Model
{
    public Graph? BarGraph { get; set; }
    public BarNetwork? BarNetworkData { get; set; }
    public SolidMesh? SolidMesh { get; set; }
    public Mesh? ShellMesh { get; set; }
    public ShellMesh? ShellMeshData { get; set; }
    public double[] ShellThickness { get; set; } = [];

    public List<FeaMaterial> Materials { get; set; } = [FeaMaterial.Steel()];
    public List<FeaSection> Sections { get; set; } = [new FeaSection()];
    public List<FeaSupport> Supports { get; set; } = [];
    public List<FeaSolidSupportRegion> SolidSupportRegions { get; set; } = [];
    public List<FeaLoad> Loads { get; set; } = [];
    public int[] MaterialAssignment { get; set; } = [];
    public int[] SectionAssignment { get; set; } = [];
    public double Tolerance { get; set; } = 0.001;
    public bool Valid { get; private set; }
    public string ValidationMessage { get; private set; } = "";

    public List<int> SupportNodeIndices { get; private set; } = [];
    public Dictionary<int, int> SupportToNodeMap { get; private set; } = [];
    private List<FeaSupport> ExpandedSupports { get; set; } = [];

    public bool HasBarElements => BarGraph != null && BarGraph.Ne > 0;
    public bool HasSolidElements => SolidMesh != null && SolidMesh.Ne > 0;
    public bool HasShellElements => ShellMesh != null && ShellMesh.Faces.Count > 0;
    public bool IsCoupled => (HasBarElements ? 1 : 0) + (HasSolidElements ? 1 : 0) + (HasShellElements ? 1 : 0) > 1;

    public int TotalNodes => 
        (HasBarElements ? BarGraph!.Nn : 0) + 
        (HasSolidElements ? SolidMesh!.Nn : 0) +
        (HasShellElements ? ShellMesh!.Vertices.Count : 0);

    public int TotalElements => 
        (HasBarElements ? BarGraph!.Ne : 0) + 
        (HasSolidElements ? SolidMesh!.Ne : 0) +
        (HasShellElements ? ShellMesh!.Faces.Count : 0);

    public FEA_Model() { }

    public FEA_Model(Graph graph, List<FeaSupport> supports,
        List<FeaMaterial> materials, List<FeaSection> sections,
        int[] matAssign, int[] secAssign, double tolerance, List<FeaLoad>? loads = null, List<FeaSolidSupportRegion>? solidSupportRegions = null)
    {
        BarGraph = graph;
        Supports = supports;
        SolidSupportRegions = solidSupportRegions != null ? new List<FeaSolidSupportRegion>(solidSupportRegions) : [];
        Materials = materials;
        Sections = sections;
        Loads = loads != null ? new List<FeaLoad>(loads) : [];
        Tolerance = tolerance;

        int ne = graph.Ne;
        MaterialAssignment = ExpandAssignment(matAssign, ne, materials.Count);
        SectionAssignment = ExpandAssignment(secAssign, ne, sections.Count);

        MatchSupportsToNodes();
        Validate();
    }

    public FEA_Model(Mesh shellMesh, double[] shellThickness, List<FeaSupport> supports,
        List<FeaMaterial> materials, int[] matAssign, double tolerance, List<FeaLoad>? loads = null)
    {
        ShellMesh = shellMesh;
        ShellThickness = shellThickness;
        Supports = supports;
        Materials = materials;
        Loads = loads != null ? new List<FeaLoad>(loads) : [];
        Tolerance = tolerance;

        int ne = shellMesh.Faces.Count;
        MaterialAssignment = ExpandAssignment(matAssign, ne, materials.Count);

        MatchSupportsToNodes();
        Validate();
    }

    public FEA_Model(FEA_Model other)
    {
        BarGraph = other.BarGraph != null ? new Graph(other.BarGraph) : null;
        BarNetworkData = other.BarNetworkData;
        SolidMesh = other.SolidMesh != null ? new SolidMesh(other.SolidMesh) : null;
        ShellMesh = other.ShellMesh?.DuplicateMesh();
        ShellMeshData = other.ShellMeshData;
        ShellThickness = other.ShellThickness.Length > 0 ? (double[])other.ShellThickness.Clone() : [];
        Supports = new List<FeaSupport>(other.Supports);
        SolidSupportRegions = other.SolidSupportRegions
            .Select(region => new FeaSolidSupportRegion
            {
                NodeIndices = new List<int>(region.NodeIndices),
                FaceIndices = new List<int>(region.FaceIndices),
                FixX = region.FixX,
                FixY = region.FixY,
                FixZ = region.FixZ,
            })
            .ToList();
        Loads = new List<FeaLoad>(other.Loads);
        Materials = new List<FeaMaterial>(other.Materials);
        Sections = new List<FeaSection>(other.Sections);
        MaterialAssignment = (int[])other.MaterialAssignment.Clone();
        SectionAssignment = (int[])other.SectionAssignment.Clone();
        Tolerance = other.Tolerance;
        SupportNodeIndices = new List<int>(other.SupportNodeIndices);
        SupportToNodeMap = new Dictionary<int, int>(other.SupportToNodeMap);
        ExpandedSupports = new List<FeaSupport>(other.ExpandedSupports);
        Valid = other.Valid;
        ValidationMessage = other.ValidationMessage;
    }

    public void Initialize()
    {
        MatchSupportsToNodes();
        Validate();
    }

    private void MatchSupportsToNodes()
    {
        ExpandedSupports = BuildExpandedSupports();
        SupportNodeIndices.Clear();
        SupportToNodeMap.Clear();
        double tol2 = Tolerance * Tolerance;

        var nodePositions = new List<Point3d>();
        if (HasBarElements)
            nodePositions.AddRange(BarGraph!.Nodes.Select(n => n.Value));
        if (HasSolidElements)
            nodePositions.AddRange(SolidMesh!.TetNodes);
        if (HasShellElements)
            nodePositions.AddRange(ShellMesh!.Vertices.Select(v => (Point3d)v));

        for (int s = 0; s < ExpandedSupports.Count; s++)
        {
            var sup = ExpandedSupports[s];
            int bestNode = sup.NodeIndex ?? -1;
            double bestDist = 0.0;

            if (!sup.NodeIndex.HasValue)
            {
                bestDist = double.MaxValue;
                for (int n = 0; n < nodePositions.Count; n++)
                {
                    double dist = nodePositions[n].DistanceToSquared(sup.Location);
                    if (dist < bestDist)
                    {
                        bestDist = dist;
                        bestNode = n;
                    }
                }
            }

            if (bestNode >= 0 &&
                (sup.NodeIndex.HasValue || bestDist <= tol2))
            {
                SupportNodeIndices.Add(bestNode);
                SupportToNodeMap[s] = bestNode;
            }
        }
    }

    private void Validate()
    {
        Valid = true;
        ValidationMessage = "";

        if (!HasBarElements && !HasSolidElements && !HasShellElements)
        {
            Valid = false;
            ValidationMessage = "Model has no elements (provide Edges for bar, SolidMesh for solid, or Mesh for shell).";
            return;
        }

        if (Supports.Count + SolidSupportRegions.Count == 0)
        {
            Valid = false;
            ValidationMessage = "At least one support is required.";
            return;
        }

        if (!HasSolidElements && SolidSupportRegions.Count > 0)
        {
            Valid = false;
            ValidationMessage = "Solid support regions require a solid mesh.";
            return;
        }

        if (HasSolidElements)
        {
            for (int i = 0; i < SolidSupportRegions.Count; i++)
            {
                var region = SolidSupportRegions[i];
                if (!region.FixX && !region.FixY && !region.FixZ)
                {
                    Valid = false;
                    ValidationMessage = $"Solid support region {i} has no constrained DOFs.";
                    return;
                }

                foreach (int nodeIndex in region.NodeIndices)
                {
                    if (nodeIndex < 0 || nodeIndex >= SolidMesh!.Nn)
                    {
                        Valid = false;
                        ValidationMessage = $"Solid support region {i} references invalid solid node index {nodeIndex}.";
                        return;
                    }
                }
            }
        }

        if (SupportNodeIndices.Count != ExpandedSupports.Count)
        {
            Valid = false;
            ValidationMessage = $"Only {SupportNodeIndices.Count} of {ExpandedSupports.Count} supports matched to nodes within tolerance {Tolerance}.";
            return;
        }

        if (SupportNodeIndices.Distinct().Count() != SupportNodeIndices.Count)
        {
            Valid = false;
            ValidationMessage = "Multiple supports mapped to the same node. Reduce tolerance, refine the mesh, or remove duplicate supports.";
            return;
        }

        if (Materials.Count == 0)
        {
            Valid = false;
            ValidationMessage = "At least one material is required.";
            return;
        }

        int ne = HasBarElements ? BarGraph!.Ne
               : HasSolidElements ? SolidMesh!.Ne
               : ShellMesh!.Faces.Count;
        for (int e = 0; e < ne && e < MaterialAssignment.Length; e++)
        {
            if (MaterialAssignment[e] >= Materials.Count)
            {
                Valid = false;
                ValidationMessage = $"Element {e} references material index {MaterialAssignment[e]} but only {Materials.Count} materials exist.";
                return;
            }
        }

        if (HasBarElements)
        {
            if (Sections.Count == 0)
            {
                Valid = false;
                ValidationMessage = "At least one section is required for bar elements.";
                return;
            }
            for (int e = 0; e < BarGraph!.Ne && e < SectionAssignment.Length; e++)
            {
                if (SectionAssignment[e] >= Sections.Count)
                {
                    Valid = false;
                    ValidationMessage = $"Element {e} references section index {SectionAssignment[e]} but only {Sections.Count} sections exist.";
                    return;
                }
            }
        }

        if (HasSolidElements || HasShellElements)
        {
            foreach (var mat in Materials)
            {
                if (mat.Nu <= -1.0 || mat.Nu >= 0.5)
                {
                    Valid = false;
                    string elemType = HasSolidElements ? "solid" : "shell";
                    ValidationMessage = $"Poisson's ratio {mat.Nu} is out of valid range (-1, 0.5) for {elemType} elements.";
                    return;
                }
            }
        }

        for (int i = 0; i < Loads.Count; i++)
        {
            if (Loads[i].NodeIndex < 0 || Loads[i].NodeIndex >= TotalNodes)
            {
                Valid = false;
                ValidationMessage = $"Load {i} references node index {Loads[i].NodeIndex}, but the model only has {TotalNodes} nodes.";
                return;
            }
        }
    }

    public List<object> GetSupportDefinitions()
    {
        var result = new List<object>(Supports.Count + SolidSupportRegions.Count);
        result.AddRange(Supports);
        result.AddRange(SolidSupportRegions);
        return result;
    }

    public List<FeaSupport> GetExpandedSupports() => new(ExpandedSupports);

    public Point3d GetNodePosition(int nodeIndex)
    {
        if (nodeIndex < 0 || nodeIndex >= TotalNodes)
            throw new ArgumentOutOfRangeException(nameof(nodeIndex));

        int offset = 0;

        if (HasBarElements)
        {
            if (nodeIndex < offset + BarGraph!.Nn)
                return BarGraph.Nodes[nodeIndex - offset].Value;
            offset += BarGraph.Nn;
        }

        if (HasSolidElements)
        {
            if (nodeIndex < offset + SolidMesh!.Nn)
                return SolidMesh.TetNodes[nodeIndex - offset];
            offset += SolidMesh.Nn;
        }

        if (HasShellElements)
            return (Point3d)ShellMesh!.Vertices[nodeIndex - offset];

        throw new InvalidOperationException($"Node index {nodeIndex} could not be resolved.");
    }

    public Node CreateNode(int nodeIndex)
    {
        return new Node
        {
            Value = GetNodePosition(nodeIndex),
            Index = nodeIndex,
        };
    }

    public bool[] GetConstrainedDofs(int nodeIndex)
    {
        var constrained = new bool[3];

        foreach (var kvp in SupportToNodeMap)
        {
            if (kvp.Value != nodeIndex)
                continue;

            var support = ExpandedSupports[kvp.Key];
            constrained[0] |= support.FixX;
            constrained[1] |= support.FixY;
            constrained[2] |= support.FixZ;
        }

        return constrained;
    }

    /// <summary>
    /// Creates an FEA_Network from the bar-element portion of this model.
    /// Requires <see cref="HasBarElements"/> to be true.
    /// </summary>
    public FEA_Network ToBarNetwork()
    {
        if (!HasBarElements)
            throw new InvalidOperationException("Model has no bar elements.");

        return new FEA_Network(
            BarGraph!,
            Supports,
            Materials,
            Sections,
            MaterialAssignment,
            SectionAssignment,
            Tolerance);
    }

    /// <summary>
    /// Creates a ShellSolverData for the shell portion of this model.
    /// Uses pre-triangulated ShellMeshData when available, otherwise triangulates on the fly.
    /// </summary>
    internal ShellSolverData ToShellData()
    {
        if (!HasShellElements)
            throw new InvalidOperationException("Model has no shell elements.");

        int nn = ShellMesh!.Vertices.Count;

        var nodePos = new double[nn * 3];
        for (int i = 0; i < nn; i++)
        {
            var v = ShellMesh.Vertices[i];
            nodePos[i * 3] = v.X;
            nodePos[i * 3 + 1] = v.Y;
            nodePos[i * 3 + 2] = v.Z;
        }

        var elements = new List<int>();
        var triMatAssignment = new List<int>();
        var triSecAssignment = new List<int>();

        if (ShellMeshData != null)
        {
            foreach (var tri in ShellMeshData.Elements)
            {
                elements.Add(tri[0]);
                elements.Add(tri[1]);
                elements.Add(tri[2]);
            }
            for (int e = 0; e < ShellMeshData.Ne; e++)
            {
                int origFace = ShellMeshData.TriToFaceMap[e];
                int matIdx = origFace < MaterialAssignment.Length ? MaterialAssignment[origFace] : 0;
                int secIdx = origFace < SectionAssignment.Length ? SectionAssignment[origFace] : 0;
                triMatAssignment.Add(matIdx);
                triSecAssignment.Add(secIdx);
            }
        }
        else
        {
            int faceCount = ShellMesh.Faces.Count;
            for (int i = 0; i < faceCount; i++)
            {
                var f = ShellMesh.Faces[i];
                int matIdx = i < MaterialAssignment.Length ? MaterialAssignment[i] : 0;
                int secIdx = i < SectionAssignment.Length ? SectionAssignment[i] : 0;

                elements.Add(f.A);
                elements.Add(f.B);
                elements.Add(f.C);
                triMatAssignment.Add(matIdx);
                triSecAssignment.Add(secIdx);

                if (f.IsQuad)
                {
                    elements.Add(f.A);
                    elements.Add(f.C);
                    elements.Add(f.D);
                    triMatAssignment.Add(matIdx);
                    triSecAssignment.Add(secIdx);
                }
            }
        }

        int ne = triMatAssignment.Count;
        var numNodesPerElement = new int[ne];
        for (int i = 0; i < ne; i++)
            numNodesPerElement[i] = 3;

        var matData = new double[Materials.Count * 4];
        for (int i = 0; i < Materials.Count; i++)
        {
            matData[i * 4] = Materials[i].E;
            matData[i * 4 + 1] = Materials[i].Nu;
            matData[i * 4 + 2] = Materials[i].Density;
            matData[i * 4 + 3] = Materials[i].YieldStress;
        }

        var expanded = BuildExpandedSupports();
        double tol2 = Tolerance * Tolerance;
        var nodePositions = ShellMesh.Vertices.Select(v => (Point3d)v).ToList();

        var shellSupports = new List<int>();
        foreach (var sup in expanded)
        {
            int bestNode = -1;
            double bestDist = double.MaxValue;
            for (int n = 0; n < nn; n++)
            {
                double dist = nodePositions[n].DistanceToSquared(sup.Location);
                if (dist < bestDist)
                {
                    bestDist = dist;
                    bestNode = n;
                }
            }
            if (bestNode >= 0 && bestDist <= tol2)
            {
                shellSupports.Add(bestNode);
                shellSupports.Add(sup.FixX ? 1 : 0);
                shellSupports.Add(sup.FixY ? 1 : 0);
                shellSupports.Add(sup.FixZ ? 1 : 0);
                shellSupports.Add(sup.FixRX ? 1 : 0);
                shellSupports.Add(sup.FixRY ? 1 : 0);
                shellSupports.Add(sup.FixRZ ? 1 : 0);
            }
        }

        var loadForces = new double[Loads.Count * 3];
        var loadMoments = new double[Loads.Count * 3];
        var loadNodes = new int[Loads.Count];
        for (int i = 0; i < Loads.Count; i++)
        {
            loadNodes[i] = Loads[i].NodeIndex;
            loadForces[i * 3] = Loads[i].Force.X;
            loadForces[i * 3 + 1] = Loads[i].Force.Y;
            loadForces[i * 3 + 2] = Loads[i].Force.Z;
            loadMoments[i * 3] = Loads[i].Moment.X;
            loadMoments[i * 3 + 1] = Loads[i].Moment.Y;
            loadMoments[i * 3 + 2] = Loads[i].Moment.Z;
        }

        var secData = new double[Sections.Count];
        for (int i = 0; i < Sections.Count; i++)
            secData[i] = Sections[i].Offset;

        var ep = new int[ne * 2];
        for (int e = 0; e < ne; e++)
        {
            ep[e * 2] = triMatAssignment[e];
            ep[e * 2 + 1] = triSecAssignment[e];
        }

        return new ShellSolverData(
            nn, ne, nodePos, ShellThickness, elements.ToArray(), numNodesPerElement,
            Materials.Count, matData, Sections.Count, secData, ep,
            shellSupports.Count / 7, shellSupports.ToArray(),
            Loads.Count, loadForces, loadMoments, loadNodes);
    }

    private List<FeaSupport> BuildExpandedSupports()
    {
        var expanded = new Dictionary<int, FeaSupport>();
        var floating = new List<FeaSupport>();

        foreach (var support in Supports)
        {
            if (support.NodeIndex.HasValue)
            {
                expanded[support.NodeIndex.Value] = MergeSupport(expanded.GetValueOrDefault(support.NodeIndex.Value), support);
            }
            else
            {
                floating.Add(support);
            }
        }

        if (HasSolidElements)
        {
            foreach (var region in SolidSupportRegions)
            {
                foreach (int nodeIndex in region.NodeIndices.Distinct())
                {
                    if (nodeIndex < 0 || nodeIndex >= SolidMesh!.Nn)
                        continue;

                    var support = new FeaSupport
                    {
                        Location = SolidMesh.TetNodes[nodeIndex],
                        NodeIndex = nodeIndex,
                        FixX = region.FixX,
                        FixY = region.FixY,
                        FixZ = region.FixZ,
                    };
                    expanded[nodeIndex] = MergeSupport(expanded.GetValueOrDefault(nodeIndex), support);
                }
            }
        }

        var result = new List<FeaSupport>(floating);
        result.AddRange(expanded.OrderBy(kvp => kvp.Key).Select(kvp => kvp.Value));
        return result;
    }

    private static FeaSupport MergeSupport(FeaSupport? existing, FeaSupport incoming)
    {
        if (existing == null)
            return incoming;

        return new FeaSupport
        {
            Location = incoming.Location,
            NodeIndex = incoming.NodeIndex ?? existing.NodeIndex,
            FixX = existing.FixX || incoming.FixX,
            FixY = existing.FixY || incoming.FixY,
            FixZ = existing.FixZ || incoming.FixZ,
            FixRX = existing.FixRX || incoming.FixRX,
            FixRY = existing.FixRY || incoming.FixRY,
            FixRZ = existing.FixRZ || incoming.FixRZ,
        };
    }

    internal static int[] ExpandAssignment(int[] input, int count, int maxIdx)
    {
        if (input == null || input.Length == 0)
            return Enumerable.Repeat(0, count).ToArray();

        var result = new int[count];
        for (int i = 0; i < count; i++)
            result[i] = i < input.Length ? input[i] : input[^1];
        return result;
    }
}
