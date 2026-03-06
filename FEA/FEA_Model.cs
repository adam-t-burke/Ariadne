namespace Ariadne.FEA;

using System;
using System.Collections.Generic;
using System.Linq;
using Ariadne.Graphs;
using Rhino.Geometry;

public class FEA_Model
{
    public Graph? BarGraph { get; set; }
    public SolidMesh? SolidMesh { get; set; }

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
    public bool IsCoupled => HasBarElements && HasSolidElements;

    public int TotalNodes => (HasBarElements ? BarGraph!.Nn : 0) + (HasSolidElements ? SolidMesh!.Nn : 0);
    public int TotalElements => (HasBarElements ? BarGraph!.Ne : 0) + (HasSolidElements ? SolidMesh!.Ne : 0);

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

    public FEA_Model(SolidMesh solidMesh, List<FeaSupport> supports,
        List<FeaMaterial> materials, int[] matAssign, double tolerance, List<FeaLoad>? loads = null, List<FeaSolidSupportRegion>? solidSupportRegions = null)
    {
        SolidMesh = solidMesh;
        Supports = supports;
        SolidSupportRegions = solidSupportRegions != null ? new List<FeaSolidSupportRegion>(solidSupportRegions) : [];
        Materials = materials;
        Loads = loads != null ? new List<FeaLoad>(loads) : [];
        Tolerance = tolerance;

        int ne = solidMesh.Ne;
        MaterialAssignment = ExpandAssignment(matAssign, ne, materials.Count);

        MatchSupportsToNodes();
        Validate();
    }

    public FEA_Model(FEA_Model other)
    {
        BarGraph = other.BarGraph != null ? new Graph(other.BarGraph) : null;
        SolidMesh = other.SolidMesh != null ? new SolidMesh(other.SolidMesh) : null;
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

    private void MatchSupportsToNodes()
    {
        ExpandedSupports = BuildExpandedSupports();
        SupportNodeIndices.Clear();
        SupportToNodeMap.Clear();
        double tol2 = Tolerance * Tolerance;

        var nodePositions = new List<Point3d>();
        if (HasBarElements)
            nodePositions.AddRange(BarGraph!.Nodes.Select(n => n.Value));
        else if (HasSolidElements)
            nodePositions.AddRange(SolidMesh!.TetNodes);

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

        if (!HasBarElements && !HasSolidElements)
        {
            Valid = false;
            ValidationMessage = "Model has no elements (provide Edges for bar or SolidMesh for solid).";
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

        int ne = HasBarElements ? BarGraph!.Ne : SolidMesh!.Ne;
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

        if (HasSolidElements)
        {
            foreach (var mat in Materials)
            {
                if (mat.Nu <= -1.0 || mat.Nu >= 0.5)
                {
                    Valid = false;
                    ValidationMessage = $"Poisson's ratio {mat.Nu} is out of valid range (-1, 0.5) for solid elements.";
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

        if (HasBarElements && nodeIndex < BarGraph!.Nn)
            return BarGraph.Nodes[nodeIndex].Value;

        int solidNodeIndex = HasBarElements ? nodeIndex - BarGraph!.Nn : nodeIndex;
        return SolidMesh!.TetNodes[solidNodeIndex];
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
    /// Creates an FEA_SolidNetwork from the solid-element portion of this model.
    /// Requires <see cref="HasSolidElements"/> to be true.
    /// Uses the same FeaSupport objects -- node matching happens in FEA_SolidNetwork constructor.
    /// </summary>
    public FEA_SolidNetwork ToSolidNetwork()
    {
        if (!HasSolidElements)
            throw new InvalidOperationException("Model has no solid elements.");

        return new FEA_SolidNetwork(
            SolidMesh!,
            Materials,
            MaterialAssignment,
            GetExpandedSupports(),
            Tolerance);
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
        };
    }

    private static int[] ExpandAssignment(int[] input, int count, int maxIdx)
    {
        if (input == null || input.Length == 0)
            return Enumerable.Repeat(0, count).ToArray();

        var result = new int[count];
        for (int i = 0; i < count; i++)
            result[i] = i < input.Length ? input[i] : input[^1];
        return result;
    }
}
