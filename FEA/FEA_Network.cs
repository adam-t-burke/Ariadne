namespace Ariadne.FEA;

using System;
using System.Collections.Generic;
using System.Linq;
using Ariadne.Graphs;
using Rhino.Geometry;

/// <summary>
/// FEA network: graph topology + materials + sections + supports + assignments.
/// Analogous to <see cref="Ariadne.FDM.FDM_Network"/> for the FDM solver.
/// </summary>
public class FEA_Network
{
    public Graph Graph { get; set; }
    public List<FeaSupport> Supports { get; set; } = [];
    public List<FeaMaterial> Materials { get; set; } = [FeaMaterial.Steel()];
    public List<FeaSection> Sections { get; set; } = [new FeaSection()];
    public int[] MaterialAssignment { get; set; } = [];
    public int[] SectionAssignment { get; set; } = [];
    public double Tolerance { get; set; } = 0.001;
    public bool Valid { get; private set; }
    public string ValidationMessage { get; private set; } = "";

    /// <summary>Node indices that have at least one constrained DOF.</summary>
    public List<int> SupportNodeIndices { get; private set; } = [];

    /// <summary>Maps support index to graph node index.</summary>
    public Dictionary<int, int> SupportToNodeMap { get; private set; } = [];

    public FEA_Network() { Graph = new Graph(); }

    public FEA_Network(Graph graph, List<FeaSupport> supports,
        List<FeaMaterial> materials, List<FeaSection> sections,
        int[] matAssign, int[] secAssign, double tolerance)
    {
        Graph = graph;
        Supports = supports;
        Materials = materials;
        Sections = sections;
        Tolerance = tolerance;

        int ne = graph.Ne;
        MaterialAssignment = ExpandAssignment(matAssign, ne, materials.Count);
        SectionAssignment = ExpandAssignment(secAssign, ne, sections.Count);

        MatchSupportsToNodes();
        Validate();
    }

    public FEA_Network(FEA_Network other)
    {
        Graph = new Graph(other.Graph);
        Supports = new List<FeaSupport>(other.Supports);
        Materials = new List<FeaMaterial>(other.Materials);
        Sections = new List<FeaSection>(other.Sections);
        MaterialAssignment = (int[])other.MaterialAssignment.Clone();
        SectionAssignment = (int[])other.SectionAssignment.Clone();
        Tolerance = other.Tolerance;
        SupportNodeIndices = new List<int>(other.SupportNodeIndices);
        SupportToNodeMap = new Dictionary<int, int>(other.SupportToNodeMap);
        Valid = other.Valid;
        ValidationMessage = other.ValidationMessage;
    }

    private void MatchSupportsToNodes()
    {
        SupportNodeIndices.Clear();
        SupportToNodeMap.Clear();
        double tol2 = Tolerance * Tolerance;

        for (int s = 0; s < Supports.Count; s++)
        {
            var sup = Supports[s];
            int bestNode = -1;
            double bestDist = double.MaxValue;

            for (int n = 0; n < Graph.Nn; n++)
            {
                double dist = Graph.Nodes[n].Value.DistanceToSquared(sup.Location);
                if (dist < bestDist)
                {
                    bestDist = dist;
                    bestNode = n;
                }
            }

            if (bestNode >= 0 && bestDist <= tol2)
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

        if (Graph.Nn == 0)
        {
            Valid = false;
            ValidationMessage = "Network has no nodes.";
            return;
        }

        if (Supports.Count == 0)
        {
            Valid = false;
            ValidationMessage = "At least one support is required.";
            return;
        }

        if (SupportNodeIndices.Count != Supports.Count)
        {
            Valid = false;
            ValidationMessage = $"Only {SupportNodeIndices.Count} of {Supports.Count} supports matched to nodes within tolerance {Tolerance}.";
            return;
        }

        if (Materials.Count == 0)
        {
            Valid = false;
            ValidationMessage = "At least one material is required.";
            return;
        }

        if (Sections.Count == 0)
        {
            Valid = false;
            ValidationMessage = "At least one section is required.";
            return;
        }

        // Check that all material/section assignments are in range
        for (int e = 0; e < Graph.Ne; e++)
        {
            if (MaterialAssignment[e] >= Materials.Count)
            {
                Valid = false;
                ValidationMessage = $"Element {e} references material index {MaterialAssignment[e]} but only {Materials.Count} materials exist.";
                return;
            }
            if (SectionAssignment[e] >= Sections.Count)
            {
                Valid = false;
                ValidationMessage = $"Element {e} references section index {SectionAssignment[e]} but only {Sections.Count} sections exist.";
                return;
            }
        }
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
