namespace Ariadne.FEA;

using System.Collections.Generic;
using System.Linq;
using Ariadne.Graphs;

public enum BeamFormulation
{
    EulerBernoulli = 0,
    Timoshenko = 1,
}

/// <summary>
/// Linear element network: wraps a <see cref="Graph"/> with per-edge material
/// and section assignments, a bending flag (truss vs frame), and a beam
/// formulation choice (Euler-Bernoulli vs Timoshenko).
/// </summary>
public class BarNetwork
{
    public Graph Graph { get; }
    public bool Bending { get; }
    public BeamFormulation Formulation { get; }
    public List<FeaMaterial> Materials { get; }
    public List<FeaSection> Sections { get; }
    public int[] MaterialAssignment { get; }
    public int[] SectionAssignment { get; }

    public int Nn => Graph.Nn;
    public int Ne => Graph.Ne;

    public BarNetwork(
        Graph graph, bool bending, BeamFormulation formulation,
        List<FeaMaterial> materials, List<FeaSection> sections,
        int[] matAssign, int[] secAssign)
    {
        Graph = graph;
        Bending = bending;
        Formulation = formulation;
        Materials = materials;
        Sections = sections;

        int ne = graph.Ne;
        MaterialAssignment = ExpandAssignment(matAssign, ne, materials.Count);
        SectionAssignment = ExpandAssignment(secAssign, ne, sections.Count);
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

    public override string ToString()
    {
        string type = Bending ? $"Frame ({Formulation})" : "Truss";
        return $"BarNetwork ({type}, {Nn} nodes, {Ne} edges)";
    }
}
