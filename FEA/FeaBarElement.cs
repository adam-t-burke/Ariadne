namespace Ariadne.FEA;

using Ariadne.Graphs;

/// <summary>
/// A bar/truss FEA element that extends <see cref="Edge"/> with material, section,
/// and result properties. Since it inherits from Edge (which inherits GH_Curve),
/// it can be passed through Grasshopper wires and visualized natively.
/// Existing <see cref="FeaSolverContext.ResolveEdgeIndices"/> works unchanged
/// because FeaBarElement IS an Edge.
/// </summary>
public class FeaBarElement : Edge
{
    public int MaterialIndex { get; set; }
    public int SectionIndex { get; set; }
    public FeaMaterial Material { get; set; } = null!;
    public FeaSection Section { get; set; } = null!;
    public double AxialForce { get; set; }
    public double Stress { get; set; }
    public double Strain { get; set; }
    public double Utilization { get; set; }

    public FeaBarElement() { }

    public FeaBarElement(Edge edge, int matIdx, int secIdx,
        FeaMaterial material, FeaSection section)
    {
        Start = edge.Start;
        End = edge.End;
        Q = edge.Q;
        Value = edge.Value;
        ReferenceID = edge.ReferenceID;
        MaterialIndex = matIdx;
        SectionIndex = secIdx;
        Material = material;
        Section = section;
    }

    /// <summary>
    /// Create a FeaBarElement from an edge with result data populated.
    /// </summary>
    public static FeaBarElement FromEdgeWithResults(
        Edge edge, int matIdx, int secIdx,
        FeaMaterial material, FeaSection section,
        double axialForce, double stress, double strain, double utilization)
    {
        return new FeaBarElement(edge, matIdx, secIdx, material, section)
        {
            AxialForce = axialForce,
            Stress = stress,
            Strain = strain,
            Utilization = utilization,
        };
    }
}
