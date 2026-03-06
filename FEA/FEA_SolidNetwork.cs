namespace Ariadne.FEA;

using System.Collections.Generic;
using System.Linq;
using Rhino.Geometry;

/// <summary>
/// Full physics model for solid FEA: mesh + materials + supports + assignments.
/// Uses <see cref="FeaSupport"/> (point-based) with node matching at construction time,
/// matching the bar FEA pattern. Loads are not stored here -- they are a solve-time input.
/// </summary>
public class FEA_SolidNetwork
{
    public SolidMesh Mesh { get; set; }
    public List<FeaMaterial> Materials { get; set; } = [FeaMaterial.Steel()];
    public List<FeaSupport> Supports { get; set; } = [];

    public int[] MaterialAssignment { get; set; } = [];

    public bool Valid { get; private set; }
    public string ValidationMessage { get; private set; } = "";

    public List<int> SupportNodeIndices { get; private set; } = [];
    public Dictionary<int, int> SupportToNodeMap { get; private set; } = [];
    public double Tolerance { get; set; } = 0.001;

    public FEA_SolidNetwork() { Mesh = new SolidMesh(); }

    public FEA_SolidNetwork(SolidMesh mesh, List<FeaMaterial> materials,
        int[] matAssign, List<FeaSupport> supports, double tolerance = 0.001)
    {
        Mesh = mesh;
        Materials = materials;
        Supports = supports;
        Tolerance = tolerance;
        MaterialAssignment = ExpandAssignment(matAssign, mesh.Ne, materials.Count);

        MatchSupportsToNodes();
        Validate();
    }

    public FEA_SolidNetwork(FEA_SolidNetwork other)
    {
        Mesh = new SolidMesh(other.Mesh);
        Materials = new List<FeaMaterial>(other.Materials);
        Supports = new List<FeaSupport>(other.Supports);
        MaterialAssignment = (int[])other.MaterialAssignment.Clone();
        SupportNodeIndices = new List<int>(other.SupportNodeIndices);
        SupportToNodeMap = new Dictionary<int, int>(other.SupportToNodeMap);
        Tolerance = other.Tolerance;
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
            int bestNode = sup.NodeIndex ?? -1;
            double bestDist = 0.0;

            if (!sup.NodeIndex.HasValue)
            {
                bestDist = double.MaxValue;
                for (int n = 0; n < Mesh.Nn; n++)
                {
                    double dist = Mesh.TetNodes[n].DistanceToSquared(sup.Location);
                    if (dist < bestDist)
                    {
                        bestDist = dist;
                        bestNode = n;
                    }
                }
            }

            if (bestNode >= 0 && (sup.NodeIndex.HasValue || bestDist <= tol2))
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

        if (Mesh.Nn == 0) { Valid = false; ValidationMessage = "Mesh has no nodes."; return; }
        if (Mesh.Ne == 0) { Valid = false; ValidationMessage = "Mesh has no elements."; return; }
        if (Supports.Count == 0) { Valid = false; ValidationMessage = "At least one support is required."; return; }

        if (SupportNodeIndices.Count != Supports.Count)
        {
            Valid = false;
            ValidationMessage = $"Only {SupportNodeIndices.Count} of {Supports.Count} supports matched to nodes within tolerance {Tolerance}.";
            return;
        }

        if (SupportNodeIndices.Distinct().Count() != SupportNodeIndices.Count)
        {
            Valid = false;
            ValidationMessage = "Multiple supports mapped to the same node. Reduce tolerance, refine the mesh, or remove duplicate supports.";
            return;
        }

        if (Materials.Count == 0) { Valid = false; ValidationMessage = "At least one material is required."; return; }

        for (int e = 0; e < Mesh.Ne; e++)
        {
            if (MaterialAssignment[e] >= Materials.Count)
            {
                Valid = false;
                ValidationMessage = $"Element {e} references material index {MaterialAssignment[e]} but only {Materials.Count} materials exist.";
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
