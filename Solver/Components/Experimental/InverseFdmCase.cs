using System;
using System.Collections.Generic;
using System.Linq;
using System.Text.Json.Serialization;
using Ariadne.FDM;
using Ariadne.Graphs;
using Grasshopper.Kernel.Data;
using Rhino.Geometry;

namespace Ariadne.Solver.Components.Experimental;

/// <summary>
/// One external inverse-FDM benchmark case as written by <c>bench/external/export_*.py</c>
/// and consumed by the Rust <c>warm_start_bench external</c> loader. Per-node arrays cover
/// every node (0-based); per-edge arrays are parallel to <see cref="Edges"/>.
/// Sign convention: <c>q &gt; 0</c> tension.
/// </summary>
public sealed class InverseFdmCase
{
    [JsonPropertyName("name")] public string Name { get; set; } = "";
    [JsonPropertyName("source")] public string? Source { get; set; }
    [JsonPropertyName("description")] public string? Description { get; set; }
    [JsonPropertyName("nodes")] public double[][] Nodes { get; set; } = [];
    [JsonPropertyName("edges")] public int[][] Edges { get; set; } = [];
    [JsonPropertyName("fixed")] public int[] Fixed { get; set; } = [];
    [JsonPropertyName("loads")] public double[][] Loads { get; set; } = [];
    [JsonPropertyName("q_ref")] public double[] QRef { get; set; } = [];
    [JsonPropertyName("target")] public double[][] Target { get; set; } = [];
    [JsonPropertyName("signs")] public double[]? Signs { get; set; }
    /// <summary>Loose box (the case's own, or the suite's ×/÷100 box); <c>null</c> = none.</summary>
    [JsonPropertyName("bounds")] public InverseFdmCaseBounds? Bounds { get; set; }
    /// <summary>Synthetic suite only: the snug ×/÷1 box around <c>q_ref</c>.</summary>
    [JsonPropertyName("bounds_snug")] public InverseFdmCaseBounds? BoundsSnug { get; set; }
    /// <summary>Designer target per node; <c>null</c> entries mean the node carries no goal.</summary>
    [JsonPropertyName("target_original")] public double[]?[]? TargetOriginal { get; set; }
    /// <summary>Synthetic suite only: perturbed full-node targets keyed by variant (<c>jit2pctd</c>, <c>bump10pctd</c>).</summary>
    [JsonPropertyName("target_variants")] public Dictionary<string, double[][]>? TargetVariants { get; set; }
    /// <summary>Characteristic length used to normalise geometric error (err / L).</summary>
    [JsonPropertyName("extent")] public double? Extent { get; set; }

    /// <summary>Mirrors the Rust <c>Case::validate</c> checks.</summary>
    public void Validate()
    {
        int n = Nodes.Length;
        int m = Edges.Length;
        if (n == 0) throw new InvalidOperationException($"Case '{Name}' has no nodes.");
        if (m == 0) throw new InvalidOperationException($"Case '{Name}' has no edges.");
        if (Loads.Length != n || Target.Length != n)
            throw new InvalidOperationException($"Case '{Name}': loads/target length != {n} nodes.");
        if (QRef.Length != m)
            throw new InvalidOperationException($"Case '{Name}': q_ref length != {m} edges.");
        if (Signs is not null && Signs.Length != m)
            throw new InvalidOperationException($"Case '{Name}': signs length != {m} edges.");
        if (Bounds is not null && (Bounds.Lo.Length != m || Bounds.Hi.Length != m))
            throw new InvalidOperationException($"Case '{Name}': bounds length != {m} edges.");
        if (BoundsSnug is not null && (BoundsSnug.Lo.Length != m || BoundsSnug.Hi.Length != m))
            throw new InvalidOperationException($"Case '{Name}': bounds_snug length != {m} edges.");
        if (TargetOriginal is not null && TargetOriginal.Length != n)
            throw new InvalidOperationException($"Case '{Name}': target_original length != {n} nodes.");
        if (Nodes.Any(p => p.Length != 3) || Loads.Any(p => p.Length != 3) || Target.Any(p => p.Length != 3))
            throw new InvalidOperationException($"Case '{Name}': nodes/loads/target entries must have 3 components.");
        if (TargetOriginal is not null && TargetOriginal.Any(p => p is not null && p.Length != 3))
            throw new InvalidOperationException($"Case '{Name}': target_original entries must have 3 components.");
        if (TargetVariants is not null)
        {
            foreach (var (key, variant) in TargetVariants)
            {
                if (variant.Length != n || variant.Any(p => p.Length != 3))
                    throw new InvalidOperationException($"Case '{Name}': target_variants['{key}'] must have {n} xyz entries.");
            }
        }
        if (Edges.Any(e => e.Length != 2 || e[0] < 0 || e[1] < 0 || e[0] >= n || e[1] >= n)
            || Fixed.Any(f => f < 0 || f >= n))
            throw new InvalidOperationException($"Case '{Name}': edge or support index out of range.");
    }
}

/// <summary>Per-edge box; <c>null</c> entries mean that side is unbounded.</summary>
public sealed class InverseFdmCaseBounds
{
    [JsonPropertyName("lo")] public double?[] Lo { get; set; } = [];
    [JsonPropertyName("hi")] public double?[] Hi { get; set; } = [];
}

/// <summary>
/// A case turned into Grasshopper-ready inputs. All free-node lists are in
/// <see cref="FDM_Network.Free"/> order; per-edge arrays follow <see cref="Graph.Edges"/>.
/// </summary>
public sealed class InverseFdmCaseNetwork
{
    public required string Name { get; init; }
    public required FDM_Network Network { get; init; }
    /// <summary>Target position of every free node.</summary>
    public required List<Point3d> Target { get; init; }
    /// <summary>Load on every free node (zeros for pure prestress cases).</summary>
    public required List<Vector3d> Loads { get; init; }
    /// <summary>Current position of every free node, parallel to <see cref="Loads"/>.</summary>
    public required List<Point3d> LoadNodes { get; init; }
    /// <summary>+1 tension, -1 compression, 0 free.</summary>
    public required int[] Signs { get; init; }
    /// <summary>Lower bound on q per edge; <c>null</c> = unbounded below.</summary>
    public required double?[] Lower { get; init; }
    /// <summary>Upper bound on q per edge; <c>null</c> = unbounded above.</summary>
    public required double?[] Upper { get; init; }
    /// <summary>Reference force densities that reproduce the case geometry (oracle, not seeded on the network).</summary>
    public required double[] QRef { get; init; }
    /// <summary>Target actually used: <see cref="InverseFdmCaseNetworkBuilder.ExactTarget"/>, <see cref="InverseFdmCaseNetworkBuilder.DesignerTarget"/>, or a variant key.</summary>
    public required string TargetKind { get; init; }
    /// <summary>Target kinds this case can serve.</summary>
    public required List<string> AvailableTargets { get; init; }
    /// <summary>True when the snug box was requested and the case carries one.</summary>
    public required bool UsedSnugBox { get; init; }
    public required List<string> Warnings { get; init; }

    public bool UsedDesignerTarget => TargetKind == InverseFdmCaseNetworkBuilder.DesignerTarget;

    /// <summary>
    /// Bound side as a dense per-edge array for a Grasshopper tree, or <c>null</c> when every
    /// entry is unbounded (emit an empty tree). Mixed sides fill unbounded entries with <paramref name="fill"/>.
    /// </summary>
    public static double[]? DenseBounds(double?[] side, double fill)
    {
        if (side.All(v => v is null)) return null;
        return side.Select(v => v ?? fill).ToArray();
    }
}

/// <summary>
/// Builds an <see cref="FDM_Network"/> from a case while preserving JSON node indices.
/// Anchors are assigned by index (not geometric matching) so coincident nodes are never merged,
/// and nodes are partitioned free-first exactly as <see cref="FDM_Network"/> does.
/// </summary>
public static class InverseFdmCaseNetworkBuilder
{
    /// <summary>Geometric tolerance stored on the network; nodes are matched by index, so this is only used downstream.</summary>
    public const double Tolerance = 1e-6;

    /// <summary>Target kind: the case's exact target (<c>target</c>, equal to <c>nodes</c> in every export).</summary>
    public const string ExactTarget = "exact";
    /// <summary>Target kind: the upstream example's own design target (<c>target_original</c>).</summary>
    public const string DesignerTarget = "designer";

    /// <summary>Target kinds a case can serve, in dropdown order: exact, designer (when complete), then variants.</summary>
    public static List<string> AvailableTargets(InverseFdmCase c)
    {
        var kinds = new List<string> { ExactTarget };
        if (c.TargetOriginal is not null)
        {
            var isFixed = new HashSet<int>(c.Fixed);
            bool complete = Enumerable.Range(0, c.Nodes.Length)
                .Where(i => !isFixed.Contains(i))
                .All(i => c.TargetOriginal[i] is not null);
            if (complete) kinds.Add(DesignerTarget);
        }
        if (c.TargetVariants is not null)
            kinds.AddRange(c.TargetVariants.Keys.OrderBy(k => k, StringComparer.Ordinal));
        return kinds;
    }

    /// <summary>
    /// Builds the network at the requested target. For anything other than <see cref="ExactTarget"/>
    /// the free nodes of the network are moved onto that target too, so Theseus objectives that
    /// read current positions fit the same geometry Inverse FDM receives.
    /// </summary>
    /// <param name="c">The parsed, validated case.</param>
    /// <param name="targetKind"><see cref="ExactTarget"/>, <see cref="DesignerTarget"/>, or a <c>target_variants</c> key.</param>
    /// <param name="snugBox">Use <c>bounds_snug</c> when the case carries it; otherwise <c>bounds</c> with a warning.</param>
    public static InverseFdmCaseNetwork Build(InverseFdmCase c, string targetKind = ExactTarget, bool snugBox = false)
    {
        ArgumentNullException.ThrowIfNull(c);
        c.Validate();

        int n = c.Nodes.Length;
        var warnings = new List<string>();

        var isFixed = new bool[n];
        foreach (int f in c.Fixed) isFixed[f] = true;
        var freeIdx = Enumerable.Range(0, n).Where(i => !isFixed[i]).ToList();
        var fixedIdx = Enumerable.Range(0, n).Where(i => isFixed[i]).ToList();

        var available = AvailableTargets(c);
        string kind = string.IsNullOrWhiteSpace(targetKind) ? ExactTarget : targetKind.Trim();
        double[]?[]? chosen = null;
        if (kind == ExactTarget)
        {
            // handled below: target = c.Target, positions = c.Nodes
        }
        else if (kind == DesignerTarget)
        {
            if (c.TargetOriginal is null)
                warnings.Add($"Case '{c.Name}' has no designer target (target_original); using the exact target.");
            else if (freeIdx.Any(i => c.TargetOriginal[i] is null))
            {
                int missing = freeIdx.Count(i => c.TargetOriginal[i] is null);
                warnings.Add($"Case '{c.Name}' designer target covers only {freeIdx.Count - missing} of {freeIdx.Count} free nodes; using the exact target.");
            }
            else chosen = c.TargetOriginal;
        }
        else if (c.TargetVariants is not null && c.TargetVariants.TryGetValue(kind, out var variant))
        {
            chosen = variant;
        }
        else
        {
            warnings.Add($"Case '{c.Name}' has no target '{kind}' (available: {string.Join(", ", available)}); using the exact target.");
        }
        if (chosen is null) kind = ExactTarget;

        var positions = new Point3d[n];
        var target = new Point3d[n];
        for (int i = 0; i < n; i++)
        {
            positions[i] = ToPoint(c.Nodes[i]);
            target[i] = ToPoint(c.Target[i]);
            if (chosen is not null && !isFixed[i])
            {
                target[i] = ToPoint(chosen[i]!);
                positions[i] = target[i];
            }
        }

        var nodes = new Node[n];
        for (int i = 0; i < n; i++)
            nodes[i] = new Node { Value = positions[i], Anchor = isFixed[i] };

        var edges = new List<Edge>(c.Edges.Length);
        var edgeInputMap = new List<(int branchIndex, int itemIndex)>(c.Edges.Length);
        var edgePaths = new List<GH_Path>(c.Edges.Length);
        for (int e = 0; e < c.Edges.Length; e++)
        {
            var start = nodes[c.Edges[e][0]];
            var end = nodes[c.Edges[e][1]];
            edges.Add(new Edge
            {
                Start = start,
                End = end,
                Value = new LineCurve(start.Value, end.Value),
            });
            start.Neighbors.Add(end);
            end.Neighbors.Add(start);
            edgeInputMap.Add((0, e));
            edgePaths.Add(new GH_Path(0));
        }

        var free = freeIdx.Select(i => nodes[i]).ToList();
        var fixedNodes = fixedIdx.Select(i => nodes[i]).ToList();
        var graph = new Graph
        {
            Tolerance = Tolerance,
            Nodes = free.Concat(fixedNodes).ToList(),
            Edges = edges,
            EdgeInputMap = edgeInputMap,
            EdgeInputPaths = edgePaths,
        };
        graph.UpdateNodeIndices();

        bool valid = fixedNodes.Count >= 2;
        if (!valid)
            warnings.Add($"Case '{c.Name}' has {fixedNodes.Count} support(s); Ariadne needs at least 2 anchors.");

        var network = new FDM_Network
        {
            Graph = graph,
            Anchors = fixedNodes.Select(nd => nd.Value).ToList(),
            ATol = Tolerance,
            ETol = Tolerance,
            Free = free,
            Fixed = fixedNodes,
            FreeNodes = free.Select(nd => nd.Index).ToList(),
            FixedNodes = fixedNodes.Select(nd => nd.Index).ToList(),
            Valid = valid,
        };

        int m = c.Edges.Length;
        var signs = new int[m];
        for (int e = 0; e < m; e++)
        {
            double s = c.Signs is not null ? c.Signs[e] : Math.Sign(c.QRef[e]);
            signs[e] = Math.Sign(s);
        }

        bool useSnug = false;
        InverseFdmCaseBounds? box = c.Bounds;
        if (snugBox)
        {
            if (c.BoundsSnug is not null)
            {
                box = c.BoundsSnug;
                useSnug = true;
            }
            else
                warnings.Add($"Case '{c.Name}' has no snug box (bounds_snug); using its {(c.Bounds is null ? "unbounded" : "loose")} bounds.");
        }

        var lower = new double?[m];
        var upper = new double?[m];
        if (box is not null)
        {
            for (int e = 0; e < m; e++)
            {
                lower[e] = box.Lo[e];
                upper[e] = box.Hi[e];
            }
        }

        return new InverseFdmCaseNetwork
        {
            Name = c.Name,
            Network = network,
            Target = freeIdx.Select(i => target[i]).ToList(),
            Loads = freeIdx.Select(i => ToVector(c.Loads[i])).ToList(),
            LoadNodes = free.Select(nd => nd.Value).ToList(),
            Signs = signs,
            Lower = lower,
            Upper = upper,
            QRef = (double[])c.QRef.Clone(),
            TargetKind = kind,
            AvailableTargets = available,
            UsedSnugBox = useSnug,
            Warnings = warnings,
        };
    }

    private static Point3d ToPoint(double[] xyz) => new(xyz[0], xyz[1], xyz[2]);
    private static Vector3d ToVector(double[] xyz) => new(xyz[0], xyz[1], xyz[2]);
}
