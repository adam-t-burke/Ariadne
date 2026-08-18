namespace Ariadne.Solver;

using System;
using System.Collections.Generic;
using System.Linq;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;

internal sealed record QTreeMappingResult(
    IReadOnlyList<double>? Values,
    string? Error,
    string? Warning)
{
    public bool Success => Error is null && Values is not null;
}

internal static class QTreeMapper
{
    public static QTreeMappingResult Map(
        GH_Structure<GH_Number> qTree,
        IReadOnlyList<GH_Path> edgePaths)
    {
        ArgumentNullException.ThrowIfNull(qTree);
        ArgumentNullException.ThrowIfNull(edgePaths);

        int edgeCount = edgePaths.Count;
        if (edgeCount == 0)
            return Failure("The network has no edges, so force densities cannot be assigned.");

        var qBranches = new List<QBranch>(qTree.PathCount);
        int totalQCount = 0;
        for (int branchIndex = 0; branchIndex < qTree.PathCount; branchIndex++)
        {
            GH_Path path = qTree.Paths[branchIndex];
            var branch = qTree.get_Branch(path);
            var values = new List<double>(branch.Count);
            foreach (var item in branch)
            {
                if (item is GH_Number number)
                    values.Add(number.Value);
            }

            totalQCount += values.Count;
            qBranches.Add(new QBranch(path.ToString(), values));
        }

        if (totalQCount == 0)
            return Failure("Initial force densities cannot be empty.");

        if (totalQCount == 1)
        {
            double value = qBranches.SelectMany(branch => branch.Values).Single();
            return Success(Enumerable.Repeat(value, edgeCount).ToArray());
        }

        var edgeIndicesByPath = new Dictionary<string, List<int>>(StringComparer.Ordinal);
        for (int edgeIndex = 0; edgeIndex < edgePaths.Count; edgeIndex++)
        {
            string path = edgePaths[edgeIndex].ToString();
            if (!edgeIndicesByPath.TryGetValue(path, out var indices))
            {
                indices = [];
                edgeIndicesByPath.Add(path, indices);
            }
            indices.Add(edgeIndex);
        }

        if (qBranches.Count == 1 && totalQCount == edgeCount)
        {
            string? warning = edgeIndicesByPath.Count > 1
                ? $"q has {edgeCount} values matching the total edge count; values were applied 1:1 " +
                  "in flattened edge order, not per edge-tree branch. Graft/match q to the edge " +
                  "tree to assign per branch."
                : null;
            return Success(qBranches[0].Values.ToArray(), warning);
        }

        var qByPath = new Dictionary<string, IReadOnlyList<double>>(StringComparer.Ordinal);
        foreach (var qBranch in qBranches)
        {
            if (!qByPath.TryAdd(qBranch.Path, qBranch.Values))
                return Failure($"The q tree contains duplicate path {qBranch.Path}.");
        }

        foreach (var (path, edgeIndices) in edgeIndicesByPath)
        {
            if (!qByPath.ContainsKey(path))
            {
                return Failure(
                    $"No q values were supplied for edge branch {path}, which contains " +
                    $"{edgeIndices.Count} edge(s). Supply q values for every edge branch.");
            }
        }

        foreach (var (path, values) in qByPath)
        {
            if (values.Count > 0 && !edgeIndicesByPath.ContainsKey(path))
                return Failure($"q branch {path} does not match any edge branch.");
        }

        var result = new double[edgeCount];
        foreach (var (path, edgeIndices) in edgeIndicesByPath)
        {
            IReadOnlyList<double> values = qByPath[path];
            if (values.Count == 1)
            {
                foreach (int edgeIndex in edgeIndices)
                    result[edgeIndex] = values[0];
                continue;
            }

            if (values.Count != edgeIndices.Count)
            {
                return Failure(
                    $"q at {path} has {values.Count} value(s), but that branch has " +
                    $"{edgeIndices.Count} edge(s). Supply 1 value to broadcast or " +
                    $"{edgeIndices.Count} values to assign per edge.");
            }

            for (int itemIndex = 0; itemIndex < edgeIndices.Count; itemIndex++)
                result[edgeIndices[itemIndex]] = values[itemIndex];
        }

        return Success(result);
    }

    private static QTreeMappingResult Success(
        IReadOnlyList<double> values,
        string? warning = null) =>
        new(values, null, warning);

    private static QTreeMappingResult Failure(string error) =>
        new(null, error, null);

    private sealed record QBranch(string Path, IReadOnlyList<double> Values);
}
