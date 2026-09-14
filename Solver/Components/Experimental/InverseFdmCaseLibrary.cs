using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text.Json;

namespace Ariadne.Solver.Components.Experimental;

/// <summary>
/// Locates and parses external inverse-FDM benchmark cases. The JSON files under
/// <c>bench/external/cases/</c> are embedded in the assembly at build time; a folder
/// override loads <c>{name}.json</c> from disk instead (for regenerated exports).
/// </summary>
public static class InverseFdmCaseLibrary
{
    /// <summary>Logical-name prefix set in <c>Ariadne.csproj</c> for the embedded case files.</summary>
    public const string ResourcePrefix = "Ariadne.Bench.Cases.";

    private static readonly JsonSerializerOptions Options = new()
    {
        AllowTrailingCommas = true,
        ReadCommentHandling = JsonCommentHandling.Skip,
    };

    private static readonly Lazy<IReadOnlyList<string>> Bundled = new(() =>
        typeof(InverseFdmCaseLibrary).Assembly
            .GetManifestResourceNames()
            .Where(r => r.StartsWith(ResourcePrefix, StringComparison.Ordinal)
                        && r.EndsWith(".json", StringComparison.OrdinalIgnoreCase))
            .Select(r => r.Substring(ResourcePrefix.Length, r.Length - ResourcePrefix.Length - ".json".Length))
            .OrderBy(n => n, StringComparer.Ordinal)
            .ToList());

    /// <summary>Names of the cases embedded in the plugin, sorted.</summary>
    public static IReadOnlyList<string> BundledNames => Bundled.Value;

    /// <summary>Parses one case and validates its array lengths and indices.</summary>
    public static InverseFdmCase Parse(string json)
    {
        var c = JsonSerializer.Deserialize<InverseFdmCase>(json, Options)
            ?? throw new InvalidOperationException("Case JSON is empty.");
        c.Validate();
        return c;
    }

    /// <summary>Loads an embedded case by name (filename stem, e.g. <c>jaxfdm_arch_loadpath</c>).</summary>
    public static InverseFdmCase LoadBundled(string name)
    {
        string resource = ResourcePrefix + name + ".json";
        using Stream? stream = typeof(InverseFdmCaseLibrary).Assembly.GetManifestResourceStream(resource);
        if (stream is null)
            throw new FileNotFoundException(
                $"No bundled case named '{name}'. Available: {string.Join(", ", BundledNames)}.");
        using var reader = new StreamReader(stream);
        return Parse(reader.ReadToEnd());
    }

    /// <summary>Loads <c>{name}.json</c> (or <paramref name="name"/> verbatim if it already ends in .json) from a folder.</summary>
    public static InverseFdmCase LoadFromFolder(string folder, string name)
    {
        if (!Directory.Exists(folder))
            throw new DirectoryNotFoundException($"Case folder not found: {folder}");
        string file = name.EndsWith(".json", StringComparison.OrdinalIgnoreCase) ? name : name + ".json";
        string path = Path.Combine(folder, file);
        if (!File.Exists(path))
            throw new FileNotFoundException($"Case file not found: {path}");
        return Parse(File.ReadAllText(path));
    }

    /// <summary>Folder override when non-empty, otherwise the embedded case.</summary>
    public static InverseFdmCase Load(string name, string? folder)
    {
        ArgumentException.ThrowIfNullOrWhiteSpace(name);
        return string.IsNullOrWhiteSpace(folder)
            ? LoadBundled(name.Trim())
            : LoadFromFolder(folder.Trim(), name.Trim());
    }
}
