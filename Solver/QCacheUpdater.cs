namespace Ariadne.Solver;

using System.Collections.Generic;
using System.Linq;

internal static class QCacheUpdater
{
    public static void Apply(
        IDictionary<EdgeKey, double> cache,
        IReadOnlyList<EdgeKey> keys,
        IReadOnlyList<double> forceDensities,
        IReadOnlySet<EdgeKey> ambiguousKeys)
    {
        var seen = new HashSet<EdgeKey>();
        for (int i = 0; i < keys.Count; i++)
        {
            var key = keys[i];
            if (ambiguousKeys.Contains(key))
                continue;

            seen.Add(key);
            if (i < forceDensities.Count && double.IsFinite(forceDensities[i]))
                cache[key] = forceDensities[i];
        }

        foreach (var key in cache.Keys.Where(key => !seen.Contains(key)).ToList())
            cache.Remove(key);
    }
}
