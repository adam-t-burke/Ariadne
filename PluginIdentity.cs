namespace Ariadne;

/// <summary>
/// Shared Grasshopper / Yak identity for the development channel.
/// Stable public releases keep the historical "Ariadne" package name and assembly GUID.
/// Component constructors use the same ribbon label as <see cref="GhCategory"/>.
/// </summary>
public static class PluginIdentity
{
    /// <summary>Grasshopper ribbon tab for every component in this build.</summary>
    public const string GhCategory = "Ariadne Dev";

    /// <summary>Yak Package Manager package name (must not collide with stable <c>Ariadne</c>).</summary>
    public const string YakPackageName = "ariadne-dev";
}
