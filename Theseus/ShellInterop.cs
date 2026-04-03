using System;
using System.Runtime.InteropServices;

namespace Theseus.Interop;

/// <summary>
/// Raw P/Invoke declarations for the shell FEA functions in theseus.dll.
/// Maps 1:1 to <c>extern "C"</c> symbols in <c>rust/src/shell_ffi.rs</c>.
/// </summary>
internal static class ShellInterop
{
    const string DLL = "theseus";

    // ── Handle lifecycle ─────────────────────────────────────

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern IntPtr theseus_shell_create(
        nuint numNodes, nuint numElements,
        double[] nodePositions, double[] nodeThicknesses, int[] elements, int[] numNodesPerElement,
        nuint numMaterials, double[] materials,
        nuint numSections, double[] sections,
        int[] elementProps,
        nuint numSupports, int[] supports,
        nuint numLoads, double[] loadForces, double[] loadMoments, int[] loadNodes,
        byte includeSelfWeight, double[] gravity);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void theseus_shell_destroy(IntPtr handle);

    // ── Forward solve ────────────────────────────────────────

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_shell_solve_forward(
        IntPtr handle,
        double[] outDisplacements, double[] outReactions, double[] outPrincipalStresses,
        double[] outTopStresses, double[] outBottomStresses, double[] outVonMises);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_shell_spr_recover(
        nuint numNodes, nuint numElements,
        double[] nodePositions, int[] elements, int[] numNodesPerElement,
        double[] membraneStresses, double[] topStresses, double[] bottomStresses,
        double[] outNodalMembrane, double[] outNodalTop, double[] outNodalBottom,
        double[] outElementErrors,
        double[] outPrincipalValuesTop, double[] outPrincipalValuesBot,
        double[] outPrincipalVectorsTop, double[] outPrincipalVectorsBot,
        double[] outVonMisesTop, double[] outVonMisesBot);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_shell_last_error(byte[] buf, nuint bufLen);
}
