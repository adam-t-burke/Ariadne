using System;
using System.Runtime.InteropServices;

namespace Theseus.Interop;

/// <summary>
/// Raw P/Invoke declarations for the FEA functions in theseus.dll.
/// Maps 1:1 to <c>extern "C"</c> symbols in <c>rust/src/fea_ffi.rs</c>.
/// </summary>
internal static class FeaInterop
{
    const string DLL = "theseus";

    // ── Handle lifecycle ─────────────────────────────────────

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern IntPtr theseus_fea_create(
        nuint num_nodes, nuint num_elements,
        double[] node_positions, nuint[] edge_nodes,
        nuint num_materials, double[] materials,
        nuint num_sections, double[] sections,
        nuint[] element_props,
        nuint num_supports, nuint[] supports,
        nuint num_loads, double[] loads_data, double[] load_moments, nuint[] load_nodes,
        int include_self_weight, double[] gravity,
        byte beam_formulation);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void theseus_fea_free(IntPtr handle);

    // ── Objective registration ───────────────────────────────

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_fea_add_compliance(IntPtr handle, double weight);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_fea_add_max_displacement(
        IntPtr handle, double weight, nuint[] node_indices, nuint num_nodes);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_fea_add_target_displacement(
        IntPtr handle, double weight, nuint[] node_indices, nuint num_nodes, double[] targets);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_fea_add_min_weight(IntPtr handle, double weight);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_fea_add_max_stress(
        IntPtr handle, double weight, nuint[] edge_indices, nuint num_edges,
        double[] thresholds, double sharpness);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_fea_add_target_geometry(
        IntPtr handle, double weight, nuint[] node_indices, nuint num_nodes, double[] targets);

    // ── Solver options ───────────────────────────────────────

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_fea_set_solver_options(
        IntPtr handle, nuint max_iterations, double abs_tol, double rel_tol,
        double barrier_weight, double barrier_sharpness);

    // ── Variable mask ────────────────────────────────────────

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_fea_set_variable_mask(
        IntPtr handle,
        int opt_node_positions, nuint[] node_indices, nuint num_node_indices,
        int opt_areas, nuint[] area_indices, nuint num_area_indices,
        int opt_support_positions, nuint[] support_node_indices, nuint num_support_indices);

    // ── Progress callback ────────────────────────────────────

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    public delegate byte FeaNativeProgressCallback(
        nuint iteration, double loss, IntPtr xyz, nuint numNodes);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_fea_set_progress_callback(
        IntPtr handle, FeaNativeProgressCallback? callback, nuint frequency);

    // ── Forward solve ────────────────────────────────────────

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_fea_solve_forward(
        IntPtr handle,
        double[] out_displacements, double[] out_reactions,
        double[] out_axial_forces, double[] out_stresses,
        double[] out_strains, double[] out_deformed_xyz,
        double[] out_utilization,
        double[] out_internal_forces);

    // ── Optimization ─────────────────────────────────────────

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_fea_optimize(
        IntPtr handle,
        double[] out_displacements, double[] out_reactions,
        double[] out_axial_forces, double[] out_stresses,
        double[] out_strains, double[] out_deformed_xyz,
        double[] out_utilization,
        double[] out_node_positions, double[] out_areas,
        ref nuint out_iterations, ref byte out_converged);
}
