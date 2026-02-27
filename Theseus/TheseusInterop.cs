using System;
using System.Runtime.InteropServices;

namespace Theseus.Interop;

/// <summary>
/// Raw P/Invoke declarations for theseus.dll.
///
// / Every function here maps 1:1 to an <c>extern "C"</c> symbol
/// in <c>rust/src/ffi.rs</c>.  Prefer the managed <see cref="TheseusSolver"/>
/// wrapper for production use.
/// </summary>
internal static class TheseusInterop
{
    const string DLL = "theseus";

    // ── Error reporting ──────────────────────────────────────

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_last_error(byte[] buf, nuint buf_len);

    // ── Handle lifecycle ─────────────────────────────────────

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern IntPtr theseus_create(
        nuint num_edges, nuint num_nodes, nuint num_free,
        nuint[] coo_rows, nuint[] coo_cols, double[] coo_vals, nuint coo_nnz,
        nuint[] free_node_indices, nuint[] fixed_node_indices, nuint num_fixed,
        double[] loads, double[] fixed_positions,
        double[] q_init, double[] lower_bounds, double[] upper_bounds);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern void theseus_free(IntPtr handle);

    // ── Objective registration ───────────────────────────────

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_add_target_xyz(
        IntPtr handle, double weight,
        nuint[] node_indices, nuint num_nodes,
        double[] target_xyz);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_add_target_xy(
        IntPtr handle, double weight,
        nuint[] node_indices, nuint num_nodes,
        double[] target_xy);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_add_target_plane(
        IntPtr handle, double weight,
        nuint[] node_indices, nuint num_nodes,
        double[] target_xyz,
        double[] origin, double[] x_axis, double[] y_axis);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_add_planar_constraint_along_direction(
        IntPtr handle, double weight,
        nuint[] node_indices, nuint num_nodes,
        double[] origin, double[] x_axis, double[] y_axis, double[] direction);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_add_target_length(
        IntPtr handle, double weight,
        nuint[] edge_indices, nuint num_edges,
        double[] targets);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_add_length_variation(
        IntPtr handle, double weight,
        nuint[] edge_indices, nuint num_edges,
        double sharpness);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_add_force_variation(
        IntPtr handle, double weight,
        nuint[] edge_indices, nuint num_edges,
        double sharpness);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_add_sum_force_length(
        IntPtr handle, double weight,
        nuint[] edge_indices, nuint num_edges);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_add_min_length(
        IntPtr handle, double weight,
        nuint[] edge_indices, nuint num_edges,
        double[] thresholds, double sharpness);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_add_max_length(
        IntPtr handle, double weight,
        nuint[] edge_indices, nuint num_edges,
        double[] thresholds, double sharpness);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_add_min_force(
        IntPtr handle, double weight,
        nuint[] edge_indices, nuint num_edges,
        double[] thresholds, double sharpness);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_add_max_force(
        IntPtr handle, double weight,
        nuint[] edge_indices, nuint num_edges,
        double[] thresholds, double sharpness);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_add_rigid_set_compare(
        IntPtr handle, double weight,
        nuint[] node_indices, nuint num_nodes,
        double[] target_xyz);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_add_reaction_direction(
        IntPtr handle, double weight,
        nuint[] anchor_indices, nuint num_anchors,
        double[] target_dirs);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_add_reaction_direction_magnitude(
        IntPtr handle, double weight,
        nuint[] anchor_indices, nuint num_anchors,
        double[] target_dirs, double[] target_mags);

    // ── Solver options ───────────────────────────────────────

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_set_solver_options(
        IntPtr handle,
        nuint max_iterations, double abs_tol, double rel_tol,
        double barrier_weight, double barrier_sharpness);

    // ── Progress callback ────────────────────────────────────

    [UnmanagedFunctionPointer(CallingConvention.Cdecl)]
    public delegate byte NativeProgressCallback(
        nuint iteration, double loss, IntPtr xyz, nuint numNodes,
        IntPtr q, nuint numEdges);

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_set_progress_callback(
        IntPtr handle,
        NativeProgressCallback? callback,
        nuint frequency);

    // ── Optimisation ─────────────────────────────────────────

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_optimize(
        IntPtr handle,
        double[] out_xyz, double[] out_lengths, double[] out_forces,
        double[] out_q, double[] out_reactions,
        ref nuint out_iterations, ref byte out_converged);

    // ── Forward solve ────────────────────────────────────────

    [DllImport(DLL, CallingConvention = CallingConvention.Cdecl)]
    public static extern int theseus_solve_forward(
        IntPtr handle,
        double[] out_xyz, double[] out_lengths, double[] out_forces,
        double[] out_q, double[] out_reactions);
}
