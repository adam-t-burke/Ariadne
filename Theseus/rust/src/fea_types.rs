//! FEA-specific data structures for linear elastic truss analysis.
//!
//! Phase 1: 3D bar/truss elements (3 DOF per node, axial only).
//! Phase 2 will extend to beam elements (6 DOF per node).

use crate::sparse::SparseColMatOwned;
use crate::types::{Factorization, TheseusError, find_nz_index};
use ndarray::Array2;
use std::fmt;

// ─────────────────────────────────────────────────────────────
//  Material
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FeaMaterial {
    pub e: f64,
    pub density: f64,
    pub yield_stress: f64,
}

impl FeaMaterial {
    pub fn default_steel() -> Self {
        Self {
            e: 210e9,
            density: 7850.0,
            yield_stress: 250e6,
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Section
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FeaSection {
    pub area: f64,
}

impl Default for FeaSection {
    fn default() -> Self {
        Self { area: 0.01 }
    }
}

// ─────────────────────────────────────────────────────────────
//  Per-element property assignment
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FeaElementProps {
    pub material_idx: usize,
    pub section_idx: usize,
}

// ─────────────────────────────────────────────────────────────
//  Support (boundary condition)
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FeaSupport {
    pub node_idx: usize,
    pub fixed_dofs: [bool; 3],
}

impl FeaSupport {
    pub fn pinned(node_idx: usize) -> Self {
        Self { node_idx, fixed_dofs: [true, true, true] }
    }

    pub fn roller_x(node_idx: usize) -> Self {
        Self { node_idx, fixed_dofs: [true, false, false] }
    }

    pub fn roller_y(node_idx: usize) -> Self {
        Self { node_idx, fixed_dofs: [false, true, false] }
    }

    pub fn roller_z(node_idx: usize) -> Self {
        Self { node_idx, fixed_dofs: [false, false, true] }
    }
}

// ─────────────────────────────────────────────────────────────
//  Load
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FeaLoad {
    pub node_idx: usize,
    pub force: [f64; 3],
}

// ─────────────────────────────────────────────────────────────
//  DOF mapping
// ─────────────────────────────────────────────────────────────

/// Maps between global DOFs and free/constrained DOF partitions.
#[derive(Debug, Clone)]
pub struct DofMap {
    pub num_nodes: usize,
    pub num_total_dofs: usize,
    pub num_free_dofs: usize,
    pub num_fixed_dofs: usize,
    /// For each global DOF, its index in the free partition (None if constrained).
    pub global_to_free: Vec<Option<usize>>,
    /// Ordered list of free global DOF indices.
    pub free_dofs: Vec<usize>,
    /// Ordered list of constrained global DOF indices.
    pub fixed_dofs: Vec<usize>,
}

impl DofMap {
    pub fn from_supports(num_nodes: usize, supports: &[FeaSupport]) -> Self {
        let num_total = num_nodes * 3;
        let mut is_fixed = vec![false; num_total];
        for sup in supports {
            for d in 0..3 {
                if sup.fixed_dofs[d] {
                    is_fixed[sup.node_idx * 3 + d] = true;
                }
            }
        }

        let mut free_dofs = Vec::new();
        let mut fixed_dofs = Vec::new();
        let mut global_to_free = vec![None; num_total];

        for i in 0..num_total {
            if is_fixed[i] {
                fixed_dofs.push(i);
            } else {
                global_to_free[i] = Some(free_dofs.len());
                free_dofs.push(i);
            }
        }

        Self {
            num_nodes,
            num_total_dofs: num_total,
            num_free_dofs: free_dofs.len(),
            num_fixed_dofs: fixed_dofs.len(),
            global_to_free,
            free_dofs,
            fixed_dofs,
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Sparsity mapping: element contributions → K.data[] indices
// ─────────────────────────────────────────────────────────────

/// Pre-computed contribution of element `e` to the CSC nzval array of K.
/// For each element: list of (nz_index_in_K_data, local_row, local_col)
/// where local_row/col are in [0..6) for a bar element.
#[derive(Debug, Clone)]
pub struct ElementToNz {
    pub entries: Vec<Vec<(usize, usize, usize)>>,
}

// ─────────────────────────────────────────────────────────────
//  Problem definition
// ─────────────────────────────────────────────────────────────

#[derive(Debug)]
pub struct FeaProblem {
    pub num_nodes: usize,
    pub num_elements: usize,
    pub materials: Vec<FeaMaterial>,
    pub sections: Vec<FeaSection>,
    pub element_props: Vec<FeaElementProps>,
    pub supports: Vec<FeaSupport>,
    pub loads: Vec<FeaLoad>,
    pub node_positions: Array2<f64>,
    pub edge_nodes: Vec<(usize, usize)>,
    pub dof_map: DofMap,
    pub include_self_weight: bool,
    pub gravity: [f64; 3],
    pub objectives: Vec<Box<dyn FeaObjectiveTrait>>,
    pub solver: FeaSolverOptions,
}

// ─────────────────────────────────────────────────────────────
//  Solver options
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FeaSolverOptions {
    pub absolute_tolerance: f64,
    pub relative_tolerance: f64,
    pub max_iterations: usize,
    pub report_frequency: usize,
    pub barrier_weight: f64,
    pub barrier_sharpness: f64,
}

impl Default for FeaSolverOptions {
    fn default() -> Self {
        Self {
            absolute_tolerance: 1e-6,
            relative_tolerance: 1e-6,
            max_iterations: 500,
            report_frequency: 1,
            barrier_weight: 10.0,
            barrier_sharpness: 10.0,
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  Variable mask (which design variables to optimize)
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FeaVariableMask {
    pub node_positions: bool,
    pub node_indices: Vec<usize>,
    pub cross_section_areas: bool,
    pub area_element_indices: Vec<usize>,
    pub support_positions: bool,
    pub support_node_indices: Vec<usize>,
}

impl Default for FeaVariableMask {
    fn default() -> Self {
        Self {
            node_positions: true,
            node_indices: Vec::new(),
            cross_section_areas: false,
            area_element_indices: Vec::new(),
            support_positions: false,
            support_node_indices: Vec::new(),
        }
    }
}

// ─────────────────────────────────────────────────────────────
//  FEA Cache (mutable workspace)
// ─────────────────────────────────────────────────────────────

/// All mutable workspace for the forward solve, adjoint, and gradient
/// accumulation. Built once from a FeaProblem, reused across iterations.
pub struct FeaCache {
    // ── Sparse system ──────────────────────────────────────
    pub k_matrix: SparseColMatOwned,
    pub factorization: Option<Factorization>,
    pub element_to_nz: ElementToNz,

    // ── DOF mapping ────────────────────────────────────────
    pub dof_map: DofMap,

    // ── Element geometry (updated when node positions change) ──
    pub elem_lengths: Vec<f64>,
    pub elem_cos: Vec<[f64; 3]>,

    // ── Primal buffers ─────────────────────────────────────
    pub displacements: Vec<f64>,
    pub rhs: Vec<f64>,
    pub lambda: Vec<f64>,

    // ── Gradient buffers ───────────────────────────────────
    pub grad_u: Vec<f64>,

    // ── Post-processed results ─────────────────────────────
    pub axial_forces: Vec<f64>,
    pub stresses: Vec<f64>,
    pub strains: Vec<f64>,
    pub reactions: Vec<f64>,
    pub utilization: Vec<f64>,
}

impl fmt::Debug for FeaCache {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FeaCache(n_free_dof={}, n_elem={})",
            self.dof_map.num_free_dofs, self.elem_lengths.len())
    }
}

impl FeaCache {
    /// Build a fully pre-allocated cache from a FeaProblem.
    pub fn new(problem: &FeaProblem) -> Result<Self, TheseusError> {
        let ne = problem.num_elements;
        let n_free = problem.dof_map.num_free_dofs;
        let n_total = problem.dof_map.num_total_dofs;

        // Build K's sparsity pattern from element connectivity
        let mut triplets: Vec<(u32, u32, f64)> = Vec::new();
        for e in 0..ne {
            let (ni, nj) = problem.edge_nodes[e];
            let global_dofs = [ni * 3, ni * 3 + 1, ni * 3 + 2,
                               nj * 3, nj * 3 + 1, nj * 3 + 2];
            for &gi in &global_dofs {
                if let Some(fi) = problem.dof_map.global_to_free[gi] {
                    for &gj in &global_dofs {
                        if let Some(fj) = problem.dof_map.global_to_free[gj] {
                            triplets.push((fi as u32, fj as u32, 0.0));
                        }
                    }
                }
            }
        }
        // Add tiny diagonal to ensure all free DOFs appear in the pattern
        for i in 0..n_free {
            triplets.push((i as u32, i as u32, 0.0));
        }

        let k_matrix = SparseColMatOwned::from_triplets(n_free, n_free, &triplets)
            .map_err(|e| TheseusError::Shape(e))?;

        // Build element_to_nz mapping
        let mut elem_entries: Vec<Vec<(usize, usize, usize)>> = vec![Vec::new(); ne];
        for e in 0..ne {
            let (ni, nj) = problem.edge_nodes[e];
            let global_dofs = [ni * 3, ni * 3 + 1, ni * 3 + 2,
                               nj * 3, nj * 3 + 1, nj * 3 + 2];
            for (li, &gi) in global_dofs.iter().enumerate() {
                if let Some(fi) = problem.dof_map.global_to_free[gi] {
                    for (lj, &gj) in global_dofs.iter().enumerate() {
                        if let Some(fj) = problem.dof_map.global_to_free[gj] {
                            let nz_idx = find_nz_index(
                                &k_matrix.col_ptrs,
                                &k_matrix.row_indices,
                                fi, fj,
                            ).ok_or(TheseusError::SparsityMismatch {
                                edge: e, row: fi, col: fj,
                            })?;
                            elem_entries[e].push((nz_idx, li, lj));
                        }
                    }
                }
            }
        }

        // Compute initial element geometry
        let mut elem_lengths = vec![0.0; ne];
        let mut elem_cos = vec![[0.0; 3]; ne];
        for e in 0..ne {
            let (ni, nj) = problem.edge_nodes[e];
            let dx = problem.node_positions[[nj, 0]] - problem.node_positions[[ni, 0]];
            let dy = problem.node_positions[[nj, 1]] - problem.node_positions[[ni, 1]];
            let dz = problem.node_positions[[nj, 2]] - problem.node_positions[[ni, 2]];
            let l = (dx * dx + dy * dy + dz * dz).max(0.0).sqrt();
            elem_lengths[e] = l;
            if l > f64::EPSILON {
                elem_cos[e] = [dx / l, dy / l, dz / l];
            }
        }

        Ok(Self {
            k_matrix,
            factorization: None,
            element_to_nz: ElementToNz { entries: elem_entries },
            dof_map: problem.dof_map.clone(),
            elem_lengths,
            elem_cos,
            displacements: vec![0.0; n_free],
            rhs: vec![0.0; n_free],
            lambda: vec![0.0; n_free],
            grad_u: vec![0.0; n_free],
            axial_forces: vec![0.0; ne],
            stresses: vec![0.0; ne],
            strains: vec![0.0; ne],
            reactions: vec![0.0; n_total],
            utilization: vec![0.0; ne],
        })
    }
}

// ─────────────────────────────────────────────────────────────
//  FEA Result
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FeaResult {
    pub displacements: Vec<f64>,
    pub deformed_positions: Array2<f64>,
    pub reactions: Vec<f64>,
    pub axial_forces: Vec<f64>,
    pub stresses: Vec<f64>,
    pub strains: Vec<f64>,
    pub utilization: Vec<f64>,
}

// ─────────────────────────────────────────────────────────────
//  FEA Solver Result (from optimization)
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FeaSolverResult {
    pub fea_result: FeaResult,
    pub node_positions: Array2<f64>,
    pub areas: Vec<f64>,
    pub loss_trace: Vec<f64>,
    pub iterations: usize,
    pub converged: bool,
    pub termination_reason: String,
}

// ─────────────────────────────────────────────────────────────
//  FEA Geometry Snapshot (read-only view after forward solve)
// ─────────────────────────────────────────────────────────────

pub struct FeaSnapshot<'a> {
    pub displacements: &'a [f64],
    pub deformed_xyz: &'a Array2<f64>,
    pub axial_forces: &'a [f64],
    pub stresses: &'a [f64],
    pub strains: &'a [f64],
    pub reactions: &'a [f64],
    pub utilization: &'a [f64],
    pub elem_lengths: &'a [f64],
    pub node_positions: &'a Array2<f64>,
    pub areas: &'a [f64],
    pub densities: &'a [f64],
    pub section_indices: &'a [usize],
    pub edge_nodes: &'a [(usize, usize)],
    pub elem_cos: &'a [[f64; 3]],
    pub materials_e: &'a [f64],
}

// ─────────────────────────────────────────────────────────────
//  FEA Objective trait
// ─────────────────────────────────────────────────────────────

pub trait FeaObjectiveTrait: fmt::Debug + Send + Sync {
    fn loss(&self, snap: &FeaSnapshot) -> f64;
    fn accumulate_grad_u(&self, grad_u: &mut [f64], snap: &FeaSnapshot);
    fn weight(&self) -> f64;
    /// Whether this is a compliance objective (needs special dJ/du = 2f handling).
    fn is_compliance(&self) -> bool { false }
    /// Whether this objective has explicit dJ/d(area) contributions.
    fn has_area_gradient(&self) -> bool { false }
    fn accumulate_grad_area(&self, _grad_a: &mut [f64], _snap: &FeaSnapshot) {}
    /// Whether this objective has explicit dJ/d(node_pos) contributions.
    fn has_position_gradient(&self) -> bool { false }
    fn accumulate_grad_pos(&self, _grad_pos: &mut [f64], _snap: &FeaSnapshot) {}
}

// ─────────────────────────────────────────────────────────────
//  Optimization state
// ─────────────────────────────────────────────────────────────

#[derive(Debug, Clone)]
pub struct FeaOptimizationState {
    pub node_positions: Array2<f64>,
    pub areas: Vec<f64>,
    pub support_positions: Vec<f64>,
    pub loss_trace: Vec<f64>,
    pub iterations: usize,
}

impl FeaOptimizationState {
    pub fn new(node_positions: Array2<f64>, areas: Vec<f64>) -> Self {
        Self {
            node_positions,
            areas,
            support_positions: Vec::new(),
            loss_trace: Vec::new(),
            iterations: 0,
        }
    }
}
