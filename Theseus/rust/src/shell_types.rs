//! Data types for shell FEA: triangle/quad elements with 6 DOFs per node.
//! Prioritized for future Cosserat/Micropolar expansion.

use crate::sparse::SparseColMatOwned;
use crate::types::{Factorization, TheseusError, find_nz_index};
use ndarray::Array2;
use std::fmt;

#[derive(Debug, Clone)]
pub struct ShellMaterial {
    pub e: f64,
    pub nu: f64,
    pub density: f64,
    pub yield_stress: f64,
}

#[derive(Debug, Clone)]
pub struct ShellSection {
    pub offset: f64,
}

#[derive(Debug, Clone)]
pub struct ShellElementProps {
    pub material_idx: usize,
    pub section_idx: usize,
}

#[derive(Debug, Clone)]
pub struct ShellLoad {
    pub node_idx: usize,
    pub force: [f64; 3],
    pub moment: [f64; 3],
}

#[derive(Debug, Clone)]
pub struct ShellSupport {
    pub node_idx: usize,
    pub fixed_dofs: [bool; 6], // 3 translations, 3 rotations
}

#[derive(Debug, Clone)]
pub struct ShellDofMap {
    pub num_nodes: usize,
    pub num_total_dofs: usize,
    pub num_free_dofs: usize,
    pub num_fixed_dofs: usize,
    pub global_to_free: Vec<Option<usize>>,
    pub free_dofs: Vec<usize>,
    pub fixed_dofs: Vec<usize>,
}

impl ShellDofMap {
    pub fn from_supports(num_nodes: usize, supports: &[ShellSupport]) -> Self {
        let num_total = num_nodes * 6;
        let mut is_fixed = vec![false; num_total];
        for sup in supports {
            for d in 0..6 {
                if sup.fixed_dofs[d] {
                    is_fixed[sup.node_idx * 6 + d] = true;
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

#[derive(Debug, Clone)]
pub struct ShellElementToNz {
    pub entries: Vec<Vec<(usize, usize, usize)>>,
}

#[derive(Debug)]
pub struct ShellProblem {
    pub num_nodes: usize,
    pub num_elements: usize,
    pub materials: Vec<ShellMaterial>,
    pub sections: Vec<ShellSection>,
    pub element_props: Vec<ShellElementProps>,
    pub supports: Vec<ShellSupport>,
    pub loads: Vec<ShellLoad>,
    pub node_positions: Array2<f64>,
    pub node_thicknesses: Vec<f64>,
    pub elements: Vec<Vec<usize>>,
    pub dof_map: ShellDofMap,
    pub include_self_weight: bool,
    pub gravity: [f64; 3],
}

pub struct ShellCache {
    pub k_matrix: SparseColMatOwned,
    pub factorization: Option<Factorization>,
    pub element_to_nz: ShellElementToNz,
    pub displacements: Vec<f64>,
    pub rhs: Vec<f64>,
    pub reactions: Vec<f64>,
    pub dof_map: ShellDofMap,
}

impl fmt::Debug for ShellCache {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ShellCache(n_free_dof={}, n_elem={})",
            self.dof_map.num_free_dofs, self.element_to_nz.entries.len())
    }
}

impl ShellCache {
    pub fn new(problem: &ShellProblem) -> Result<Self, TheseusError> {
        let ne = problem.num_elements;
        let n_free = problem.dof_map.num_free_dofs;
        let n_total = problem.dof_map.num_total_dofs;

        // Build K's sparsity pattern from element connectivity
        let mut triplets: Vec<(u32, u32, f64)> = Vec::new();
        for e in 0..ne {
            let nodes = &problem.elements[e];
            let num_nodes = nodes.len();
            let num_dofs = num_nodes * 6;
            let mut global_dofs = Vec::with_capacity(num_dofs);
            for &ni in nodes {
                for d in 0..6 {
                    global_dofs.push(ni * 6 + d);
                }
            }
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
            let nodes = &problem.elements[e];
            let num_nodes = nodes.len();
            let num_dofs = num_nodes * 6;
            let mut global_dofs = Vec::with_capacity(num_dofs);
            for &ni in nodes {
                for d in 0..6 {
                    global_dofs.push(ni * 6 + d);
                }
            }
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

        Ok(Self {
            k_matrix,
            factorization: None,
            element_to_nz: ShellElementToNz { entries: elem_entries },
            displacements: vec![0.0; n_free],
            rhs: vec![0.0; n_free],
            reactions: vec![0.0; n_total],
            dof_map: problem.dof_map.clone(),
        })
    }
}

#[derive(Debug, Clone)]
pub struct ShellResult {
    pub displacements: Vec<f64>,
    pub reactions: Vec<f64>,
    pub principal_stresses: Vec<[f64; 2]>,
    pub membrane_stresses_global: Vec<[f64; 6]>,
    pub top_stresses_global: Vec<[f64; 6]>,
    pub bottom_stresses_global: Vec<[f64; 6]>,
    pub von_mises: Vec<f64>,
}
