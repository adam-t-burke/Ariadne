//! Data types for solid (continuum) FEA: tet4 elements.

use crate::sparse::SparseColMatOwned;
use crate::types::{Factorization, TheseusError, find_nz_index};
use ndarray::Array2;
use std::fmt;

#[derive(Debug, Clone)]
pub struct SolidMaterial {
    pub e: f64,
    pub nu: f64,
    pub density: f64,
    pub yield_stress: f64,
}

#[derive(Debug, Clone)]
pub struct SolidElementProps {
    pub material_idx: usize,
}

#[derive(Debug, Clone)]
pub struct SolidLoad {
    pub node_idx: usize,
    pub force: [f64; 3],
}

#[derive(Debug, Clone)]
pub struct SolidSupport {
    pub node_idx: usize,
    pub fixed_dofs: [bool; 3],
}

#[derive(Debug, Clone)]
pub struct SolidDofMap {
    pub num_nodes: usize,
    pub num_total_dofs: usize,
    pub num_free_dofs: usize,
    pub num_fixed_dofs: usize,
    pub global_to_free: Vec<Option<usize>>,
    pub free_dofs: Vec<usize>,
    pub fixed_dofs: Vec<usize>,
}

impl SolidDofMap {
    pub fn from_supports(num_nodes: usize, supports: &[SolidSupport]) -> Self {
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

#[derive(Debug, Clone)]
pub struct SolidElementToNz {
    pub entries: Vec<Vec<(usize, usize, usize)>>,
}

#[derive(Debug)]
pub struct SolidProblem {
    pub num_nodes: usize,
    pub num_elements: usize,
    pub materials: Vec<SolidMaterial>,
    pub element_props: Vec<SolidElementProps>,
    pub supports: Vec<SolidSupport>,
    pub loads: Vec<SolidLoad>,
    pub node_positions: Array2<f64>,
    pub elements: Vec<[usize; 4]>,
    pub dof_map: SolidDofMap,
    pub include_self_weight: bool,
    pub gravity: [f64; 3],
}

pub struct SolidCache {
    pub k_matrix: SparseColMatOwned,
    pub factorization: Option<Factorization>,
    pub element_to_nz: SolidElementToNz,
    pub dof_map: SolidDofMap,
    pub displacements: Vec<f64>,
    pub rhs: Vec<f64>,
    pub stresses: Vec<[f64; 6]>,
    pub strains: Vec<[f64; 6]>,
    pub von_mises: Vec<f64>,
    pub reactions: Vec<f64>,
}

impl fmt::Debug for SolidCache {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "SolidCache(n_free_dof={}, n_elem={})",
            self.dof_map.num_free_dofs, self.stresses.len())
    }
}

impl SolidCache {
    pub fn new(problem: &SolidProblem) -> Result<Self, TheseusError> {
        let ne = problem.num_elements;
        let n_free = problem.dof_map.num_free_dofs;
        let n_total = problem.dof_map.num_total_dofs;

        let mut triplets: Vec<(u32, u32, f64)> = Vec::new();
        for e in 0..ne {
            let nodes = &problem.elements[e];
            let mut global_dofs = [0usize; 12];
            for (i, &ni) in nodes.iter().enumerate() {
                for d in 0..3 {
                    global_dofs[i * 3 + d] = ni * 3 + d;
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
        for i in 0..n_free {
            triplets.push((i as u32, i as u32, 0.0));
        }

        let k_matrix = SparseColMatOwned::from_triplets(n_free, n_free, &triplets)
            .map_err(|e| TheseusError::Shape(e))?;

        let mut elem_entries: Vec<Vec<(usize, usize, usize)>> = vec![Vec::new(); ne];
        for e in 0..ne {
            let nodes = &problem.elements[e];
            let mut global_dofs = [0usize; 12];
            for (i, &ni) in nodes.iter().enumerate() {
                for d in 0..3 {
                    global_dofs[i * 3 + d] = ni * 3 + d;
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
            element_to_nz: SolidElementToNz { entries: elem_entries },
            dof_map: problem.dof_map.clone(),
            displacements: vec![0.0; n_free],
            rhs: vec![0.0; n_free],
            stresses: vec![[0.0; 6]; ne],
            strains: vec![[0.0; 6]; ne],
            von_mises: vec![0.0; ne],
            reactions: vec![0.0; n_total],
        })
    }
}

#[derive(Debug, Clone)]
pub struct SolidResult {
    pub displacements: Vec<f64>,
    pub deformed_positions: Array2<f64>,
    pub reactions: Vec<f64>,
    pub stresses: Vec<[f64; 6]>,
    pub strains: Vec<[f64; 6]>,
    pub von_mises: Vec<f64>,
}
