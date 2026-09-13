# compas_cem export summary

- tool: compas_cem 0.8.6 (git 0964b764783a), paper data cem_ad_cad 2b2fec3ceb3c
- `paper` = case reproduced from the CAD 2023 paper's numerical validation notebooks
- `tension/compression` counts edges by the sign of the reference force; `aux` = auxiliary trail edges (force ~ 0 after optimisation)
- `fdm_residual`: max |C^T Q C x - p| on free nodes, and relative error of the FDM re-solve of the free nodes
- `opt`: nlopt algorithm, wall time, number of objective evaluations, loss at start -> end, nlopt status

| case | paper | n_nodes | n_edges | n_fixed | tension/compression | cem_edges | params | goals | ff_time_s | opt | fdm_residual |
|---|---|---|---|---|---|---|---|---|---|---|---|
| cem_braced_tower_2d | no | 6 | 9 | 2 | 3/6 | 4 trail / 5 dev (3 indirect) | - | - | 0.0009 | - | 1.3e-15 abs / solve err 1.3e-15 rel |
| cem_bridge_2d | no | 8 | 9 | 2 | 8/1 | 6 trail / 3 dev (1 indirect) | 9 (6 TrailEdgeParameter, 3 DeviationEdgeParameter) | 2 PointGoal | 0.0004 | nlopt LD_SLSQP, 0.1s, 32 evals, loss 668 -> 2.64e-07, NLOPT_EPSVAL_REACHED | 8.9e-16 abs / solve err 5.4e-16 rel |
| cem_curved_bridge_10_3d | yes | 40 | 58 | 14 | 19/39 (10 aux) | 26 trail / 32 dev (16 indirect) | 36 (32 DeviationEdgeParameter, 4 TrailEdgeParameter) | 10 TrailEdgeForceGoal, 4 LineGoal | 0.0025 | nlopt LD_SLSQP, 1.7s, 49 evals, loss 10.4 -> 9.63e-07, NLOPT_EPSVAL_REACHED | 1.3e-14 abs / solve err 2.6e-11 rel |
| cem_curved_bridge_22_3d | yes | 88 | 130 | 26 | 45/85 (22 aux) | 62 trail / 68 dev (40 indirect) | 72 (68 DeviationEdgeParameter, 4 TrailEdgeParameter) | 22 TrailEdgeForceGoal, 4 LineGoal | 0.0054 | nlopt LD_SLSQP, 4.7s, 62 evals, loss 55.4 -> 1.42e-05, NLOPT_FTOL_REACHED | 3.2e-14 abs / solve err 3.2e-11 rel |
| cem_tensegrity_wheel_2d | no | 16 | 24 | 3 | 16/8 | 0 trail / 24 dev (0 indirect) | 24 (24 DeviationEdgeParameter) | 16 TrailEdgeForceGoal | 0.0013 | nlopt LD_LBFGS, 0.0s, 3 evals, loss 0.193 -> 4.74e-31, NLOPT_EPSVAL_REACHED | 9.2e-16 abs / solve err 2.0e-15 rel |
| cem_tensegrity_wheel_64_paper | yes | 64 | 96 | 3 | 64/32 | 0 trail / 96 dev (0 indirect) | 96 (96 DeviationEdgeParameter) | 64 TrailEdgeForceGoal | 0.0047 | nlopt LD_LBFGS, 0.1s, 3 evals, loss 52.1 -> 1.11e-28, NLOPT_EPSVAL_REACHED | 4.4e-15 abs / solve err 6.1e-15 rel |
| cem_tree_2d | no | 4 | 4 | 2 | 1/3 | 1 trail / 3 dev (0 indirect) | 3 (3 DeviationEdgeParameter) | 4 TrailEdgeForceGoal | 0.0004 | nlopt LD_SLSQP, 0.0s, 9 evals, loss 1.55 -> 9.86e-32, NLOPT_EPSVAL_REACHED | 4.4e-16 abs / solve err 2.2e-16 rel |
| cem_tree_canopy_3d | yes | 92 | 144 | 46 | 48/96 (46 aux) | 46 trail / 98 dev (0 indirect) | 186 (98 DeviationEdgeParameter, 44 OriginNodeZParameter, 44 OriginNodeYParameter) | 44 TrailEdgeForceGoal | 0.0039 | nlopt LD_SLSQP, 8.2s, 153 evals, loss 60.1 -> 9.23e-07, NLOPT_EPSVAL_REACHED | 1.4e-14 abs / solve err 2.9e-13 rel |

## CEM parameterisation vs. FDM

CEM parameters are signed trail-edge lengths (or plane offsets), deviation-edge force magnitudes, and optionally origin-node coordinates / node loads; FDM has one force density per edge.

| case | n_edges (FDM q) | CEM design params in example | trail lengths | deviation forces | origin coords | converged | stop criterion hit | tolerance |
|---|---|---|---|---|---|---|---|---|
| cem_braced_tower_2d | 9 | (form-finding only: 4 trail lengths + 5 deviation forces prescribed) | 4 of 4 | 5 of 5 | 0 | - | - | eta 1e-6 |
| cem_bridge_2d | 9 | 9 | 6 of 6 | 3 of 3 | 0 | yes | NLOPT_EPSVAL_REACHED | stopval 1e-06, ftol_abs 1e-08, max 100 evals |
| cem_curved_bridge_10_3d | 58 | 36 | 4 of 26 | 32 of 32 | 0 | yes | NLOPT_EPSVAL_REACHED | stopval 1e-06, ftol_abs 1e-08, max 100 evals |
| cem_curved_bridge_22_3d | 130 | 72 | 4 of 62 | 68 of 68 | 0 | no (loss 1.4e-05 > eps) | NLOPT_FTOL_REACHED | stopval 1e-06, ftol_abs 1e-08, max 100 evals |
| cem_tensegrity_wheel_2d | 24 | 24 | 0 of 0 | 24 of 24 | 0 | yes | NLOPT_EPSVAL_REACHED | stopval 1e-06, ftol_abs 1e-08, max 100 evals |
| cem_tensegrity_wheel_64_paper | 96 | 96 | 0 of 0 | 96 of 96 | 0 | yes | NLOPT_EPSVAL_REACHED | stopval 1e-06, ftol_abs 1e-08, max 1000 evals |
| cem_tree_2d | 4 | 3 | 0 of 1 | 3 of 3 | 0 | yes | NLOPT_EPSVAL_REACHED | stopval 1e-06, ftol_abs 1e-08, max 100 evals |
| cem_tree_canopy_3d | 144 | 186 | 0 of 46 | 98 of 98 | 88 | yes | NLOPT_EPSVAL_REACHED | stopval 1e-06, ftol_abs 1e-08, max 500 evals |

## Dropped / not reproduced

- spiral_staircase (CAD 2023 Sec. 5): only available as a binary Grasshopper definition (cem_ad_cad/case_study/staircase.gh, examples/ghpython/spiral_staircase.ghx); no Python source or geometry export
- examples/ghpython/{bridge_3d,dome,tensegrity_jessen}.ghx: Grasshopper-only definitions, geometry internalised in Rhino data
- historical examples/07_arching_tower_3d.py, 06_surface_structure.py: depend on compas_singular (CoarseQuadMesh densification), not on PyPI and COMPAS 1 only
- examples/01_quick_start.py: compression-only three-bar chain (no tension members)

## Notes

- Sign convention: compas_cem stores forces positive in tension and negative in compression for both trail and deviation edges (`Diagram.edge_force`); `q_ref = force / length` needs no sign flip. At support nodes the FDM residual `C^T Q C x - p` equals compas_cem's `reaction_force` vector.
- The reference geometry is re-solved with `static_equilibrium(eta<=1e-10, tmax=10000)`; with the examples' own `eta=1e-6` the indirect deviation edges leave nodal residuals of ~1e-7..1e-9.
- `Optimizer.solve` leaves the topology at nlopt's *last evaluated* parameters; the export re-applies `x_opt` before computing the reference equilibrium.
- Auxiliary trails (paper's extension) are exported as ordinary edges whose far node is fixed; after optimisation their forces are only as small as the loss threshold allows (`sum f^2 < 1e-6`).
