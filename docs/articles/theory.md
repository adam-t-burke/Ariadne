# Force density method and optimization

## Reference

**Schek, H.-J.** "The force density method for form finding and computation of general networks." *Computer Methods in Applied Mechanics and Engineering*, Vol. 3, No. 1, January 1974, pp. 115–134.

This is the classic reference for the force density method (FDM). The formulation below is a short summary; see the paper for full derivation and applications (including the Munich Olympic structures).

---

## Mathematical introduction to the force density method

### Equilibrium of a cable network

Consider a network of **branches** (cables) and **nodes**. Each branch has a force *F* and a length *L*. At each **free node**, equilibrium under applied loads requires that the vector sum of branch forces equals the external load. For **fixed nodes**, positions are given; for **free nodes**, we solve for positions.

### Force density

The **force density** of a branch is the ratio

$$q = \frac{F}{L}.$$

For a branch connecting nodes *i* and *j*, the force vector on node *i* from that branch is proportional to the direction *(x_j − x_i)* and the force density *q*. Summing over all branches and equating to the load at each free node leads to a **linear system** in the node coordinates.

### Matrix formulation

- **C** — Branch–node incidence matrix (each row is a branch; entries +1/−1 for the two incident nodes, 0 otherwise). Splitting into free and fixed columns gives **C_f** (free) and **C_x** (fixed).
- **Q** — Diagonal matrix of force densities (one per branch).
- **x** — Vector of free-node coordinates (e.g. x, y, z stacked).
- **f** — Vector of loads on free nodes (and terms from fixed-node positions).

Equilibrium for the free nodes is

$$(C_f^T Q C_f)\, x = f.$$

So for a given topology (C), fixed positions, force densities **q** (in Q), and loads **f**, the free-node positions **x** are found by solving this **linear system**. The matrix is symmetric and, under usual assumptions, positive definite, so the solution is unique and efficient to compute.

### What the force densities represent

- **q** defines the “stiffness” of each branch in the equilibrium: higher *q* means more force per unit length.
- For a given geometry and loads, **q** and branch lengths determine the member forces *F = q L*.
- In **form finding**, we often fix topology and loads and choose **q** (or geometry) to get a desired shape or distribution of forces.

---

## Optimization over parameters

In Ariadne/Theseus, the main **design variables** are the force densities **q**. The solver can operate in two modes:

### Forward solve

- **Given:** topology (nodes, branches, fixed nodes), force densities **q**, and loads **f**.
- **Compute:** equilibrium geometry (free-node positions **x**) by solving *(C_f^T Q C_f) x = f*.
- No optimization: a single linear solve.

### Inverse optimization (form finding)

- **Given:** topology, loads, bounds on **q** (e.g. *q_min*, *q_max*), and **objectives** to minimize.
- **Compute:** force densities **q** that minimize the objectives (e.g. target lengths, target forces, reaction targets, or performance metrics) subject to equilibrium and bounds.
- The solver repeatedly solves the equilibrium system for candidate **q** and updates **q** using gradient-based optimization. Objectives can involve lengths *L*, forces *F = q L*, reactions at supports, or other derived quantities.

So “optimizing over parameters” here means: **find the force densities (and thus the equilibrium shape and internal forces) that best meet your design goals**, such as matching target member lengths, limiting forces, or steering reactions, while keeping the structure in equilibrium.

In the plugin, you supply **initial q** and (for optimization) an **Opt Config** with objectives and bounds; Theseus Solve performs either the forward solve or the iterative optimization and streams the result (and optionally intermediate results) to the outputs.
