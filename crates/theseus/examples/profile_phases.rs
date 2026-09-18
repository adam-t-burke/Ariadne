//! Per-phase wall-clock profile of one fused objective/gradient evaluation on
//! square grid networks, from a few thousand to ~100k edges.
//!
//!   cargo run --release -p theseus --example profile_phases -- [grid sizes...]
//!
//! Default grid sizes: 72 (10k edges), 160 (51k), 224 (100k). Set
//! `RAYON_NUM_THREADS` to compare threading configurations.

#[path = "../tests/support/grid.rs"]
#[allow(dead_code)]
mod grid;

use ndarray::Array2;
use std::time::Instant;
use theseus::types::*;

pub fn ms(start: Instant) -> f64 {
    start.elapsed().as_secs_f64() * 1e3
}

fn main() {
    let sizes: Vec<usize> = std::env::args()
        .skip(1)
        .map(|a| a.parse().expect("grid size"))
        .collect();
    let sizes = if sizes.is_empty() {
        vec![72, 160, 224]
    } else {
        sizes
    };
    let threads = std::env::var("RAYON_NUM_THREADS").unwrap_or_else(|_| "default".into());
    println!("rayon threads: {threads}");
    println!(
        "{:>5} {:>7} {:>9} | {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} | {:>9}",
        "grid",
        "edges",
        "nodes",
        "asm_A",
        "asm_rhs",
        "factor",
        "solve",
        "geom",
        "loss",
        "expl",
        "adjoint",
        "impl",
        "total"
    );
    for &n in &sizes {
        let problem = grid::make_grid_problem(n);
        let ne = problem.topology.num_edges;
        let nn = problem.topology.num_nodes;
        let anchors = Array2::zeros((0, 3));
        let q = vec![1.0; ne];
        let mut cache = FdmCache::new(&problem).unwrap();
        // Warm up: symbolic analysis and first numeric factorization.
        theseus::fdm::solve_fdm(&mut cache, &q, &problem, &anchors, 1e-12).unwrap();

        let reps = if ne > 60_000 { 3 } else { 5 };
        let mut acc = [0.0f64; 9];
        let mut total = 0.0;
        for _ in 0..reps {
            let t_all = Instant::now();
            cache.q.copy_from_slice(&q);
            let t = Instant::now();
            theseus::fdm::assemble_a(&mut cache);
            acc[0] += ms(t);
            theseus::fdm::update_fixed_positions(&mut cache, &problem, &anchors);
            let t = Instant::now();
            theseus::fdm::assemble_rhs(&mut cache, &problem);
            acc[1] += ms(t);
            // factor + solve are fused in factor_and_solve; time the refactor separately.
            let t = Instant::now();
            {
                let fac = cache.factorization.as_mut().unwrap();
                fac.update(&cache.a_matrix, &mut cache.factor_stack)
                    .unwrap();
            }
            acc[2] += ms(t);
            let t = Instant::now();
            {
                let fac = cache.factorization.as_ref().unwrap();
                fac.solve_into(
                    &cache.rhs,
                    &mut cache.x,
                    &mut cache.solve_workspace,
                    &mut cache.solve_stack,
                )
                .unwrap();
            }
            acc[3] += ms(t);
            for (i, &node) in problem.topology.free_node_indices.iter().enumerate() {
                for d in 0..3 {
                    cache.nf[[node, d]] = cache.x[[i, d]];
                }
            }
            let t = Instant::now();
            theseus::fdm::compute_geometry(&mut cache, &problem);
            acc[4] += ms(t);
            let t = Instant::now();
            let snap = GeometrySnapshot {
                xyz_full: &cache.nf,
                member_lengths: &cache.member_lengths,
                member_forces: &cache.member_forces,
                reactions: &cache.reactions,
            };
            let loss = theseus::objectives::total_loss(&problem.objectives, &snap);
            acc[5] += ms(t);
            std::hint::black_box(loss);
            let t = Instant::now();
            cache.grad_q.fill(0.0);
            cache.grad_nf.fill(0.0);
            theseus::gradients::accumulate_explicit_gradients(&mut cache, &problem);
            acc[6] += ms(t);
            let t = Instant::now();
            theseus::gradients::solve_adjoint(&mut cache).unwrap();
            acc[7] += ms(t);
            let t = Instant::now();
            theseus::gradients::accumulate_implicit_gradients(&mut cache, &problem);
            acc[8] += ms(t);
            total += ms(t_all);
        }
        let r = reps as f64;
        println!(
            "{:>5} {:>7} {:>9} | {:>8.2} {:>8.2} {:>8.2} {:>8.2} {:>8.2} {:>8.2} {:>8.2} {:>8.2} {:>8.2} | {:>9.2}",
            n, ne, nn, acc[0] / r, acc[1] / r, acc[2] / r, acc[3] / r, acc[4] / r, acc[5] / r, acc[6] / r, acc[7] / r, acc[8] / r, total / r
        );
    }
}
