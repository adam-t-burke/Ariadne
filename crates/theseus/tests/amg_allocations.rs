//! After the first `solve`, `AmgSolver::solve` must not allocate on the heap
//! (program plan §4.1 / §7.3 WS-C). This binary deliberately holds one test
//! and a process-wide counting allocator: the counter is global (not
//! thread-local) so that allocations made on rayon worker threads during
//! the solve are attributed as well.

#[path = "support/fixtures/mod.rs"]
mod fixtures;

use fixtures::grid;
use ndarray::Array2;
use std::alloc::{GlobalAlloc, Layout, System};
use std::sync::atomic::{AtomicBool, AtomicUsize, Ordering};
use theseus::amg::AmgSolver;
use theseus::backend::CpuBackend;
use theseus::linear_solver::{
    IterativeSolverOptions, LinearSystemSolver, Precision, SolveRequest, TolerancePolicy,
};
use theseus::types::*;

struct CountingAllocator;

static ENABLED: AtomicBool = AtomicBool::new(false);
static ALLOCATIONS: AtomicUsize = AtomicUsize::new(0);

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        if ENABLED.load(Ordering::Relaxed) {
            ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        }
        unsafe { System.alloc(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        if ENABLED.load(Ordering::Relaxed) {
            ALLOCATIONS.fetch_add(1, Ordering::Relaxed);
        }
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

fn count_allocations(f: impl FnOnce()) -> usize {
    ALLOCATIONS.store(0, Ordering::SeqCst);
    ENABLED.store(true, Ordering::SeqCst);
    f();
    ENABLED.store(false, Ordering::SeqCst);
    ALLOCATIONS.load(Ordering::SeqCst)
}

#[test]
fn warm_solver_allocates_nothing_during_solve() {
    // 96 × 96: above the rayon threshold on level 0, several AMG levels.
    let problem = grid::make_recoverable_grid_problem(96);
    let q = fixtures::smooth_q_star(problem.topology.num_edges);
    let mut cache = FdmCache::new(&problem).unwrap();
    theseus::fdm::solve_fdm(&mut cache, &q, &problem, &Array2::zeros((0, 3)), 0.0).unwrap();
    let rhs = cache.rhs.as_slice().unwrap().to_vec();
    let n3 = rhs.len();

    // Warm up the rayon pool (its lazy initialisation allocates) before the
    // measured region.
    rayon::broadcast(|_| {});

    for precision in [Precision::F64, Precision::F32] {
        let options = IterativeSolverOptions {
            tolerance: TolerancePolicy::Fixed(1e-10),
            coarsest_size: 200,
            precondition_precision: Some(precision),
            ..IterativeSolverOptions::default()
        };
        let mut solver: AmgSolver<CpuBackend> =
            AmgSolver::cpu(&problem.topology, &problem.bounds, &options).unwrap();
        solver.update(&q).unwrap();
        assert!(
            solver.level_sizes().len() >= 3,
            "{:?}",
            solver.level_sizes()
        );

        let mut x = vec![0.0; n3];
        let x0 = vec![0.0; n3];
        let request = || SolveRequest {
            rhs: &rhs,
            x0: Some(&x0),
            tolerance: 1e-10,
            max_iterations: 500,
            cancel: None,
        };
        let warm = solver.solve(request(), &mut x).unwrap();
        assert!(warm.converged);

        // Same system again: identical work, no allocation.
        let allocations = count_allocations(|| {
            let stats = solver.solve(request(), &mut x).unwrap();
            assert!(stats.converged);
            assert_eq!(stats.iterations, warm.iterations);
        });
        assert_eq!(allocations, 0, "warm solve allocated ({precision:?})");

        // After a numeric update the buffers are still the same size: the
        // solve must stay allocation-free too.
        let q2: Vec<f64> = q.iter().map(|v| v * 1.02).collect();
        solver.update(&q2).unwrap();
        assert_eq!(solver.counters().setups, 1);
        let allocations = count_allocations(|| {
            let stats = solver.solve(request(), &mut x).unwrap();
            assert!(stats.converged);
        });
        assert_eq!(
            allocations, 0,
            "solve after update allocated ({precision:?})"
        );
    }
}
