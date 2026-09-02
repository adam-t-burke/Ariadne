use ariadne_lbfgsb::{Backend, Bounds, Control, Options, Solver};
use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
use std::convert::Infallible;

struct CountingAllocator;

thread_local! {
    static ENABLED: Cell<bool> = const { Cell::new(false) };
    static ALLOCATIONS: Cell<usize> = const { Cell::new(0) };
}

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        ENABLED.with(|enabled| {
            if enabled.get() {
                ALLOCATIONS.with(|count| count.set(count.get() + 1));
            }
        });
        unsafe { System.alloc(layout) }
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        unsafe { System.dealloc(ptr, layout) }
    }

    unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        ENABLED.with(|enabled| {
            if enabled.get() {
                ALLOCATIONS.with(|count| count.set(count.get() + 1));
            }
        });
        unsafe { System.realloc(ptr, layout, new_size) }
    }
}

fn extended_rosenbrock(x: &[f64], gradient: &mut [f64]) -> Result<f64, Infallible> {
    let mut value = 0.25 * (x[0] - 1.0).powi(2);
    for i in 1..x.len() {
        let residual = x[i] - x[i - 1].powi(2);
        value += residual * residual;
    }
    value *= 4.0;

    let mut residual = x[1] - x[0].powi(2);
    gradient[0] = 2.0 * (x[0] - 1.0) - 16.0 * x[0] * residual;
    for i in 1..x.len() - 1 {
        let previous = residual;
        residual = x[i + 1] - x[i].powi(2);
        gradient[i] = 8.0 * previous - 16.0 * x[i] * residual;
    }
    gradient[x.len() - 1] = 8.0 * residual;
    Ok(value)
}

/// This integration-test binary deliberately contains one test. Counting is
/// also thread-local, so allocations from Cargo's test harness and unrelated
/// test processes cannot be attributed to the solve.
#[test]
fn warm_solver_allocates_nothing_during_solve() {
    const N: usize = 200;
    let mut lower = [-100.0; N];
    let upper = [100.0; N];
    for i in (0..N).step_by(2) {
        lower[i] = 1.0;
    }
    #[cfg(feature = "faer-backend")]
    let backends = &[Backend::Deterministic, Backend::Faer];
    #[cfg(not(feature = "faer-backend"))]
    let backends = &[Backend::Deterministic];

    for &backend in backends {
        for history in [2, 10, 20] {
            let bounds = Bounds::new(&lower, &upper, N).unwrap();
            let options = Options::new()
                .with_history_size(history)
                .unwrap()
                .with_projected_gradient_tolerance(1.0e-10)
                .unwrap()
                .with_backend(backend)
                .unwrap();
            let mut solver = Solver::new(options);
            let mut x = [3.0; N];

            let warm = solver
                .minimize_with_callback(&mut x, bounds, extended_rosenbrock, |_| Control::Continue)
                .unwrap();
            assert!(warm.stats.accepted_updates > options.history_size());

            x.fill(3.0);
            ALLOCATIONS.with(|count| count.set(0));
            ENABLED.with(|enabled| enabled.set(true));
            let report = solver
                .minimize_with_callback(&mut x, bounds, extended_rosenbrock, |_| Control::Continue)
                .unwrap();
            ENABLED.with(|enabled| enabled.set(false));
            let allocations = ALLOCATIONS.with(Cell::get);

            assert!(report.stats.accepted_updates > options.history_size());
            assert!(report.stats.line_search_probes >= report.stats.iterations);
            assert_eq!(
                report.stats.evaluations,
                report.stats.line_search_probes + 1
            );
            assert_eq!(
                report.stats.active_variables + report.stats.free_variables,
                N
            );
            assert_eq!(
                allocations, 0,
                "warm same-thread solver-owned allocations for {backend:?}, m={history}"
            );
        }
    }
}
