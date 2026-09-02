# Third-party notices

## L-BFGS-B 3.0 reference implementation

The test fixtures under `tests/reference/fixtures` were generated with the
official L-BFGS-B 3.0 Fortran distribution by Ciyou Zhu, Richard Byrd, Jorge
Nocedal, and José Luis Morales. The upstream source is downloaded and compiled
only by the opt-in fixture generator; it is not included in or linked into the
production crate.

The authoritative project page states:

> This software is freely available, but we expect that all publications
> describing work using this software, or all commercial products using it,
> quote at least one of the references given below. This software is released
> under the "New BSD License" (aka "Modified BSD License" or "3-clause
> license").

Upstream: <http://users.iems.northwestern.edu/~nocedal/lbfgsb.html>

The downloaded distribution contains `License.txt`, preserved verbatim as
`UPSTREAM_LICENSE.txt`. The file is an incomplete historical BSD license
template: its copyright line and organization/holder fields are blank. This
project does not infer or assign an upstream copyright holder. The complete
BSD-3-Clause redistribution conditions and disclaimer governing this Rust
crate are in the crate-local `LICENSE`.

Requested references:

R. H. Byrd, P. Lu, and J. Nocedal, "A Limited Memory Algorithm for Bound
Constrained Optimization," SIAM Journal on Scientific and Statistical
Computing 16(5), 1995, pp. 1190-1208.

C. Zhu, R. H. Byrd, and J. Nocedal, "Algorithm 778: L-BFGS-B: Fortran
Subroutines for Large-Scale Bound-Constrained Optimization," ACM Transactions
on Mathematical Software 23(4), 1997, pp. 550-560.

J. L. Morales and J. Nocedal, "Remark on Algorithm 778: L-BFGS-B: Fortran
Subroutines for Large-Scale Bound Constrained Optimization," ACM Transactions
on Mathematical Software 38(1), 2011, pp. 7:1-7:4.

Archive identity, objective definitions, build flags, and generated-file hashes
are recorded in `tests/reference/PROVENANCE.md`.


