# Convergent slender-body quadrature (CSBQ) --- research implementation

authors: **Dhairya Malhotra**, with slender-body theory comparisons by
**Alex Barnett**

This is a high-performance parallel C++ implementation of a high-order
Nystrom quadrature for the boundary integral equations arising
in 3D Laplace and Stokes Dirichlet and rigid mobility boundary-value problems
for closed loop filaments of arbitrarily small circular cross-section.
Its cost is independent of the slenderness parameter.

This repository also contains MATLAB codes implementing the classical
slender-body theory asymptotic approximation,
and solving its linear inverse problem as needed for a mobility solve.

It is research software used as part of investigations of PDE
solvers at the Center for Computational Mathematics at the Flatiron Institute.

Use at your own risk.


### Minimum requirements to compile:

C++ compiler that supports C++11 standard


### Optional libraries for better performance:

BLAS, LAPACK, FFTW, MPI (distributed memory parallelism), PVFMM

### Instructions:

```bash
git clone https://github.com/dmalhotra/CSBQ.git
cd CSBQ
git submodule init
git submodule update
# update Makefile as necessary
make
./bin/demo1-geometry
```

For visualization do (after installing PARAVIEW):

```bash
make && ./bin/demo1-geometry && paraview vis/ring.pvtu
```

For slender-body numerical implementations in MATLAB, see `SBT/` directory.


### Overview:
This is a header-only C++ library, you only need to `#include <csbq.hpp>` in your code.
It relies on the [SCTL](https://github.com/dmalhotra/SCTL) library (included as a submodule) and the path to `sctl.hpp` must be provided to the C++ compiler.
The included Makefile may be used as a template for new projects.
Everything is contained within `sctl::` namespace.
Demo codes for learning to use the library are provided in `tutorial/`.
Precomputed quadrature tables for Laplace and Stokes kernels and some geometry files are provided in `data/`.
Test codes used to generate the results in the paper are provided in `test/`, along with SLURM scripts in `scripts/`.
