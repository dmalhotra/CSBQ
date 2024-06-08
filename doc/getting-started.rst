.. _csbq_getting-started:

Getting Started
===============

.. note::

    CSBQ requires a C++11 compliant compiler with OpenMP 4.0 support. It has been tested with GCC-9 and newer on macOS and Linux platforms.


Downloading CSBQ
----------------

To get started, download the latest version of CSBQ from the `CSBQ GitHub repository <https://github.com/dmalhotra/CSBQ>`_:

.. code-block:: bash

   git clone --recurse-submodules https://github.com/dmalhotra/CSBQ.git


Using CSBQ in Your Projects
---------------------------

CSBQ is built on top of the `SCTL library <https://sctl.readthedocs.io/>`_, which is included as a submodule.
Everything is contained within the ``sctl::`` namespace.
As a header-only library, CSBQ does not require compilation or installation.
Simply include the ``csbq.hpp`` header file in your C++ project and start using the provided classes and functions:

.. code-block:: cpp

   #include <csbq.hpp>

Ensure the compiler can locate the header file by providing the path to ``CSBQ_ROOT/include`` using the flag ``-I ${CSBQ_ROOT}/include``.
Since CSBQ relies on SCTL, the path to ``sctl.hpp`` must also be provided to the compiler.
See the `SCTL usage instructions <https://sctl.readthedocs.io/>`_ for more details.

The following compiler flags are **required** to compile code with the CSBQ library:

- ``-std=c++11``: Enables C++11 standard.
- ``-fopenmp``: Enables OpenMP support for parallelization.

Additional **required** flags for the CSBQ library:

- ``-DSCTL_DATA_PATH=$(SCTL_DATA_PATH)``: Specifies the path for quadrature tables, where ``SCTL_DATA_PATH`` defaults to ``./data``.
- ``-DSCTL_QUAD_T=__float128``: Enables quadruple precision for stored quadrature tables.


Optional Compiler Flags
------------------------

Several optional compiler flags can be used to enable additional features or optimizations:

- ``-O3``: Maximum optimization.
- ``-march=native``: Enables SIMD vectorization and generates code optimized for the local machine.
- ``-DSCTL_MEMDEBUG``: Enables memory checks specific to CSBQ. This incurs a significant performance penalty, only use it for debugging.
- ``-DSCTL_PROFILE=5``: Enables profiling with a level of 5.
- ``-DSCTL_VERBOSE``: Enables verbose output.


Optional Libraries
------------------

The CSBQ library can be optionally configured to use several additional libraries for enhanced functionality:

- **BLAS**: Enable by defining ``SCTL_HAVE_BLAS``.
- **LAPACK**: Enable by defining ``SCTL_HAVE_LAPACK``.
- **libmvec**: Enable by defining ``SCTL_HAVE_LIBMVEC``.
- **Intel SVML**: Enable by defining ``SCTL_HAVE_SVML``.
- **MPI**: Enable by defining ``SCTL_HAVE_MPI`` (see `Comm <https://sctl.readthedocs.io/en/latest/doxygen/comm.html>`_).
- `FFTW <https://www.fftw.org>`_: Enable double precision by defining ``SCTL_HAVE_FFTW``, single precision by defining ``SCTL_HAVE_FFTWF``, or long double precision by defining ``SCTL_HAVE_FFTWL`` (see `FFT <https://sctl.readthedocs.io/en/latest/doxygen/fft_wrapper.html>`_).
- `PVFMM <http://pvfmm.org>`_: Enable by defining ``SCTL_HAVE_PVFMM`` (requires MPI, see `ParticleFMM <https://sctl.readthedocs.io/en/latest/doxygen/fmm-wrapper.html>`_).

To enable support for any of these libraries, define the corresponding flag during compilation. For example, to enable MPI support, use ``-DSCTL_HAVE_MPI``.


Compiling and Running Demo Codes
--------------------------------

Demo codes for learning to use the library are provided in the ``tutorial/`` directory and precomputed quadrature tables for Laplace and Stokes kernels and some geometry files are provided in the ``data/`` directory.
In the CSBQ directory, modify the ``Makefile`` as necessary for your environment, particularly to set the C++ compiler.
Several optional flags can also be enabled as documented above.
The included ``Makefile`` can also be used as a template for new projects.
Then, compile and run the demo codes as follows:

.. code-block:: bash

    make && ./bin/demo1-geometry

For visualization (after installing ParaView):

.. code-block:: bash

    make && ./bin/demo1-geometry && paraview vis/ring.pvtu

In addition, test codes used to generate the results in the `paper <http://dx.doi.org/10.1016/j.jcp.2024.112855>`_ are provided in the ``test/`` directory, along with SLURM scripts in the ``scripts/`` directory.
For slender-body theory numerical implementations in MATLAB, see the ``SBT/`` directory.


Features and Capabilities
-------------------------

The underlying `SCTL library <https://sctl.readthedocs.io/>`_ provides several useful features, such as
containers
(`Vector <https://sctl.readthedocs.io/en/latest/tutorial/vector.html>`_,
`Matrix <https://sctl.readthedocs.io/en/latest/tutorial/matrix.html>`_),
utilities
(`profiling <https://sctl.readthedocs.io/en/latest/tutorial/profile.html>`_,
`SIMD vectorization <https://sctl.readthedocs.io/en/latest/tutorial/vec.html>`_,
`writing VTK visualizations <https://sctl.readthedocs.io/en/latest/tutorial/vtudata.html>`_, etc.),
`GMRES <https://sctl.readthedocs.io/en/latest/tutorial/gmres.html>`_ for solving linear systems,
`SDC (Spectral Deferred Correction) <https://sctl.readthedocs.io/en/latest/tutorial/sdc.html>`_ for high order time-stepping,
`Kernel functions <https://sctl.readthedocs.io/en/latest/doxygen/kernel_functions.html>`_,
`ParticleFMM <https://sctl.readthedocs.io/en/latest/tutorial/fmm.html>`_ and
`BoundaryIntegralOp <https://sctl.readthedocs.io/en/latest/tutorial/boundaryintegralop.html>`_ for building boundary integral operators.

The `CSBQ library <https://github.com/dmalhotra/CSBQ>`_ extends these features with:

- :ref:`SlenderElemList class <demo1>` for discretizing slender fibers, which can be used with SCTL's `BoundaryIntegralOp <https://sctl.readthedocs.io/en/latest/tutorial/boundaryintegralop.html>`_ class to evaluate :ref:`layer-potentials <demo2>` and solve :ref:`boundary integral equations <demo3>`.

- :ref:`RigidBodyList and Mobility class <mobility_hpp>` for representing rigid slender fibers and solving Stokes mobility problems on these geometries.

- :ref:`Utilities <utils_hpp>` such as the `CubeVolumeVis class` for writing VTK files to visualize 3D volume data, routines for constructing certain predefined geometries, and for parsing command-line options.


