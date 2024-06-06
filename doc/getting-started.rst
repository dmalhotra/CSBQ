.. _csbq_getting-started:

Getting Started with CSBQ
=========================

.. note::

    CSBQ requires a C++11 compliant compiler with OpenMP 4.0 support.

Downloading CSBQ
----------------

To get started, download the latest version of CSBQ from the `CSBQ GitHub repository <https://github.com/dmalhotra/CSBQ>`_:

.. code-block:: bash

   git clone --recurse-submodules https://github.com/dmalhotra/CSBQ.git

Using CSBQ
----------

CSBQ is built on top of the `SCTL library <https://sctl.readthedocs.io/>`_, which is included as a submodule. Everything is contained within the ``sctl::`` namespace. As a header-only library, CSBQ does not require compilation or installation. Simply include the ``csbq.hpp`` header file in your C++ project and start using the provided classes and functions:

.. code-block:: cpp

   #include <csbq.hpp>

Ensure the compiler can locate the header file by providing the path to ``CSBQ_ROOT/include`` using the flag ``-I ${CSBQ_ROOT}/include``. Since CSBQ relies on SCTL, the path to ``sctl.hpp`` must also be provided to the compiler. See the `SCTL usage instructions <https://sctl.readthedocs.io/>`_ for more details.

Optional Dependencies
---------------------

CSBQ can utilize the following libraries for better performance:

- BLAS
- LAPACK
- FFTW
- MPI (for distributed memory parallelism)
- PVFMM

Refer to the instructions on how to `enable these libraries in SCTL <https://sctl.readthedocs.io/>`_.

Compiling and Running Demo Codes
--------------------------------

In the CSBQ directory, modify the ``Makefile`` as necessary for your environment, particularly to set the C++ compiler. Several optional flags can also be enabled. Then, compile and run the demo codes as follows:

.. code-block:: bash

    make
    ./bin/demo1-geometry

For visualization (after installing ParaView):

.. code-block:: bash

    make && ./bin/demo1-geometry && paraview vis/ring.pvtu

Demo codes for learning to use the library are provided in the ``tutorial/`` directory.
Precomputed quadrature tables for Laplace and Stokes kernels and some geometry files are provided in the ``data/`` directory.
Test codes used to generate the results in the paper are provided in the ``test/`` directory, along with SLURM scripts in the ``scripts/`` directory.

Additional Resources
--------------------

- The included ``Makefile`` can be used as a template for new projects.
- For slender-body numerical implementations in MATLAB, see the ``SBT/`` directory.

Tutorials
---------

- :ref:`SlenderElemList <tutorial-slenderelemlist>`: Discretization of slender fibers and computing boundary integrals.
- :ref:`BoundaryIntegralOp <tutorial-boundaryintegralop>`: Generic class for instantiating layer-potential operators.


