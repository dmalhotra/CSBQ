.. _demo2:

Demo 2: Boundary Integral Operator
==================================

The code (``tutorial/demo2-bio.cpp``) demonstrates the basic usage of the ``BoundaryIntegralOp`` class.
The ``BoundaryIntegralOp`` class (in `SCTL <https://sctl.readthedocs.io>`_) is designed to construct and evaluate boundary integral operators, such as the Stokes single-layer and double-layer potential operators.
For more advanced usage and additional features, please refer to the class API in `boundary_integral.hpp <https://sctl.readthedocs.io/en/latest/doxygen/boundary_integral.html>`_.

Building a Boundary Integral Operator
-------------------------------------

To build a boundary integral operator, follow these steps:

1. Initialize the Boundary Integral Operator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Initialize the boundary integral operator with the desired kernel and communication context.
In this example, we use the Stokes double-layer kernel (``Stokes3D_DxU``).

.. code-block:: cpp

   const Stokes3D_DxU ker;
   BoundaryIntegralOp<double,Stokes3D_DxU> BIOp(ker, false, comm);

2. Set Quadrature Accuracy
~~~~~~~~~~~~~~~~~~~~~~~~~~

Set the quadrature accuracy for numerical integration.

.. code-block:: cpp

   BIOp.SetAccuracy(1e-10);

3. Add Geometry to the Operator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Add the geometry of the surface to the boundary integral operator. This includes specifying the discretization of the surface.
In this example the boundary is given by an instance of :ref:`SlenderElemList class <slender_element_hpp>` with geometry data loaded from a file.

.. code-block:: cpp

   SlenderElemList<double> elem_lst;
   elem_lst.Read<double>("data/loop.geom", comm);
   BIOp.AddElemList(elem_lst);

4. Set Evaluation Points
~~~~~~~~~~~~~~~~~~~~~~~~

The target points can be specified as follows.
If not set or ``Xt`` is empty, then the default target points are the surface discretization nodes.

.. code-block:: cpp

   BIOp.SetTargetCoord(Xt);

Evaluating Potentials
----------------------

Once the boundary integral operator is constructed, you can evaluate potentials at on- and off-surface target points.

1. Compute the Potential
~~~~~~~~~~~~~~~~~~~~~~~~

Compute the potential using the boundary integral operator and a density function defined on the surface.

.. code-block:: cpp

   Vector<double> sigma(Ninput);
   sigma = 1; // Set the density function sigma at each surface discretization node
   Vector<double> U;
   BIOp.ComputePotential(U, sigma); // Compute the potential U

2. Visualize the Results
~~~~~~~~~~~~~~~~~~~~~~~~~

Visualize the geometry and the computed potential.

.. code-block:: cpp

   elem_lst.WriteVTK("vis/U", U, comm); // Write visualization data to VTK file


Compiling and Running the Code
------------------------------

- **Without MPI**: Navigate to the project root directory and run:

   .. code-block:: bash
   
      make bin/demo2-bio && ./bin/demo2-bio

- **Without MPI**: In the project root directory, edit ``Makefile`` to set ``CXX=mpicxx`` and
  uncomment ``CXXFLAGS += -DSCTL_HAVE_MPI``. Then, run:

   .. code-block:: bash
   
      make bin/demo2-bio && mpirun -n 2 --map-by slot:pe=1 ./bin/demo2-bio

This will evaluate the Stokes double-layer potential on a loop geometry and write the VTK visualization to the file ``vis/U.pvtu`` which can be opened in ParaView.


Complete Example Code
---------------------

.. raw:: html

   <div style="border-top: 1px solid"></div>
   <br>

.. literalinclude:: ../../tutorial/demo2-bio.cpp
   :language: c++
  
