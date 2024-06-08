.. _csbq_documentation:

Convergent slender-body quadrature (CSBQ)
===========================================

.. image:: https://badgen.net/github/tag/dmalhotra/CSBQ
   :alt: Stable Version
   :target: https://github.com/dmalhotra/CSBQ/tags

.. image:: https://img.shields.io/github/v/release/dmalhotra/CSBQ?color=%233D9970
   :alt: Latest Release
   :target: https://github.com/dmalhotra/CSBQ/releases

.. image:: https://zenodo.org/badge/363025186.svg
   :alt: DOI
   :target: https://zenodo.org/doi/10.5281/zenodo.10456743

Principal author **Dhairya Malhotra**, with additional code by **Alex Barnett**.
This work was done at the Center for Computational Mathematics at the Flatiron Institute, NY, NY.

`CSBQ <https://github.com/dmalhotra/CSBQ>`_ is a high-performance parallel C++ implementation of a high-order
adaptive Nyström quadrature for the boundary integral equations arising
in 3D Laplace and Stokes Dirichlet and rigid mobility boundary-value problems
for closed loop filaments of arbitrarily small circular cross-section.
Its quadrature setup cost is independent of the slenderness parameter, being around 20000 unknowns/sec per core, at 6-digit accuracy, away from close-to-touching regions.
Close-to-touching geometries may be handled to close to machine accuracy using adaptivity.
Open-ended fibers with rounded ends are possible and will be added soon.

This repository also contains MATLAB codes implementing the classical
slender-body theory asymptotic approximation,
and solving its linear inverse problem as needed for a mobility solve.

It is research software; use at your own risk. The following figures show some of the capabilities of the code (see the publication below for details).

.. figure:: ../pics/tangle-stokes-streamlines_sm.png
   :class: with-border
   :align: center
   :width: 400

   Stokes flow solution around rigid slender fiber with aspect ratio :math:`10^3`, max error :math:`10^{-10}`.

.. figure:: ../pics/close-to-touching-streamlines_sm.png
   :class: with-border
   :align: center
   :width: 400
   
   Stokes flow solution near close-to-touching rings, max error :math:`10^{-11}`.

.. figure:: ../pics/sed512-117_sm.png
   :class: with-border
   :align: center
   :width: 400

   Sedimentation of 512 rings each of aspect ratio 20, timestepped to 7-digit accuracy on 160 cores.


Citing this work
-----------------

If you find this code useful in your research, please cite our publication:

- Dhairya Malhotra and Alex Barnett, "Efficient Convergent Boundary Integral Methods for Slender Bodies," *Journal of Computational Physics*, vol. 503, p. 112855, Apr. 2024. DOI: [10.1016/j.jcp.2024.112855](http://dx.doi.org/10.1016/j.jcp.2024.112855)



.. toctree::
   :maxdepth: 1
   :hidden:

   Introduction <self>
   getting-started
   tutorial/index

.. toctree::
   :caption: API Reference
   :maxdepth: 1
   :hidden:

   doxygen/index

