.. _utils_hpp:

utils.hpp
=========

This header file provides a collection of utility functions and classes.

Geometry Generation Functions
-----------------------------

Functions to generate various predefined geometries for use in slender body simulations.

- ``GenericGeom``: Builds geometry from a given geometric function, allowing customization of panel lengths, Fourier orders, and other parameters.

- ``GenericGeomAdap``: Builds adaptively refined geometry from a given geometric function, with tolerance for adaptive refinement.

- ``GeomEllipse``: Builds a torus geometry with specified radii and thickness, customizable with panel lengths and Fourier orders.

- ``GeomTouchingTori``: Builds geometry with two nearly touching tori with a specified separation, enabling the creation of complex configurations.

- ``GeomTangle``: Generates AB tangle geometry, facilitating the simulation of tangled filament networks.

- ``GeomSphere``: Generates geometry of a sphere with specified radius, panel lengths, and Fourier orders.

..

Solver Functions
----------------

Functions to verify layer-potential quadratures using Green's identity and to solve exterior Dirichlet boundary value problems (BVPs) using different integral equation formulations.

- ``double_layer_test``: Performs a double-layer identity test to validate the implementation of a given kernel type.

- ``test_greens_identity``: Tests Green's identity for Laplace and Stokes kernels, validating their correctness.

- ``bvp_solve``: Solves an exterior Dirichlet boundary value problem (BVP) using various kernel types, providing the solution at off-surface points.

- ``bvp_solve_combined``: Solves an exterior Dirichlet BVP using a combined kernel, providing the solution at off-surface points.

..

Visualization Class
-------------------

- ``CubeVolumeVis``:  For uniformly discretizing a cube volume and generating visualization data in VTK format.

..

Command-Line Option Parsing Functions
-------------------------------------

Utility functions to parse command-line options, including starting and ending option parsing and retrieving option values.

- ``commandline_option_start``: Starts parsing command-line options, optionally displaying help text.

- ``commandline_option``: Parses a command-line option, retrieving its value or displaying an error message if not found.

- ``commandline_option_end``: Ends parsing command-line options after processing all arguments.

|

.. raw:: html

   <div style="border-top: 1px solid"></div>
   <br>

.. literalinclude:: ../../include/csbq/utils.hpp
   :language: c++
