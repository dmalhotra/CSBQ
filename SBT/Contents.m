% SBT
%
% Files
%   dragz_curvexy        - Compute drag of xy-plane curve under const (0,0,U) velocity, in SBT.
%   ellipse_map          - simple 3D ellipse closed curve analytic map and its perim
%   falling_torus        - Compute, in the SBT approx, z-mobility of rigid symmetric torus in xy-plane.
%   Kzz_mat              - fill zz compnt of SBT K matrix for closed fiber loop
%   LIquad_panels        - add arc-length quadrature weights to panels given 3D nodes
%   load_geom            - read in 3D centerline panel quadrature from a Dhairya geom file
%   map_pans             - use chart to add 3D quadrature to parameter panel quadrature
%   nyst_diagdiscont_sca - Nystrom discretize diag-discontinuous kernel, 1D panels
%   pan_brkpts           - sum up quadrature weights in panels and cumsum to breakpoints
%   save_geom            - write 3D centerline panel quadrature to Dhairya .geom file format
%   setup_pans           - parameter panel quadrature struct array from 1D breakpoints
%   showcurve            - plot panel quadrature for a curve in 3D
%   startup              - some project- and Alex-specific MATLAB settings for SBT.
%   test_all             - run all tests of local functions and new utilities, no elaborate drivers
%   test_spec_toy        - Check spectrum of scalar SBT operator (zz-part) on unit circle vs analytic.
