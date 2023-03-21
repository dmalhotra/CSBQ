Folder for slender body theory experiments and comparison to rigid-body BIE

All MATLAB code. Do `startup` in MATLAB before running.

Summary of low-level functions:
```matlab
%   arccoords_pans          - node and breakpt arc-length coords given 3D nodes, std speeds
%   ellipse_map             - simple 3D ellipse closed curve analytic map and its perim
%   Lambda_SBT              - fill Nystrom discretization of "local" Lambda tensor in SBT
%   LIquad_panels           - add arc-length quadrature weights, etc, given only 3D nodes
%   load_geom               - read in 3D centerline panel quadrature from a Dhairya geom file
%   map_pans                - use chart to get 3D line-integral quadrature from parameter quadr
%   nodearccoords_pans      - node arc-length coords s_j using only 3D nodes
%   nyst_diagdiscont_sca    - Nystrom discretize diag-discontinuous kernel, 1D panels
%   nyst_K_SBT              - discretize tensor SBT K operator on panel-quad closed fiber
%   nyst_K_SBT_nonself      - discretize SBT interaction between two distinct closed fibers
%   nyst_Kzz_SBT            - discretize zz-cmpnt of SBT K operator, panel-quad xy-plane fiber
%   pan_brkpts              - sum up quadrature weights in panels and cumsum to breakpoints
%   save_geom               - write 3D centerline panel quadrature to Dhairya .geom file format
%   setup_pans              - parameter panel quadrature struct array from 1D breakpoints
%   showcurve               - plot panel quadrature for a curve in 3D
%   startup                 - some project- and Alex-specific MATLAB settings for SBT.
%   test_all                - run all tests of local functions and new utilities, no elaborate drivers
```

High-level scripts/drivers and experiments:
```matlab
%   test_spec_toy           - Check spectrum of scalar SBT operator (zz-part) on unit circle vs analytic.
%   drag_torus_axial        - Compute, in the SBT approx, z-mobility of rigid symmetric torus in xy-plane.
%   drag_torus_edgewise_SBT - Use Nystrom solve for rigid SBT to check torus edgewise mobility vs Johnson-Wu'79.
%   dragz_curvexy           - Solve drag (mobility) of general xy-plane curve under const (0,0,U) velocity, in SBT.
%   drag_curve              - Solve SBT 1D IE for drag on general rigid smooth panelized 3D closed curve.
%   mutual_drag_tori        - Solve SBT for mutual mobility drag between two close tori.

```

Alex Barnett, starting late Dec 2021
