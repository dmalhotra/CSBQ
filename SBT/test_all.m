% Run all tests of local functions and non-trivial utilities.
startup

% in rough order from low- to high-level:
interpmat_1d
load_geom
save_geom
setup_pans
pan_brkpts
map_pans
arccoords_pans
nyst_diagdiscont_sca     % checks eigenmodes of inv-sin kernel
nyst_Kzz_SBT             % checks circle Kzz kernel is inv-sin kernel
nyst_K_SBT               % checks 1) rot-equiv to Kzz, 2) Fou modes convergent
drag_torus_edgewise_SBT  % validates tensor SBT vs Johnson'79
%drag_torus_axial        % science not a self-test
