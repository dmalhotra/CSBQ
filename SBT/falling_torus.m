% Compute, in the SBT approx, z-mobility of rigid symmetric torus in xy-plane.
%
% Centerline is circle radius R in xy plane, minor radius eps,
% falling at constant velocity u(s) = (0,0,U), under a constant force density
% f(s) = (0,0,f). Then total force F = 2*pi*R*f.
% This is the simplest rigid sedimentation case computed to high accuracy
% (not the SBT approx) in
% Mitchell arxiv:2102.01791, Table 2, which improves on Amarakoon '82; we
% include those numbers. For this easy case the K operator in SBT vanishes.
%
% Goal: get the Lambda_zz(s) "local" SBT term right (which is not
% really local, ie, not an extensive quantity, since its value depends on the
% global circumference in a log way),
% and compare to Dhairya's code, for various major radii R and minor radii eps.
%
% Dhairya could extract the total force F via integrating traction over a
% distant enclosing sphere, or by integrating it over the body surface.
%
% Barnett 1/5/22
clear

R = 1.0; %3.7;    % major radius of torus, ie of centerline circle, can be anything

mitchepss = [1e-1 1e-2 1e-3 1e-4 1e-5];   % Mitchell Table 2 bottom rows...
dragcoeffs = [0.7843079118776118     % F', ie force ratio to a sphere
              0.5773112010340116     % radius 1+eps = max torus radius.
              0.4410799504734741     % Note Mitchell eps is scaled for R=1 only
              0.3552543625215476
              0.2972352562536550];
ne = numel(mitchepss);

mu = 1;       % bulk viscosity (simply multiplies F, so no point in changing)
U = 1;        % const vel (0,0,U) bdry data for Stokes vel-BVP (ditto)
L = 2*pi*R;   % circumf

for e=1:ne                       % loop over eps cases
  eps = R*mitchepss(e);        % true inner rad
  fprintf('Major rad R=%g, minor rad eps=%g...\n',R,eps)
  Fsph = 6*pi*mu*U*(R+eps);   % Stokes sphere drag
  Fex = Fsph * dragcoeffs(e);     % "exact" (Mitchell claims 12 digits)
  Lambda_zz = 1 - 2*log(pi*eps/(4*L));   % Ohm general-L "local" term, scalar
  fSBT = 8*pi*mu*U / Lambda_zz;        % solve SBT eqn for f, const scalar
  FSBT = L*fSBT;                  % integrate force dens over circumf
  fprintf('\ttot F exact\t%.12g\n\ttot F SBT\t%.12g\n', Fex, FSBT)
  fprintf('\trel diff |FSBT-Fex|/Fex:\t\t%.8g\n', abs(FSBT/Fex-1))
  fprintf('\tempirical fit c.(eps/R)^2.log(R/eps):\t%.8g\n', 0.235*(eps/R)^2*log(R/eps))
end
% Observe error in SBT here *very* well fit by O((eps/R)^2 . log(R/eps))

