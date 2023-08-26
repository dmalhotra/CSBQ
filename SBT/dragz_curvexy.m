% Solve drag of general xy-plane curve under const (0,0,U) velocity, in SBT.
%
% Needs solve (0,0,f(s)) function in 1D IE SBT[f] = u, just zz-cmpnts, then
% integrate f to get total force F.
% This enables direct comparison with Stokes Dirichlet (vel) BVP.
% Radius is eps, const along curve.
% t is the analytic description parameter. s will be arc-len coord param.
%
% Barnett 1/14/22
clear; verb = 0;

mu = 1;               % fluid viscosity
eps = 1e-2;           % radius

a = 2; b = 0.5;
[Z,Zp] = ellipse_map(a,b,eye(3),zeros(3,1));       % ellipse, xy-plane
fprintf('ellipse (semiaxes %g,%g)\tradius eps=%.3g...\n',a,b,eps)

%a = 0.5;             % inverted ellipse
%zmap = @(t) exp(1i*t) ./ (1+a*exp(2i*t));   % C-plane
%Z = @(t) [real(zmap(t)); imag(zmap(t)); 0*t]; Zp=[];  % xy-plane
%a = 0.4;             % kite, a is narrowness
%Z = @(t) [a*cos(t)+(1-a)*cos(2*t); sin(t); 0*t];
%Zp = @(t) [-a*sin(t)-2*(1-a)*sin(2*t); cos(t); 0*t];

p = 12;                       % order
npans = 10:10:40;

for i=1:numel(npans)         % h-convergence study ..........
  npan = npans(i);
  tpan = 2*pi*(0:npan)'/npan;   % pan t-param breakpoints (first=0, last=2pi)
  pan = setup_pans(tpan,p);
  pan = map_pans(pan,Z,Zp);
  [pan sbrk] = arccoords_pans(pan);
  L = sbrk(end);
  s = vertcat(pan.s);                    % arc coords of nodes
  N = numel(s);                          % # total nodes
  Lambda_zz = 1 - 2*log(pi*eps/(4*L));   % "local" term, scalar, const
  Kzz = nyst_Kzz_SBT(pan,sbrk);
  A = (1/(8*pi*mu)) * ( Lambda_zz*eye(N) + Kzz );    % system mat
  u = ones(N,1);                                     % const vel U = (0,0,1)
  f = A\u;
  w = vertcat(pan.w);                    % arc-len quadr wei for nodes
  F = sum(w.*f);                         % total force
  fprintf('N=%d:\tL=%.16g\tdrag force F=%.16g\n',N,L,F)
end                          % ...............................

if verb
  figure; subplot(2,1,1); showcurve(pan,gca); title('curve in xy-plane');
  subplot(2,1,2); plot(s,f,'.-'); axis tight; title('solution f(s)');
  xlabel('s (arc-length coord)'); ylabel('f (force density)');
  print -dpng figs/dragz_curvexy_ellipse.png
  save_geom('output/ellipse_a2_b.5_p12_npan40_r1e-3.geom',pan,1e-3,10);
  figure; l=eig(A); plot(l,'+'); axis equal; title('spec A');
end
