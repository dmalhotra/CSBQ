% Solve SBT 1D IE for drag on general rigid smooth panelized 3D closed curve.
%
% Does convergence study for a parameterized (analytic) curve.
% This is the general 3D version of dragz_curvexy.m, which is a special case.
%
% Barnett 1/18/22
clear; verb = 1;

mu = 1;               % fluid viscosity
eps = 1e-3;   % radius (const) % 0.05054952 gets converged kappa > 1e8

a = 2; b = 0.5;       % semiaxes
Z = ellipse_map(a,b,eye(3),zeros(3,1));       % ellipse, xy-plane
fprintf('ellipse (semiaxes %g,%g)\tradius eps=%.3g...\n',a,b,eps)
  
U = [3;2;1]; U = U/norm(U);    % RHS data: const drag vel
fprintf('   data U=(\t%.12g\t%.12g\t%.12g)\n',U(1),U(2),U(3))

p = 12;                       % order
npans = 10:10:40;
for i=1:numel(npans)         % h-convergence study ..........
  npan = npans(i);
  N = p*npan;                   % total # nodes
  tpan = 2*pi*(0:npan)'/npan;   % pan t-param breakpoints (first=0, last=2pi)
  pan = setup_pans(tpan,p);
  pan = map_pans(pan,Z);
  [pan sbrk] = arccoords_pans(pan);
  perim = sbrk(end);
  logterm = log((4*perim)/(pi*eps));     % classical SBT, Laurel verified
  Lambda = Lambda_SBT(pan, logterm);
  K = nyst_K_SBT(pan,sbrk);
  A = (1/(8*pi*mu)) * ( Lambda + K );    % system mat
  u = kron(ones(N,1),U);                 % RHS: const vel vec U at all nodes
  f = A\u;
  f = reshape(f,[3 N]);                  % each col a 3-vec
  w = vertcat(pan.w);                    % arc-len quadr wei for nodes, col
  F = f*w;                               % use quadr for total force vec
  fprintf('N=%d:\tF=(\t%.12g\t%.12g\t%.12g)\tk(A)=%.3g\n',N,F(1),F(2),F(3),cond(A))
end                           % ..............................

if verb
  figure(1); clf; subplot(2,1,1); showcurve(pan,gca); title('curve panels');
  hold on; plot3([0 U(1)],[0 U(2)],[0 U(3)],'k.-'); text(U(1),U(2),U(3),'U');
  s = vertcat(pan.s);                    % arc coords of nodes
  subplot(2,1,2); plot(s,f,'.-'); axis tight; title(sprintf('solution f(s), for \\epsilon=%g',eps));
  legend('f_x','f_y','f_z');
  xlabel('s (arc-length coord)'); ylabel('f (force density cmpnt)');
  % print -dpng figs/drag_curve_ellipse.png
  %  save_geom('output/ellipse_a2_b.5_p12_npan40_r1e-3.geom',pan,1e-3,10);
end
if verb>1, figure; l=eig(A); plot(l,'+'); axis equal; title('spec A'); end
