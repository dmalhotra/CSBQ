% Use Nystrom solve for rigid SBT to check torus drag edgewise vs analytic SBT.
%
% Codes the tensor case of local SBT term (make a function?)
% Major radius R, minor radius eps (note eps is not the aspect ratio b/a).
% Uses const vel data u(s) = (1,0,0), solves 1D IE SBT[f] = u, then integrate
% f(s) vector to get total force F. Compares against SBT analytic solve
% by Johnson-Wu '79, p.272.
%
% Barnett 1/17/22
clear; verb = 1;

mu = 1;               % fluid viscosity (keep at 1 to get drag coeff)
R = 1.0;              % circle (major) radius
eps = 1e-4;           % fiber radius << R
[Z,Zp] = ellipse_map(R,R,eye(3),zeros(3,1));   % torus (circle rad R), xy-plane
fprintf('circle (major radius %g)\tfiber radius eps=%.3g...\n',R,eps)

p = 12;                       % order
npans = 10:10:30;

for i=1:numel(npans)         % h-convergence study ..........
  npan = npans(i);
  tpan = 2*pi*(0:npan)'/npan;   % pan t-param breakpoints (first=0, last=2pi)
  pan = setup_pans(tpan,p);
  pan = map_pans(pan,Z,Zp);
  [pan sbrk] = arccoords_pans(pan);
  s = vertcat(pan.s);                    % arc coords of nodes
  N = numel(s);                          % # total nodes
  tx = horzcat(pan.tx);                  % s-hat's at nodes, 3*N
  logterm = log(8*R/eps);                % >0
  % SBT "local" tensor: (I-3shat.shat^T) + 2*(I+shat.shat^T).logterm, not const
  Lambda = (1+2*logterm)*eye(3*N);       % start w/ isotropic part of Lambda
  for j=1:N, j3=3*j+(-2:0);              % loop over node inds, dyadic part
    shat = tx(:,j);                      % this node hat{s}, unit tangent
    Lambda(j3,j3) = Lambda(j3,j3) + (-3+2*logterm)*(shat*shat');  % 3x3 diag blk
  end
  K = nyst_K_SBT(pan,sbrk);
  A = (1/(8*pi*mu)) * ( Lambda + K );    % system mat
  u = kron(ones(N,1),[1;0;0]);           % const vel U = (1,0,0)
  f = A\u;
  w = vertcat(pan.w);                    % arc-len quadr wei per nodes, col
  f = reshape(f,[3 N]);                  % each row a Cartesian
  F = sum(w'.*f,2);                      % total force vec = rowwise wei sum
  fprintf('N=%d:\tdrag force F=(%.16g, %.3g, %.3g) \tshould be x-only\n',N,F(1),F(2),F(3))
end
FJW = mu*R * 2*pi^2*(3*logterm-17/2) / ((logterm-.5)*(logterm-2)-2);
fprintf('compare Fx_Johnson-Wu %.16g\n',FJW)

if verb
  figure; plot(s,f,'.-'); axis tight; title('solution vec f(s) cmpts');
  xlabel('s (arc-length coord)'); ylabel('f (force density cmpt)');
  legend('f_x','f_y','f_z');
  %print -dpng figs/dragz_curvexy_ellipse.png
  %figure; l=eig(A); plot(l,'+'); axis equal; title('spec A');
end
