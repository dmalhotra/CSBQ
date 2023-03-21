% Use Nystrom solve for rigid SBT to check torus drag edgewise vs Johnson-Wu'79.
%
% Codes the tensor case of local SBT term (make a function?)
% Major radius R, minor radius eps (note eps is not the aspect ratio b/a).
% Uses const vel data u(s) = (1,0,0), solves 1D IE SBT[f] = u, then integrate
% f(s) vector to get total force F. Compares against SBT analytic solve
% by Johnson-Wu '79, p.272. Align torus in xy-plane since we tested rot by now.
%
% Barnett 1/17/22
clear; verb = 0;

mu = 1;               % fluid viscosity (keep at 1 to get drag coeff)
R = 3.7;              % circle (major) radius
eps = 1e-2;           % fiber radius << R
[Z,Zp] = ellipse_map(R,R,eye(3),zeros(3,1));   % torus (circle rad R), xy-plane
fprintf('circle (major) radius %g\tfiber radius eps=%.3g...\n',R,eps)

p = 12;                       % order
npans = 10:10:30;

for i=1:numel(npans)         % h-convergence study ..........
  npan = npans(i);
  tpan = 2*pi*(0:npan)'/npan;   % pan t-param breakpoints (first=0, last=2pi)
  pan = setup_pans(tpan,p);
  pan = map_pans(pan,Z,Zp);
  [pan sbrk] = arccoords_pans(pan);
  N = p*npan;                            % # total nodes
  logterm = log(8*R/eps);                % >0, circle case
  Lambda = Lambda_SBT(pan, logterm);     % fill 2 SBT terms...
  K = nyst_K_SBT(pan,sbrk);
  A = (1/(8*pi*mu)) * ( Lambda + K );    % usual SBT prefac, system mat
  u = kron(ones(N,1),[1;0;0]);           % const edgewise vel U = (1,0,0)
  f = A\u;
  f = reshape(f,[3 N]);                  % each col a 3-vec
  w = vertcat(pan.w);                    % arc-len quadr wei per nodes, col
  F = f*w;                               % use quadr for total force vec
  fprintf('N=%d:\tdrag force F=(%.16g, %.3g, %.3g) \tcond(A)=%.3g\n',N,F(1),F(2),F(3),cond(A))  % cond growing when N>1/eps since spec denser @ 0
end                          % .............................

FJW = mu*R * 2*pi^2*(3*logterm-17/2) / ((logterm-.5)*(logterm-2)-2);   % cool!
fprintf('compare Fx Johnson-Wu %.16g\n',FJW)       % an exact solve of SBT

if verb
  s = vertcat(pan.s)';                    % arc-len of nodes, row
  figure; plot(s,f,'.-'); axis tight; title('solution vec f(s) cmpts');
  xlabel('s (arc-length coord)'); ylabel('f (force density cmpt)');
  legend('f_x','f_y','f_z');
  %figure; l=eig(A); plot(l,'+'); axis equal; title('spec A');
end
