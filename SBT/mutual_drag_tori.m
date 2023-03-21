% Testing breakdown of SBT (Nystrom solve) for mutual drag between two near tori
%
% Major radius R, minor radius eps (note eps is not the aspect ratio b/a).
% Torus 1 is stationary; torus 2 has u(s) = U2 = (1,0,0).
% We solve 1D (2x2 block for 2 objects) IE SBT[f] = u, then integrate
% f(s) vector over each object to get total forces F1 and F2.
% Sign convention: Fj is the net force externally applied to body j.
%
% Barnett 3/17/23
clear; verb = 1;

mu = 1;               % fluid viscosity (keep at 1 to get drag coeff)
R = 1.0;              % circle (major) radius
eps = 1e-3;           % fiber radius << R
delta = 0.1; %10*eps;       % surface separation (centerline sep = delta+2eps)
U1 = [0;0;0]; U2 = [0;0;1];     % imposed rigid vels of each torus
c1 = zeros(3,1); c2 = [delta+2*(R+eps);0;0];    % centers of tori for delta sep
[Z1,Zp1] = ellipse_map(R,R,eye(3),c1);     % torus 1, xy-plane
[Z2,Zp2] = ellipse_map(R,R,eye(3),c2);
fprintf('circle (major) radii R=%g, fiber radii eps=%g, surface min separation delta=%g...\n',R,eps,delta)

p = 12;                       % order
npans = 10:10:50;

for i=1:numel(npans)         % h-convergence study ..........
  npan = npans(i);
  N = p*npan;                            % # nodes per torus
  tpan = 2*pi*(0:npan)'/npan;            % pan t-param brkpts (first=0 last=2pi)
  pan = setup_pans(tpan,p);
  pan1 = map_pans(pan,Z1,Zp1);
  pan2 = map_pans(pan,Z2,Zp2);
  [pan1 sbrk1] = arccoords_pans(pan1);
  [pan2 sbrk2] = arccoords_pans(pan2);
  logterm = log(8*R/eps);                % >0, circle case, two tori same
  Lambda = Lambda_SBT([pan1 pan2], logterm);  % local SBT term (each pan alone)
  K11 = nyst_K_SBT(pan1,sbrk1);
  K22 = nyst_K_SBT(pan2,sbrk2);
  K21 = nyst_K_SBT_nonself(pan2,pan1);
  K12 = nyst_K_SBT_nonself(pan1,pan2);
  A = (1/(8*pi*mu)) * ( Lambda + [K11,K12;K21,K22] );  % SBT prefac, system mat
  u = [kron(ones(N,1),U1); kron(ones(N,1),U2)];    % imposed vel on all nodes
  f = A\u;
  w = vertcat(pan1.w);                   % 1-torus node arc-len quadr weis
  f = reshape(f,[3 2*N]);                % each col is a 3-vec
  F1 = f(:,1:N)*w;                       % total force vec on torus 1
  F2 = f(:,N+1:end)*w;                   % total force vec on torus 2
  fprintf('N=%d:\tmutual F1=(%.12g, %.12g, %.12g)   \tself F2=(%.12g, %.12g, %.12g)\n',N,F1(1),F1(2),F1(3),F2(1),F2(2),F2(3));
  if N<250, fprintf('\t\tcond(A) = %.3g\n',cond(A)); end
end                          % .............................

if verb  % show geom
  figure(1); clf; showcurve(pan,gca); title('curve panels'); view(3);axis vis3d;
  warning('off','');   % kill a dumb arrow axes limits change warning
  arrow(c1,c1+U1); text(c1(1)+U1(1),c1(2)+U1(2),c1(3)+U1(3),'U_1');
  arrow(c2,c2+U2); text(c2(1)+U2(1),c2(2)+U2(2),c2(3)+U2(3),'U_2');
  
  s = vertcat(pan.s)';                    % arc-len of nodes, row
  figure(2); clf; plot(s,f,'.-'); axis tight; title('solution vec f(s) cmpts');
  xlabel('s (arc-length coord)'); ylabel('f (force density cmpt)');
  legend('f_x','f_y','f_z');
  %figure; l=eig(A); plot(l,'+'); axis equal; title('spec A');
end
