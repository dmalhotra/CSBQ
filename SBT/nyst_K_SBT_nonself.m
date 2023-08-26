function K = nyst_K_SBT_nonself(tpan,span,eps)
% NYST_K_SBT_NONSELF  plain quadr tensor SBT K op from distinct closed fibers
%
% K = nyst_K_SBT_nonself(pant,pans,eps) fills a matrix, a high-order plain
%  Nystrom discretization of Johnson's '80 improvement to the Keller-Rubinow
%  classical SBT interaction of two distinct smooth closed fibers in R3.
%  We require arc-length parameters of the nodes and their arc-length weights.
%  Order p is assumed same for each panel.
%
% Inputs:
%  tpan, span - target, source struct arrays (one struct per panel) with fields:
%        x = nodes (3*p array)
%        w = weights for arc-length quadrature (length-p col vec)
%  eps - source fiber radius (const).
%  o - optional struct with optional fields:
% Output:
%  K - 3M*3N Nystrom matrix (where N, M = total nodes on source, target fiber,
%      respectively) that maps
%      values of a L-periodic (L being the source fiber perimeter) vector
%      function f(s) on all source nodes to values of
%
%    (Kf)(s) = ([S + (eps^2/2)D]f)(s) :=
%
%      int_0^L [(I + Rhat.Rhat^T)/R + (eps^2/2)(I - 3Rhat.Rhat^T)/R^3] f(s') ds'
%
%      ie a mixture of Stokeslet plus doublet (not DLP!).
%      Here ds' is arc-length element, I
%      is the 3x3 identity tensor, and hat indicates unit vector.
%      Rhat := R/|R| where R is the displacement x(s)-x(s').
%      f(s) represents force density, and Kf(s) scaled velocity.
%      Storage order is "interleaved", ie coords {xyz} fast (innermost),
%      nodes slow, ie K is a M*N block matrix with each entry a 3*3 matrix.
%
% See also: NYST_K_SBT which is the SBT self-interaction case
% 
% When called with no arguments, does self-test.
  
% Barnett 3/20/23
if nargin==0, test_nyst_K_SBT_nonself; return; end
if nargin<3, o=[]; end

nspan = numel(span);
ntpan = numel(tpan);
p = numel(span(1).w);   % order (number of nodes per pan, assumed same for each pan in src and targ)
N = nspan*p;              % total src nodes
M = ntpan*p;              % total targ nodes
inds = reshape(reshape(1:3*p, [p 3])',[1 3*p]); % pan ind list w/ xyz fast

K = nan(3*M,3*N);      % fill K mat
for ip=1:ntpan          % target pans (block rows)
  t=tpan(ip);
  ii3 = (1:3*p)+(ip-1)*3*p;   % this pan output (row) indices    
  for jp=1:nspan        % src pans
    s=span(jp);
    jj3 = (1:3*p)+(jp-1)*3*p;   % this src pan input (col) indices
    dx = t.x(1,:)' - s.x(1,:);    % go Cartesian, each cmpnt p-by-p
    dy = t.x(2,:)' - s.x(2,:);
    dz = t.x(3,:)' - s.x(3,:);
    invRdist = 1./sqrt(dx.^2+dy.^2+dz.^2);      % 1/R matrix, p*p
    dx = dx.*invRdist; dy = dy.*invRdist; dz = dz.*invRdist; % (dx,dy,dz) now Rhat
    IRR = [1+dx.^2, dx.*dy, dx.*dz; dy.*dx, 1+dy.^2, dy.*dz; dz.*dx, dz.*dy, 1+dz.^2];  % I + Rhat.Rhat^T but wrong ord: node inds fast, xyz inds slow
    Sker = IRR(inds,inds) .* kron(invRdist,ones(3));  % reord, S ker vals
    IRR = [1-3*dx.^2, dx.*dy, dx.*dz; dy.*dx, 1-3*dy.^2, dy.*dz; dz.*dx, dz.*dy, 1-3*dz.^2];  % I - 3Rhat.Rhat^T but wrong ord: node inds fast, xyz inds slow
    Dker = IRR(inds,inds) .* kron(invRdist.^3,ones(3));  % reord, doublet ker vals
    K(ii3,jj3) = (Sker + eps^2/2*Dker) .* kron(s.w,ones(3,1))';  % wei for ds'
  end
end


%%%%%%%%%
function test_nyst_K_SBT_nonself   % simple demo, not really a math test
verb = 1;
p = 12;                    % order
nspan = 10;                % check conv here
ntpan = 15;                % check conv here
eps = 1e-1;                % large

t = 2*pi*(0:nspan)'/nspan;   % pan param breakpoints (first=0, last=2pi)
rng(0);
t(2:nspan) = t(2:nspan) + 5*(rand(nspan-1,1)-.5)/nspan;  % unequal panels
span = setup_pans(t,p);
t = 2*pi*(0:ntpan)'/ntpan;   % pan param breakpoints (first=0, last=2pi)
t(2:ntpan) = t(2:ntpan) + 5*(rand(ntpan-1,1)-.5)/ntpan;  % unequal panels
tpan = setup_pans(t,p);
N = p*nspan; M = p*ntpan;

Q = eye(3);                     % warm-up: no rot
%[Q,~] = qr(randn(3));           % rand rot mat
[Z1,Zp1] = ellipse_map(1,1,Q,[0;0;0]);   % circle at that orientation
[Z2,Zp2] = ellipse_map(1,1,Q,[3;0;0]);   % translate

span = map_pans(span,Z1,Zp1);
tpan = map_pans(tpan,Z2,Zp2);
if verb, showcurve(span); showcurve(tpan,gca); end

K = nyst_K_SBT_nonself(tpan,span,eps);       % fill it

if verb, figure; ii=3*(1:M)-2; jj=3*(1:N)-2;  % inds of x cmpnts (in and out)
  Kxx = K(ii,jj); imagesc(Kxx); colorbar; title('Kxx'); end

% test applying Fourier modes, npan-convergence, etc...
t = vertcat(span.t)';          % row of params of nodes, src
tt = vertcat(tpan.t)';         % row of params of nodes, targs
mi = 3;                        % input mode
fprintf('For input sin mode freq %d, Fourier output modes are...\n',mi);
f = [sin(mi*t); 0*t; 0*t];     % x-only (pre-rot)
f = Q*f;
u = K*f(:);   % interlace xyz cmpts of f for correct indexing
u = Q'*reshape(u,[3 M]);       % unrot so curve is back in xy-plane
if verb, figure; plot(t,f,'.-'); hold on; plot(tt,u,'.-');
  legend('f_x','f_y','f_z','u_x','u_y','u_z');
  title('test nyst K SBT nonself hitting f_x sin mode'); end
w = vertcat(tpan.w)';          % for quadrature on targ output pans
for m=0:5       % extract cos, sin amplitudes by Euler-Fourier formula...
  am = sum(u.*w.*(ones(3,1)*cos(m*tt)),2)/pi;   % do quadrature (2a_0, a_1, ..)
  if m==0, am=am/2; end        % a_0 case
  fprintf(' u cos m=%d:\t (%19.12g%19.12g%19.12g)\n',m,am(1),am(2),am(3))
  if m>0
    bm = sum(u.*w.*(ones(3,1)*sin(m*tt)),2)/pi;
    fprintf(' u sin m=%d:\t (%19.12g%19.12g%19.12g)\n',m,bm(1),bm(2),bm(3))
  end
end


