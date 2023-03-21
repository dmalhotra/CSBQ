function K = nyst_K_SBT(pan,sbrk,o)
% NYST_K_SBT  discretize tensor SBT K operator on panel-quad closed fiber
%
% K = nyst_K_SBT(pan,sbrk,o) fills a square matrix, a high-order Nystrom
%  discretization of Keller-Rubinow classical SBT tensor
%  operator which applies in the case of smooth closed centerline in R3.
%  It is indep of the fiber radius (which only affects the Lambda(s) "local"
%  term, not included here).
%  We require arc-length parameters of the nodes and their arc-length weights.
%  Order p is assumed same for each panel.
%
% Inputs:
%  pan - a struct array (one struct per panel) with fields:
%        x = nodes (3*p array)
%        tx = unit tangents ("\hat{s}") at nodes (3*p array)
%        w = weights for arc-length quadrature (length-p col vec)
%        s = arc-length coords of nodes, wrt fixed origin on centerline (")
%        ipl, ipr = indices of panels (in pan struct array) to left and right
%  sbrk - panel-breakpoints in arc-length coords (length npan+1 vec).
%        For now, must have sbrk(1)=0 and sbrk(end)=L, setting the perimeter.
%  o - optional struct with optional fields:
%        o.paux = order of auxiliary nodes (otherwise chosen by default)
% Output:
%  K - 3N*3N Nystrom matrix (where N = total nodes in all panels), that maps
%      values of a L-periodic vector function f(s) on all nodes to values of
%
%    (Kf)(s) := int_0^L [ (I + Rhat.Rhat^T)/R(s,s') f(s')
%                        -(I + shat.shat^T)/d(s,s') f(s) ] ds'
%
%      at these same nodes, where d(s,s') = (L/pi)sin(pi|s-s'|/L) is the
%      SBT periodized arc-length distance, ds' is arc-length element, I
%      is the 3x3 identity tensor, and hat indicates unit vector.
%      Rhat := R/|R| where R is the displacement x(s)-x(s').
%      shat := x'(s)/|x'(s)| is the tangent at s (not s', strangely).
%      f(s) represents force density, and Kf(s) scaled velocity.
%      Storage order is "interleaved", ie coords {xyz} fast (innermost),
%      nodes slow, ie K is a N*N block matrix with each entry a 3*3 matrix.
%
% See also: NYST_KZZ_SBT which is the scalar special case.
% 
% When called with no arguments, does self-test.
  
% Barnett 1/17/22
if nargin==0, test_nyst_K_SBT; return; end
if nargin<3, o=[]; end

npan = numel(pan);
p = numel(pan(1).w);   % order (number of nodes per pan, assumed same for each pan)

if numel(sbrk)~=npan+1; warning('sbrk should have numel(pan)+1 elements!'); end
L = sbrk(end); if L<=0, error('period must be positive!'); end     % period
N = npan*p;              % total nodes
s = vertcat(pan.s); w = vertcat(pan.w);  % N*1, concat all nodes for easy access
x = horzcat(pan.x);      % 3*N
inpan = @(s,i) mod(s,L)>=sbrk(i) & mod(s,L)<sbrk(i+1);  % peri arc coord in pan?
inds3 = @(n) reshape(reshape(1:3*n, [n 3])',[1 3*n]);   % ind list w/ xyz fast

if ~isfield(o,'paux'), o.paux=2*p; end    % reasonable amount of oversampling
% paux = nodes either side (twice this = # aux nodes per target)
[az aw] = gauss(o.paux);    % discont kernel -> can use plain G-L as aux scheme

K = nan(3*N,3*N);      % fill K mat
Bii = nan(3*p,3*N);    % working array for the block row of the 2nd 1/sin term
for ip=1:npan          % target pans (block rows)
  ipl = pan(ip).ipl; ipr = pan(ip).ipr;             % left & right pan inds
  jfar = find((1:npan ~= ipl) & (1:npan ~= ip) & (1:npan ~= ipr));  % far pans
  jjfar = ones(p,1)*(p*(jfar-1)) + (1:p)';          % convert to far node inds
  jjfar = jjfar(:);
  nfar = numel(jjfar);
  jjfar3 = ones(3*p,1)*(3*p*(jfar-1)) + (1:3*p)';   % far input inds
  jjfar3 = jjfar3(:);
  ii = (1:p)+(ip-1)*p;   % block node indices for target pan
  ii3 = (1:3*p)+(ip-1)*3*p;   % this pan output (row) indices
  dx = x(1,ii)' - x(1,jjfar);    % go Cartesian, each cmpnt p-by-(npan-3)p
  dy = x(2,ii)' - x(2,jjfar);
  dz = x(3,ii)' - x(3,jjfar);
  invRdist = 1./sqrt(dx.^2+dy.^2+dz.^2);            % 1/R matrix, p*nfar
  dx = dx.*invRdist; dy = dy.*invRdist; dz = dz.*invRdist; % (dx,dy,dz) now Rhat
  IRR = [1+dx.^2, dx.*dy, dx.*dz; dy.*dx, 1+dy.^2, dy.*dz; dz.*dx, dz.*dy, 1+dz.^2];  % I + Rhat.Rhat^T but wrong ord: node inds fast, xyz inds slow
  ker = IRR(inds3(p),inds3(nfar)) .* kron(invRdist,ones(3));  % reord, vals
  K(ii3,jjfar3) = ker .* kron(w(jjfar),ones(3,1))';   % wei for ds', b'cast row
  % set up 2nd arc-len-only matrix B which will later be row-summed to diag...
  ds = s(ii) - s(jjfar)';                     % p * (npan-3)p
  invsdist = (pi/L) ./ abs(sin(ds*(pi/L)));
  sx = pan(ip).tx(1,:)';  % Cartesian cols of shat at the p targ nodes for pan
  sy = pan(ip).tx(2,:)';
  sz = pan(ip).tx(3,:)';
  Iss = [1+sx.^2, sx.*sy, sx.*sz; sy.*sx, 1+sy.^2, sy.*sz; sz.*sx, sz.*sy, 1+sz.^2];  % I + shat.shat^T but row ord wrong: node inds fast, xyz inds slow
  ker = kron(ones(1,nfar), Iss(inds3(p),:)) .* kron(invsdist,ones(3));   % reord, cp over src, 1/d
  Bii(:,jjfar3) = ker .* kron(w(jjfar),ones(3,1))';    % sin bit & ds wei
  % rest is 3-panel special quad for discont ker, for each targ k in pan...
  ips = {ipl, ip, ipr};        % 3 special pan inds: left-nei, self, right-nei
  s0 = sbrk(ipl); if ip==1, s0=s0-L; end      % special bottom end, unwrap
  s1 = sbrk(ipr+1); if ipr==1, s1=s1+L; end   % special top end, unwrap
  for k=1:p           % loop over nodes in ip'th targ pan
    i = (ip-1)*p + k;        % targ node, w/ 3D loc x(:,i)
    i3 = 3*(i-1)+(1:3); k3 = 3*(k-1)+(1:3);  % xyz Cartesian blk inds for 1 node
    t = s(i);                % targ arc coord
    shat = pan(ip).tx(:,k); Iss = eye(3)+shat*shat';  % hat{s} for i; 3x3 @ targ
    sc0 = (t-s0)/2; sc1 = (s1-t)/2;     % half-s-lengths to special interval
    sq = [(t+s0)/2 + sc0*az', (s1+t)/2 + sc1*az'];   % aux quad s-nodes, row
    wq = [sc0*aw, sc1*aw];                           % aux quad wei, row
    for sp=1:3               % loop over special sauce panels, fill a 1*p blk...
      jp = ips{sp};          % sauce pan index
      qq = inpan(sq,jp);     % aux inds in this src pan (<- denoted by "qq")
      sqq = sq(qq); wqq = wq(qq);   % rows aux nodes, wei, in this src pan
      sj = pan(jp).s - L*(ip==1 & sp==1) + L*(ip==npan & sp==3);  % wrap src
      I = interpmat_1d(sqq, sj);    % "aux vals from src node vals" mat
      xqq = pan(jp).x * I';  % rowwise interpolate to 3D aux locs, 3*numaux
      dr = x(:,i) - xqq;     % aux node R displacement vecs (3*numaux)
      invRdist = 1./sqrt(sum(dr.^2,1));               % row aux node 1/R to targ
      dx = dr(1,:); dy = dr(2,:); dz = dr(3,:);   % its Cartesians
      dx = dx.*invRdist; dy = dy.*invRdist; dz = dz.*invRdist;  % now Rhat
      IRR = [1+dx.^2, dx.*dy, dx.*dz; dy.*dx, 1+dy.^2, dy.*dz; dz.*dx, dz.*dy, 1+dz.^2];  % I + Rhat.Rhat^T but with col ord aux node inds fast, xyz inds slow
      ker = IRR(:,inds3(numel(wqq))) .* kron(invRdist,ones(3));   % reord, vals
      jjnr3 = (jp-1)*3*p + (1:3*p);   % col inds for src pan
      K(3*(i-1)+(1:3), jjnr3) = (ker .* kron(wqq,ones(1,3))) * kron(I,eye(3)); % aux wei & interp -> 3 rows
      ds = t - sqq;          % row aux node arc-dists from targ
      invsdist = (pi/L) ./ abs(sin(ds*(pi/L)));       % periodized s-s' inv-dist
      ker = kron(invsdist, Iss);
      Bii(k3, jjnr3) = (ker .* kron(wqq,ones(1,3))) * kron(I,eye(3));
    end
    Biiksums = reshape(sum(reshape(Bii(k3,:),[9 N]),2),[3 3]);  % sum each of 9
    K(i3,i3) = K(i3,i3) - Biiksums;  % apply 2nd term: outer prod coord sums to 3x3 diag entries
  end
end


%%%%%%%%%
function test_nyst_K_SBT   % simple test for rot-invariance wrt scalar Kzz
verb = 0;
p = 12;                    % order
npan = 10;                 % check conv here (>20 for ellipse)
tpan = 2*pi*(0:npan)'/npan;   % pan param breakpoints (first=0, last=2pi)
rng(0);
tpan(2:npan) = tpan(2:npan) + 5*(rand(npan-1,1)-.5)/npan;  % unequal panels
pan = setup_pans(tpan,p);
N = p*npan;

[Q,~] = qr(randn(3));           % rand rot mat
%Q = eye(3);                     % warm-up: no rot
%th = 2*pi*rand; s=sin(th); c=cos(th);
%Q = [c s 0; -s c 0; 0 0 1]'; i = [1 3 2]; Q = Q(i,i);     % rots in xz only

% choose either: circle tests Fourier modes m+-2; ellipse merely convergence...
%[Z,Zp] = ellipse_map(2,.5,Q,zeros(3,1));   % ellipse at that orientation
[Z,Zp] = ellipse_map(1,1,Q,zeros(3,1));   % circle at that orientation

pan = map_pans(pan,Z,Zp);
[pan sbrk] = arccoords_pans(pan);
%norm(sbrk-tpan)  % for unit circle only
if verb, showcurve(pan); end

K = nyst_K_SBT(pan,sbrk);       % do it, 3N*3N
if verb, figure; ii=3*(1:N)-2;  % inds of x cmpnts (for either in or out)
  Kii = K(ii,ii); imagesc(Kii-diag(diag(Kii))); colorbar; title('Kxx offdiag'); end  % disturbingly oscillatory weights

Kzz = nyst_Kzz_SBT(pan,sbrk);   % scalar version N*N (didn't care about Q)
d = Q(:,3);                     % vec that z vec rot to: f and u parallel to it
Pd = kron(eye(N),d);            % proj along d, from 3N Cart coords to N outputs
Kdd = Pd' * (K * Pd);           % proj input and output sides
fprintf('max diff, scalar Kzz vs proj tensor K: %.3g\n', norm(Kdd-Kzz, inf))

% test applying Fourier modes, npan-convergence, etc...
t = vertcat(pan.t)';           % row of params of nodes
mi = 3;                        % input mode
fprintf('For input sin mode freq %d, Fourier output modes are...\n',mi);
f = [sin(mi*t); 0*t; 0*t];     % x-only (pre-rot)
%f = [0*t;0*t;sin(mi*t)];       % z-only (pre-rot)
f = Q*f;
u = K*f(:);   % interlace xyz cmpts of f for correct indexing
u = Q'*reshape(u,[3 N]);       % unrot so curve is back in xy-plane
if verb, figure; plot(t,[f;u],'.-'); legend('f_x','f_y','f_z','u_x','u_y','u_z'); title('test nyst K SBT hitting f_x sin mode'); end
w = vertcat(pan.w)';           % for quadrature
for m=0:10       % extract cos, sin amplitudes by Euler-Fourier formula...
  am = sum(u.*w.*(ones(3,1)*cos(m*t)),2)/pi;   % do quadrature (2a_0, a_1, ..)
  if m==0, am=am/2; end        % a_0 case
  fprintf(' u cos m=%d:\t (%19.12g%19.12g%19.12g)\n',m,am(1),am(2),am(3))
  if m>0
    bm = sum(u.*w.*(ones(3,1)*sin(m*t)),2)/pi;
    fprintf(' u sin m=%d:\t (%19.12g%19.12g%19.12g)\n',m,bm(1),bm(2),bm(3))
  end
end


