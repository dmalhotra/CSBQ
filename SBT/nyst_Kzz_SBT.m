function K = nyst_Kzz_SBT(pan,sbrk,o)
% NYST_KZZ_SBT  discretize zz-cmpnt of SBT K operator, panel-quad xy-plane fiber
%
% K = nyst_Kzz_SBT(pan,sbrk,o) fills a square matrix, a high-order Nystrom
%  discretization of the zz-component of the Keller-Rubinow classical SBT
%  operator which applies in the case of smooth closed centerline in a z=const
%  plane. This is a scalar case, without any "Rhat" or "shat" outer prod terms.
%  It is indep of the fiber radius (which only affects the Lambda(s) "local"
%  term, not included here).
%  We require arc-length parameters of the nodes, arc-length weights,
%  arc-length chart for each panel. Order p is assumed same for each panel.
%
% Inputs:
%  pan - a struct array (one struct per panel) with fields:
%        x = nodes (3*p array)
%        w = weights for arc-length quadrature (length-p col vec)
%        s = arc-length coords of nodes, wrt fixed origin on centerline (")
%        scht = handle to chart from arc-len coord to panel in R3.
%  sbrk - panel-breakpoints in arc-length coords (length npan+1 vec).
%        For now, must have sbrk(1)=0 and sbrk(end)=L, setting the perimeter.
%  o - optional struct with optional fields:
%        o.paux = order of auxiliary nodes (otherwise chosen by default)
% Output:
%  A - N*N Nystrom matrix (where N = total nodes in all panels), that maps
%      values of a L-periodic function f(s) on all nodes to values of
%
%    (Kf)(s) := int_0^L [ f(s')/R(s,s') - f(s)/d(s,s') ] ds'
%
%      at these same nodes, where d(s,s') = (L/pi)sin(pi|s-s'|/L) is the
%      SBT periodized arc-length distance, and ds' is arc-length element.
%      Here f(s) represents z-cmpnt of force density.
% 
% When called with no arguments, does self-test.

% Barnett 1/14/22
if nargin==0, test_nyst_Kzz_SBT; return; end
if nargin<3, o=[]; end

npan = numel(pan);
p = numel(pan(1).w);   % order (assumed same for each pan)

if numel(sbrk)~=npan+1; warning('sbrk should have numel(pan)+1 elements!'); end
L = sbrk(end); if L<=0, error('period must be positive!'); end     % period
N = npan*p;              % total nodes
s = vertcat(pan.s); w = vertcat(pan.w);  % N*1, concat all nodes for easy access
x = horzcat(pan.x);      % 3*N
inpan = @(s,i) mod(s,L)>=sbrk(i) & mod(s,L)<sbrk(i+1);  % peri arc coord in pan?

if ~isfield(o,'paux'), o.paux=2*p; end    % reasonable amount of oversampling
% paux = nodes either side (twice this = # aux nodes per target)
[az aw] = gauss(o.paux);    % discont kernel -> can use plain G-L as aux scheme

K = nan(N,N);          % fill K mat...
Bii = nan(p,N);        % working array for the block row of the 2nd 1/sin term
for ip=1:npan          % target pans (block rows)
  ipl = mod(ip-2,npan)+1; ipr = mod(ip,npan)+1;     % left & right pan inds
  jfar = find((1:npan ~= ipl) & (1:npan ~= ip) & (1:npan ~= ipr));  % far pans
  jjfar = ones(p,1)*(p*(jfar-1)) + (1:p)';         % convert to far node inds
  jjfar = jjfar(:);
  ii = (1:p)+(ip-1)*p;   % block indices in A for target pan
  dx = x(1,ii)' - x(1,jjfar);    % go Cartesian, each cmpnt p-by-(npan-3)p
  dy = x(2,ii)' - x(2,jjfar);
  dz = x(3,ii)' - x(3,jjfar);
  invRdist = 1./sqrt(dx.^2+dy.^2+dz.^2);      % 1/R
  K(ii,jjfar) = invRdist .* w(jjfar)';        % wei for arc-len ds', b'cast row
  % set up 2nd arc-len-only matrix B which will later be row-summed to diag...
  ds = s(ii) - s(jjfar)';                     % p * (npan-3)p
  invsdist = (pi/L) ./ abs(sin(ds*(pi/L)));
  Bii(:,jjfar) = invsdist .* w(jjfar)';    % sin bit & ds wei
  % rest is 3-panel special quad for discont ker, for each targ k in pan...
  ips = {ipl, ip, ipr};        % 3 special pan inds: left-nei, self, right-nei
  s0 = sbrk(ipl); if ip==1, s0=s0-L; end      % special bottom end, unwrap
  s1 = sbrk(ipr+1); if ipr==1, s1=s1+L; end   % special top end, unwrap
  for k=1:p           % loop over nodes in ip'th targ pan
    i = (ip-1)*p + k;        % targ node, w/ 3D loc x(:,i)
    t = s(i);                % targ arc coord
    sc0 = (t-s0)/2; sc1 = (s1-t)/2;     % half-s-lengths to special interval
    sq = [(t+s0)/2 + sc0*az', (s1+t)/2 + sc1*az'];   % aux quad s-nodes, row
    wq = [sc0*aw, sc1*aw];                           % aux quad wei, row
    for sp=1:3               % loop over special sauce panels, fill a 1*p blk...
      jp = ips{sp};          % sauce pan index
      qq = inpan(sq,jp);     % aux inds in this src pan (<- denoted by "qq")
      sqq = sq(qq); wqq = wq(qq);   % rows aux nodes, wei, in this src pan
      sj = pan(jp).s - L*(ip==1 & sp==1) + L*(ip==npan & sp==3);  % wrap src
    %[I,Ip] = interpmat_1d(sqq, sj);  % "aux vals,ders from src node vals" mats
      I = interpmat_1d(sqq, sj);  % "aux vals from src node vals" mat
      xqq = pan(jp).x * I';  % rowwise interpolate to 3D aux locs, 3*numaux
      dr = x(:,i) - xqq;     % aux node R displacement vecs (3*numaux)
      invRdist = 1./sqrt(sum(dr.^2,1));               % row aux node 1/R to targ
      K(i, (jp-1)*p + (1:p)) = (invRdist .* wqq) * I; % aux wei & interp -> row
      ds = t - sqq;          % row aux node arc-dists from targ
      invsdist = (pi/L) ./ abs(sin(ds*(pi/L)));       % periodized s-s' inv-dist
      Bii(k, (jp-1)*p + (1:p)) = (invsdist .* wqq) * I;
    end
    K(i,i) = K(i,i) - sum(Bii(k,:));  % apply 2nd term: row sum to diag entry
  end
end


%%%%%%%%%
function test_nyst_Kzz_SBT
p = 12;                    % order
npan = 10;
tpan = 2*pi*(0:npan)'/npan;   % pan param breakpoints (first=0, last=2pi)
rng(0);
tpan(2:npan) = tpan(2:npan) + 5*(rand(npan-1,1)-.5)/npan;  % unequal panels
pan = setup_pans(tpan,p);

%[Z,Zp] = ellipse_map(1,1,eye(3),zeros(3,1));   % xy-plane unit circle
[Z,Zp] = ellipse_map(1,1);   % generic unit circle
pan = map_pans(pan,Z,Zp);
[pan sbrk] = arccoords_pans(pan);

K = nyst_Kzz_SBT(pan,sbrk);   % do it

ker = @(t,s) 0.5./sin(abs(s-t)/2);   % inv dist btw pts s,t on unit circ
A = nyst_diagdiscont_sca(pan,tpan,ker);
A = A - diag(sum(A,2));     % subtracting rowsum from diag does -f(s)/d(s,s')
norm(K-A, inf)              % we get 10 digits unif over mat els
%figure; subplot(1,3,1); imagesc(K); colorbar; subplot(1,3,2); imagesc(A); colorbar; subplot(1,3,3); imagesc(log10(abs(K./A-1))); colorbar; title('log10 diff')
