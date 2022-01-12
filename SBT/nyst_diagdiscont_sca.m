function A = nyst_diagdiscont_sca(pan,tpan,ker,o)
% NYST_DIAGDISCONT_SCA  Nystrom discretize diag-discontinuous kernel, 1D panels
%
% A = nyst_diagdiscont_sca(pan,tpan,ker) returns square matrix giving Nystrom
%  discretization of the operator with kernel function ker on the 1D
%  periodic panelization given by breakpoints tpan and panel quadratures
%  pan. For self- and neighbor-interactions between panels, a special scheme
%  uses Gauss-Legendre auxiliary nodes covering 3 panels (left neighbor, self,
%  and right neighbor). This is high-order accurate for kernels that are smooth
%  in s<t and in s>t, but may be discontinuous on the diagonal s=t. Any
%  period (called L below) may be used.
%
%  The suffix _sca refers to scalar-valued functions and operators.
%
% A = nyst_diagdiscont_sca(pan,tpan,ker,o) also controls options.
%
% Inputs:
%  pan - a struct array (one struct per panel) with fields w (column vec of
%        p weights), t (column vector of p nodes), c (center of panel), sc
%        (half-size of panel). Order p is assumed same for all panels for now.
%  tpan - list (npan+1 length) of panel breakpoints. The first must be 0
%        and the last (npan+1)th must be the period L, which it sets, for now
%        (simplifies wrapping).
%  ker - kernel function handle of (s,t). Must vectorize if s is a column
%        vector and t a row vector to give a rectangular array output.
%  o - optional struct with optional fields:
%        o.paux = order of auxiliary nodes (otherwise chosen by default)
% Output:
%  A - N*N Nystrom matrix (where N = total nodes in all panels), that maps
%      values of a L-periodic function u(t) on all nodes to values of
%      (Ku)(t) := int_0^L k(t,s) u(s) ds at these same nodes.
% 
% When called with no arguments, does self-test.

% Barnett 12/23/21. general period 1/11/22.
% Todo: variable order per panel? General breakpoints & get periodicity right?

if nargin==0, test_nyst_diagdiscont_sca; return; end
if nargin<4, o=[]; end
addpath utils

npan = numel(pan);
p = numel(pan(1).w);   % order (assumed same for each pan)
if numel(tpan)~=npan+1; warning('tpan should have numel(pan)+1 elements!'); end
L = tpan(end); if L<=0, error('period must be positive!'); end     % period
N=npan*p;
s.t = vertcat(pan.t); s.w = vertcat(pan.w);  % concat all nodes for easy access
inpan = @(t,i) mod(t,L)>=tpan(i) & mod(t,L)<tpan(i+1);  % is param in pan?

if ~isfield(o,'paux'), o.paux=2*p; end    % reasonable amount of oversampling
% paux = nodes either side (twice this = # aux nodes per target)
[az aw] = gauss(o.paux);    % discont kernel -> can use plain G-L as aux scheme

A = nan(N,N); % fill A...
for ip=1:npan          % target pans (block rows)
  ipl = mod(ip-2,npan)+1; ipr = mod(ip,npan)+1;     % left & right pan inds
  jfar = find((1:npan ~= ipl) & (1:npan ~= ip) & (1:npan ~= ipr));  % far pans
  jjfar = ones(p,1)*(p*(jfar-1)) + (1:p)';         % convert to far node inds
  jjfar = jjfar(:);
  ii = (1:p)+(ip-1)*p;   % block indices in A for target pan
  A(ii,jjfar) = ker(s.t(ii),s.t(jjfar)') .* s.w(jjfar)';   % broadcast to rows
  % rest is 3-panel special quad for discont ker, for each targ k in pan...
  ips = {ipl, ip, ipr};        % 3 special pan inds: left-nei, self, right-nei
  t0 = tpan(ipl); if ip==1, t0=t0-L; end      % special bottom end, unwrap
  t1 = tpan(ipr+1); if ipr==1, t1=t1+L; end    % special top end, unwrap
  for k=1:p           % nodes in ip'th targ pan
    i = (ip-1)*p + k;        % targ node
    t = s.t(i);              % targ param
    sc0 = (t-t0)/2; sc1 = (t1-t)/2;    % half-param-lengths to special interval
    tq = [(t+t0)/2 + sc0*az', (t1+t)/2 + sc1*az'];   % aux quad node params, row
    wq = [sc0*aw, sc1*aw];                       % aux quad wei, row
    for sp=1:3               % loop over special sauce panels
      jp = ips{sp};          % src pan ind
      qq = inpan(tq,jp);     % aux inds in this src pan
      tj = pan(jp).t - L*(ip==1 & sp==1) + L*(ip==npan & sp==3);  % wrap
      I = interpmat_1d(tq(qq), tj);     % aux vals from src node vals
      A(i, (jp-1)*p + (1:p)) = (ker(t,tq(qq)) .* wq(qq)) * I;    % aux wei
    end
  end
end


%%%%%%
function test_nyst_diagdiscont_sca
p = 12;                  % order
npan = 16;
L = 2*pi;                   % period (2pi for trig tests below)
tpan = L*(0:npan)/npan;    % pan param breakpoints (first=0, last=2pi)
rng(0);
tpan(2:npan) = tpan(2:npan) + 0.2*L*(2*rand(1,npan-1)-1)/npan;  % unequal panels
fprintf('npan=%d, p=%d, pan lengths in range [%.3g,%.3g]:\n',npan,p,min(diff(tpan)),max(diff(tpan)))
pan = setup_pans(tpan,p);
s.t = vertcat(pan.t);                % get all nodes

disp('smooth kernel...')
ker = @(t,s) cos(t-s);               % a smooth degenerate kernel (rank=2)
A = nyst_diagdiscont_sca(pan,tpan,ker);
u = sin(s.t+0.7); norm(A*u-pi*u)     % test eigfunc, eigval = pi
u = sin(3*s.t+1.1); norm(A*u)        % test a null vector
u = sin(7*s.t-1); norm(A*u)          % test more osc null vec

disp('diag-discont kernel...')
ker = @(t,s) cos(mod(s-t,2*pi)/2);   % discont kernel (half-cycle of cos)
A = nyst_diagdiscont_sca(pan,tpan,ker);
u = 1+0*s.t; norm(A*u)               % test the null vector (too easy)
u = exp(1i*s.t); norm(A*u-1i*8/3*u)  % test eigfunc, eigval = 8i/3, guessed
lam7 = 1i * 56/195;                  % empirical (guessed by repeating digits)
u = exp(7i*s.t); norm(A*u-lam7*u)    % test eigfunc

%figure; imagesc(A); axis equal tight; colorbar; title('A matrix')

