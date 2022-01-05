function [L Lp]  = interpmat_1d(t,s)
% INTERPMAT_1D   interpolation matrices: nodes in 1D to target nodes & derivs
%
% L = interpmat_1d(t,s) returns interpolation matrices taking values on nodes s
%  to values (L) and first derivs (Lp) at target nodes t. The source nodes s
%  should be good for polynomial interpolation, ie, have Chebyshev density when
%  their number p is large.
%
%  Note: Computed in Helsing style, thus stable up to about p~40.

% bits taken from qplp/interpmatrix.m from 2013.
% Barnett 7/17/16. Auto-centering & scaling for stability 12/23/21.
% deriv 1/5/22.

if nargin==0, test_interpmat_1d; return; end

cen = (max(s)+min(s))/2; hwid = (max(s)-min(s))/2;   % affine s to [-1,1]
s = (s-cen)/hwid;
t = (t-cen)/hwid;

p = numel(s); q = numel(t); s = s(:); t = t(:);       % all col vecs
n = p; % set the polynomial order we go up to (todo: check why bad if not p)
V = ones(p,n); for j=2:n, V(:,j) = V(:,j-1).*s; end   % polyval matrix on nodes
R = ones(q,n); for j=2:n, R(:,j) = R(:,j-1).*t; end   % polyval matrix on targs
L = (V'\R')'; % backwards-stable way to do it (Helsing) See corners/interpdemo.m
R = ones(q,n); for j=2:n, R(:,j) = R(:,j-1).*t; end   % deriv matrix on targs
R = [zeros(q,1), R(:,1:end-1)];  % convert R into polyderiv eval matrix
for j=3:n, R(:,j) = (j-1)*R(:,j); end
Lp = (V'\R')';
Lp = (1/hwid)*Lp;     % undo ordinate scale fac


%%%%%%
function test_interpmat_1d
off = 1.3;           % test centering and scaling
sc = 0.03;
x = sc*gauss(16) + off;
f = @(x) sin(x/sc + 0.7);
fp = @(x) (1/sc)*cos(x/sc + 0.7);
data = f(x);            % func on smooth (src) nodes
t = sc*(2*rand(1000,1) - 1) + off;    % targs cover same interval as the x lie
uex = f(t);     % col vec
upex = fp(t);     % col vec
[L Lp] = interpmat_1d(t,x);
u = L * data;
up = Lp * data;
fprintf('max abs err for interp in [a,b] : %.3g (val), %.3g (deriv)\n',max(abs(u - uex)),max(abs(up - upex)))
fprintf('interp mat inf-norm = %.3g;   max element size = %.3g\n',norm(L,inf),max(abs(L(:))))
