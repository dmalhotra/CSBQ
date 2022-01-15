function [f fp]  = interpfunc_1d(s,y)
% INTERPFUNC_1D   polynomial interpolation function handle from nodes in 1D
%
% [f fp] = interpfunc_1d(s,y) returns function handle f which evaluates the
%  unique degree-(p-1) polynomial interpolant matching the list of p values y at
%  the list of ordinates s. It also returns fp which evaluates df/ds. The
%  handles also vectorize over an array of ordinates.
%  The source nodes s should be good for polynomial interpolation, ie, have
%  Chebyshev density when their number p (ie, order) is large, and p<40 is
%  advised in all cases. They needn't be well-centered or scaled. f(t) and fp(t)
%  should only be evaluated for t in the minimum interval containing s nodes.
%
%  It is thus a simple non-adaptive "function generator" a la Stein-Blackwell.
  
% Barnett 1/14/22
if nargin==0, test_interpfunc_1d; return; end

cen = (max(s)+min(s))/2; fac = 2/(max(s)-min(s));   % affine s to [-1,1]
s = (s-cen)*fac;
p = numel(s);
co = polyfit(s,y,p-1);                              % coeffs of f
f = @(t) polyval(co,(t-cen)*fac);
cop = fac * (p-1:-1:1) .* co(1:end-1);              % coeffs of f' (rescale!)
fp = @(t) polyval(cop,(t-cen)*fac);


%%%%%%%%
function test_interpfunc_1d
off = 1.3;           % test centering and scaling
sc = 0.03;           % -> roundoff err factor 1/sc for val, 1/sc^2 for deriv
x = sc*gauss(16) + off;
f = @(x) sin(x/sc + 0.7);
fp = @(x) (1/sc)*cos(x/sc + 0.7);
data = f(x);            % func on smooth (src) nodes
t = sc*(2*rand(1000,1) - 1) + off;    % targs cover same interval as the x lie
uex = f(t);     % col vec
upex = fp(t);     % col vec
[F Fp] = interpfunc_1d(x,data);
u = F(t);
up = Fp(t);
fprintf('max abs err interpfunc_1d in [a,b] : %.3g (val), %.3g (deriv)\n',norm(u-uex,inf),norm(up-upex,inf))

% Note matlab's not shifting & rescaling:
% p=10; s=rand(p,1)*10+100; y=sin(s);
% P=polyfit(s,y,p-1)
% norm(y-polyval(P,s),inf)      % err 1e-3 :(
