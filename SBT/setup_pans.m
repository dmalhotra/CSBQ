function pan = setup_pans(tpan,p)
% SETUP_PANS   parameter panel quadrature struct array from 1D breakpoints
%
% pan = setup_pans(tpan,p) returns a struct array with fields:
%  t (column of nodes), v (column of weights), c (center),
%  sc (half-length, ie scale factor)
%  given tpan (a row or column of npan+1 breakpoints) and p (a panel order).
%  Gauss-Legendre rules are used on each panel. The number of elements of
%  pan is npan. No periodicity is assumed, but if tpan(1)=0 and tpan(end)=2*pi
%  then it is appropriate for 2pi-periodic panelization. This is just a
%  quadrature w.r.t. parameter t; no map to 3D geom is done.
%
% Without arguments, a self-test is done.

% Barnett 12/23/21. 1/12/22 changed quad-wei to v.
if nargin==0, test_setup_pans; return; end

[z w] = gauss(p);   % one std panel
npan = numel(tpan)-1;
for i=1:npan
  pan(i).sc = (tpan(i+1)-tpan(i))/2;
  pan(i).v = pan(i).sc*w(:);       % pan is struct array. v is column
  pan(i).c = (tpan(i+1)+tpan(i))/2;       % center param of pan
  pan(i).t = pan(i).c + pan(i).sc*z(:);
end

%%%%%%%
function test_setup_pans
npan = 16;
p = 10;
tpan = 2*pi*(0:npan)/npan;    % pan param breakpoints (first=0, last=2pi)
tpan(5) = tpan(5) + 0.2;      % check can handle unequal pans
pan = setup_pans(tpan,p);
s.t = vertcat(pan.t); s.w = vertcat(pan.v);     % concat all pan nodes
int1 = sum(s.w);              % integral of 1
err = int1-2*pi
assert(abs(err)<1e-14)

figure; plot([s.t s.t]',[0*s.w s.w]','b.-'); vline(tpan);
xlabel('t'); title('setup\_pans test'); axis tight;
