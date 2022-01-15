function [pan sbrk] = arccoords_pans(pan)
% ARCCOORDS_PANS   node and breakpt arc-length coords given 3D nodes, std speeds
%
% [pan sbrk] = arccoords_pans(pan) uses the 3D node coords field x, and the
%  speeds wrt std map field sp, in the pan struct array (assuming each panel
%  is a smooth map from standard panel nodes on [-1,1]) to add fields to each
%  element of pan:
%      s    - (length-p column vector) arc-length coords of nodes on this panel
%      scht - (function handle taking length-n row vec to 3*n matrix) chart from
%             arc-length coords to this panel in 3D.
%      schtp - similar except s-deriv of chart, giving hat{s} vector in 3D.
%  Here, arc-length coords are relative to a fixed origin on the curve: the
%  lower end of the first panel. Also returns sbrk, a npan+1 column vector of
%  panel breakpoints in the same arc-length coords, where npan=numel(pan).
%
% Notes: uses spectral method linear solve for antiderivative.
%
% Without arguments does self-test.

% Barnett 1/13/22. Chart handles & test 1/14/22, unused.
if nargin==0, test_arccoords_pans; return; end

p = size(pan(1).x,2);              % nodes per panel = "order"
[z q] = gauss(p);                  % std panel quadr on [-1,1]
[L,D] = interpmat_1d([-1;z],z);    % interp and diff mats to -1 and std nodes

J = horzcat(pan.sp);               % stack all Jacobians to p-by-npan matrix
A = [L(1,:); D(2:end,:)];          % lin sys mat: row1 imposes s(-1)=0, others
                                   % that Ds = a col of J to match deriv.
npan = numel(pan);
S = A \ [zeros(1,npan); J];        % overdet lin sys for antiderivs of cols of J
arclens = sum(J .* q(:),1);        % row vec of sums of arclen wei on each panel
sbrk = [0; cumsum(arclens')];
for i=1:npan
  s = sbrk(i) + S(:,i);            % write out the s node arc coords w/ offsets
  pan(i).s = s;
  [x xp]  = interpfunc_1d(s,pan(i).x(1,:));   % Cartesian handles for charts...
  [y yp]  = interpfunc_1d(s,pan(i).x(2,:));
  [z zp]  = interpfunc_1d(s,pan(i).x(3,:));
  pan(i).scht = @(s) [x(s); y(s); z(s)];
  pan(i).schtp = @(s) [xp(s); yp(s); zp(s)];
end


%%%%%%%%%%%%
function test_arccoords_pans
p = 12;                    % order
npan = 10;
tpan = 2*pi*(0:npan)'/npan;   % pan param breakpoints (first=0, last=2pi)
rng(0);
tpan(2:npan) = tpan(2:npan) + 5*(rand(npan-1,1)-.5)/npan;  % unequal panels
pan = setup_pans(tpan,p);

% basic unit circle test where s param = t param
[Z,Zp] = ellipse_map(1,1);
pan0 = map_pans(pan,Z,Zp);       % don't overwrite pan since need for below
[pan0 sbrk] = arccoords_pans(pan0);
t = vertcat(pan.t);              % true input param nodes
s = vertcat(pan0.s);             % param nodes extracted to test
fprintf('arccoord_pans:\tcircle param node reconstruction max err %.3g\n',norm(s-t,inf))
tt = 0.9 + 0.5*rand(1,10);      % param test points all in same panel
ii = sum(tt>tpan,1);            % which panels they're in
if norm(diff(ii))>0, error('test points not in same pan!'); end
i = ii(1);
fprintf('\t pointwise test err in scht %.3g, schtp %.3g\n',norm(pan0(i).scht(tt) - Z(tt),inf), norm(pan0(i).schtp(tt) - Zp(tt),inf))  % check scht on row vec

% fancier test where we don't know true arc-len param of ellipse, but eyeball it
[Z,Zp,perim_ex] = ellipse_map(1.9,0.7);
pan = map_pans(pan,Z,Zp);    % use analytic Zp
[pan sbrk] = arccoords_pans(pan);
perim = sbrk(end);              % we don't know how to test much else here
fprintf('arccoord_pans:\tellipse rel perim err %.3g\n',abs(perim-perim_ex)/perim_ex)
figure; plot(vertcat(pan.s),zeros(p*npan,1),'k.'); vline(sbrk); axis tight;
title('arccoords from 3D nodes and std speeds: should look like valid pans')
