function pan = LIquad_panels(pan)
% LIQUAD_PANELS  add arc-length quadrature weights to panels given 3D nodes
%
% pan = LIquad_panels(pan) given pan, a struct arrays of panels with
%  node fields x, assumes Gauss-Legendre quadrature (standard panel) on each,
%  uses p-order differentiation to get arc-length (s) quadrature weights w
%  and unit tangents tx, which are fields added to the output panel struct
%  array. Other fields are preserved.
%
% Without arguments does self-test (look at 2nd output of it).

% Barnett 1/12/22
if nargin==0, map_pans; return; end        % uses another self-test

% set up standard panel on [-1,1]...
p = size(pan(1).x,2);              % nodes per panel = "order"
[z q] = gauss(p);
[~,D] = interpmat_1d(z,z);         % high-order differentiation mat on std nodes

for i=1:numel(pan)
  v = pan(i).x * D';   % differentiate rows as if values on std nodes in [-1,1]
  sp = sqrt(sum(v.^2,1));          % row of speeds at nodes
  pan(i).tx = v ./ sp;
  arclen = sum(q .* sp);           % q is row
  pan(i).w = (arclen/2) * q(:);    % rescale weights to be w.r.t. arc length
end
