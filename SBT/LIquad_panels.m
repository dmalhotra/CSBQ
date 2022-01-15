function pan = LIquad_panels(pan)
% LIQUAD_PANELS  add arc-length quadrature weights, etc, given only 3D nodes
%
% pan = LIquad_panels(pan) given pan, a struct arrays of panels with
%  node fields x, assumes each panel nodes come from some smooth map from
%  Gauss-Legendre quadrature nodes for the standard panel [-1,1]. It
%  uses p-order differentiation to get arc-length (s) quadrature weights (w),
%  unit tangents (tx), and speeds (sp) wrt map from standard panel.
%  These fields added to the output panel struct array; others are preserved.
%  Then {x,w} is a good quadrature for line integrals on curve approximated
%  by the input nodes.
%
% Note: it's peculiar that sp (Jacobian w.r.t. unknown map) can be computed!
%
% Without arguments does self-test (look at 2nd output of it).

% Barnett 1/12/22
if nargin==0, map_pans; return; end        % uses another self-test

p = size(pan(1).x,2);              % nodes per panel = "order"
[z q] = gauss(p);                  % std panel quadr on [-1,1]
[~,D] = interpmat_1d(z,z);         % high-order differentiation mat on std nodes

for i=1:numel(pan)
  v = pan(i).x * D';   % differentiate rows as if values on std nodes in [-1,1]
  sp = sqrt(sum(v.^2,1));          % row of speeds (Jacobian factor) at nodes
  pan(i).tx = v ./ sp;
  pan(i).w = (q .* sp)';           % q is row. Use Jacobian for weight w.r.t. s
  pan(i).sp = sp(:);               % save Jacobians
  sp
end
