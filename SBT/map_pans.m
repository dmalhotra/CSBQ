function pan = map_pans(pan,Z,Zp)
% MAP_PANS   use chart to add 3D quadrature to parameter panel quadrature
%
% pan = map_pans(pan,Z) expects struct array pan to have column-vector fields
%  t & v giving parametric quadrature nodes & weights for panels on [0,L], and
%  Z function handle (map from row vector of points in [0,L] to matrix of
%  coords in R^3 with 3 rows). It returns pan with added fields:
%     x - 3*p coords in 3D of nodes
%     tx - 3*p coords in 3D of unit tangents at nodes
%     w - length-p column vector of arc-length weights
%  Then {x,w} will be a good quadrature for line integrals on the curve Z([0,L])
%
% pan = map_pans(pan,Z,Zp), if Zp is nonempty, assumes similar function
%  handle Zp(t) = (d/dt) Z(t) and uses this to compute weights and tangents
%  analytically.
%
% Calling without arguments does self-test of both modes (Zp and no Zp).

% Barnett 1/12/22
if nargin==0, test_map_pans; return; end
analderiv = nargin>2 && ~isempty(Zp);
for i=1:numel(pan)
  t = pan(i).t(:)';       % row vec
  pan(i).x = Z(t);
  if analderiv
    xp = Zp(t);
    speed = sqrt(sum(xp.^2,1));   % row vec
    pan(i).tx = xp ./ speed;
    pan(i).w = speed(:).*pan(i).v;
  end
end
if ~analderiv, pan = LIquad_panels(pan); end


%%%%%%%%
function test_map_pans
[Z,Zp,perim_ex] = ellipse_map(1.9,0.7);
p = 12;                    % order
npan = 10;
tpan = 2*pi*(0:npan)'/npan;   % pan param breakpoints (first=0, last=2pi)
rng(0);
tpan(2:npan) = tpan(2:npan) + 5*(rand(npan-1,1)-.5)/npan;  % unequal panels
pan = setup_pans(tpan,p);

pan = map_pans(pan,Z,Zp);    % use Zp case
showcurve(pan)
perim = sum(vertcat(pan.w));
fprintf('analytic Zp:   rel perim err %.3g\n',abs(perim-perim_ex)/perim_ex)

pan = map_pans(pan,Z);       % use only Z node coords (high-ord diff gets w)
perim = sum(vertcat(pan.w));
fprintf('LIquad_panels: rel perim err %.3g\n',abs(perim-perim_ex)/perim_ex)
