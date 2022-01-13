function showcurve(pan,ax)
% SHOWCURVE  plot panel quadrature for a curve in 3D
%
% showcurve(pan) plots in a new figure. Nodes x and (if present) normals nx and
%  tangents tx are shown. radii r are not shown.
%
% showcurve(pan,ax) adds to axes object ax
%
% For a test see: LOAD_GEOM

% Barnett 1/12/22
if nargin<2, figure; else, axes(ax); hold on; end
l = 0.1; % vector scaling

for i=1:numel(pan)
  x = pan(i).x;
  c = 0.8*rand(1,3);  % random RGB
  plot3(x(1,:),x(2,:),x(3,:),'.','markersize',10,'color',c); hold on;
  text(mean(x(1,:)),mean(x(2,:)),mean(x(3,:)),sprintf('%d',i),'color',c);
  if isfield(pan(i),'nx')
    nx = pan(i).nx;
    plot3(x(1,:)+[0;l]*nx(1,:),x(2,:)+[0;l]*nx(2,:),x(3,:)+[0;l]*nx(3,:),'-', 'color',c);
  end
  if isfield(pan(i),'tx')
    tx = pan(i).tx;
    plot3(x(1,:)+[0;l]*tx(1,:),x(2,:)+[0;l]*tx(2,:),x(3,:)+[0;l]*tx(3,:),'-', 'color',c);
  end
end
set(gca,'clipping','off'); axis equal
xlabel('x'); ylabel('y'); zlabel('z');
