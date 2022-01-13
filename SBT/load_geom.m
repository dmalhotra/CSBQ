function pan = load_geom(fnam)
% LOAD_GEOM  read in 3D centerline panel quadrature from a Dhairya geom file
%
% pan = load_geom(fnam) creates a struct array of panels from the
%  contents of file fnam. All panels have the same number of nodes p.
%  The i'th struct pan(i) has fields:
%     x : 3*p array of 3D node coords
%     r : length-p list of slender body radii at nodes
%     nx : 3*p array of some vectors normal to the curve (unused)
%
% Called with no arguments does self-test.

% Barnett 1/12/22
if nargin==0, test_load_geom; return; end

fid = fopen(fnam,'r');
a = textscan(fid,'%f','CommentStyle','#');
fclose(fid);
% done with reading; now decode the cell array a...
a = a{1};
p = a(8);   % order for all panels
ford = a(9);
ipandat = 7*p+2;
npan = numel(a)/ipandat;
if npan~=round(npan), error('file does not have expected number of entries!'); end
pan = struct([]);
for i=1:npan
  x = nan(3,p); x(:,1) = a((i-1)*ipandat + (1:3));
  for j=2:p, x(:,j) = a((i-1)*ipandat + 7*(j-1) + 2+(1:3)); end
  r = nan(1,p); r(1) = a((i-1)*ipandat + 4);
  for j=2:p, r(j) = a((i-1)*ipandat + 7*(j-1) + 2+4); end
  nx = nan(3,p); nx(:,1) = a((i-1)*ipandat + (5:7));
  for j=2:p, nx(:,j) = a((i-1)*ipandat + 7*(j-1) + 2+(5:7)); end
  pan(i).x = x;  pan(i).r = r; pan(i).nx = nx;
end


%%%%%%
function test_load_geom
pan = load_geom('../data/geom.data');
showcurve(pan); title('test load\_geom');

