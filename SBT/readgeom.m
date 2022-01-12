% read Dhairya's geom data format for Cheby panels in R3.
% Barnett 9/18/21

% Note: file format doesn't contain
%   * underlying parameter values for the panels (assumed equispaced in param)
%   * speed function at nodes ... can get by interpolation?

clear
fnam='../data/geom.data';
fnam='../data/tangle-adap.geom';
fid = fopen(fnam,'r');
a = textscan(fid,'%f','CommentStyle','#');
fclose(fid);
a = a{1};
p = a(8);
ford = a(9);
ipandat = 7*p+2;
npan = numel(a)/ipandat;
if npan~=round(npan), error('file not expected number of entries!'); end
pa = cell(npan,1);
figure(1); clf;
for i=1:npan
  x = nan(3,p); x(:,1) = a((i-1)*ipandat + (1:3));
  for j=2:p, x(:,j) = a((i-1)*ipandat + 7*(j-1) + 2+(1:3)); end
  r = nan(1,p); r(1) = a((i-1)*ipandat + 4);
  for j=2:p, r(j) = a((i-1)*ipandat + 7*(j-1) + 2+4); end
  nx = nan(3,p); nx(:,1) = a((i-1)*ipandat + (5:7));
  for j=2:p, nx(:,j) = a((i-1)*ipandat + 7*(j-1) + 2+(5:7)); end
  pa{1}.x = x;  pa{1}.r = r; pa{1}.nx = nx;
  plot3(x(1,:),x(2,:),x(3,:),'k.','markersize',20); hold on;
  l = 0.1; plot3(x(1,:)+[0;l]*nx(1,:),x(2,:)+[0;l]*nx(2,:),x(3,:)+[0;l]*nx(3,:),'b-');
  text(mean(x(1,:)),mean(x(2,:)),mean(x(3,:)),sprintf('%d',i));
end
set(gca,'clipping','off'); axis vis3d equal; drawnow
