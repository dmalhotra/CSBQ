function [Z,Zp,perim] = ellipse_map(a,b,Q,c)
% ELLIPSE_MAP   simple 3D ellipse closed curve analytic map and its perim
%
% [Z,Zp,perim] = ellipse_map(a,b) returns Z a function handle taking row
%  vectors of n parameter values t to 3*n coord matrix in 3D, and Zp = dZ/dt,
%  given a,b semiaxes of ellipse. The orientation and center is generic.
%  perim is its exact perimeter.
%
% [...] = ellipse_map(a,b,Q,c) instead uses translation column vector c and
%  rotation 3*3 matrix Q.

% Barnett 1/12/22
if b>a, [b,a] = deal(a,b); end              % so a>=b
if nargin<3, [Q,~] = qr(rand(3)); end       % rot, generic
if nargin<4, c = rand(3,1); end             % trans, generic
Z = @(t) c + Q*[a*cos(t); b*sin(t); 0*t];   % ellipse in xy plane... rot, trans
Zp = @(t) Q*[-a*sin(t); b*cos(t); 0*t];
k = sqrt(a^2-b^2)/a; [~,E] = ellipke(k^2); perim = 4*a*E;  % ellipse perim
