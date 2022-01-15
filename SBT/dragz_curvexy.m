% Drag of general xy-plane curve under const (0,0,U) velocity, in SBT.
%
% Needs solve (0,0,f(s)) function in 1D IE SBT[f] = u, just zz-cmpnts, then
% integrate f to get total force F.
% This enables direct comparison with Stokes Dirichlet (vel) BVP.
% Radius is eps, const along curve.
%
% Barnett 1/14/22
clear

mu = 1;               % fluid viscosity
eps = 1e-2;           % radius

a = 2; b = 0.5;
[Z,Zp] = ellipse_map(a,b,eye(3),zeros(3,1));       % ellipse, xy-plane


%a = 0.5;             % inverted ellipse
%zmap = @(t) exp(1i*t) ./ (1+a*exp(2i*t));   % C-plane
%Z = @(t) [real(zmap(t)); imag(zmap(t)); 0*t]; Zp=[];  % xy-plane
%a = 0.4;             % kite, a is narrowness
%Z = @(t) [a*cos(t)+(1-a)*cos(2*t); sin(t); 0*t];
%Zp = @(t) [-a*sin(t)-2*(1-a)*sin(2*t); cos(t); 0*t];

p = 12;                       % order
npan = 20;
tpan = 2*pi*(0:npan)'/npan;   % pan param breakpoints (first=0, last=2pi)
pan = setup_pans(tpan,p);
pan = map_pans(pan,Z,Zp);
[pan sbrk] = arccoords_pans(pan);
L = sbrk(end);
showcurve(pan);

Kzz = nyst_Kzz_SBT(pan,sbrk);
Lambda_zz = 1 - 2*log(pi*eps/(4*L));   % "local" term, scalar, const
