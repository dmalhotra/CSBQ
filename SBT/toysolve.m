% solve (ie "invert") toy scalar SBT IE to get f (force) given u (velocity),
% discretized w/ uniform panels in one 2pi-periodic variable t.
% The kernel is general, but the speed function is 1.
% Note: in real SBT will need to interpolate speed func.
% Barnett 12/17/21
clear; verb = 1;

npan = 16;
p = 10;
[z w] = gauss(p);   % one std panel
tpan = 2*pi*(0:npan)/npan;    % pan param breakpoints
                              % tpan(5) = tpan(5) + 0.1; % check unequal pans
for i=1:npan
  L = tpan(i+1)-tpan(i);       % param-length of this pan (later must be equal)
  sc = L/2;             % scale fac
  pan(i).w = sc*w(:);       % pan is struct array
  pan(i).c = (tpan(i+1)+tpan(i))/2;       % center param of pan
  pan(i).t = pan(i).c + sc*z(:);
  pan(i).sc = sc;
end
s.t = vertcat(pan.t); s.w = vertcat(pan.w);     % concat all pan nodes
inpan = @(t,i) mod(t,2*pi)>=tpan(i) & mod(t,2*pi)<tpan(i+1);  % is param in pan?
if verb>1, figure(1); clf; plot([s.t s.t]',[0*s.w s.w]','b.-');
  xlabel('t'); title('panel quadr');
  sum(s.w)-2*pi
end

paux = 20;    % nodes either side (twice this = # aux nodes per target)
[az aw] = gauss(paux);  % discont kernel can use plain G-L as aux scheme

% toy kernel function... (vectorizes to rect array given eg s=row and c=col)
%ker = @(t,s) 0.5./sin(abs(s-t)/2);   % inv dist btw pts s,t on unit circ
% matrix will be panelwise circulant
ker = @(t,s) cos(t-s);               % simple smooth kernel

N=npan*p; A = nan(N,N); % fill A...
for ip=1:npan          % target pans (block rows)
  ipl = mod(ip-2,npan)+1; ipr = mod(ip,npan)+1;     % left & right pan inds
  jfar = find((1:npan ~= ipl) & (1:npan ~= ip) & (1:npan ~= ipr));  % far pans
  jjfar = ones(p,1)*(p*(jfar-1)) + (1:p)';         % convert to far node inds
  jjfar = jjfar(:);
  ii = (1:p)+(ip-1)*p;   % block indices in A for target pan
  A(ii,jjfar) = ker(s.t(ii),s.t(jjfar)') .* s.w(jjfar)';   % broadcast to rows
  % now do 3-panel special quad for discont ker, separately for each targ in pan
  ips = {ipl, ip, ipr};        % 3 special pan inds: left, self, right
  t0 = tpan(ipl); if ip==1, t0=t0-2*pi; end      % special bottom end, unwrap
  t1 = tpan(ipr+1); if ipr==1, t1=t1+2*pi; end    % special top end, unwrap
  for k=1:p           % nodes in ip'th targ pan
    i = (ip-1)*p + k;        % targ node
    t = s.t(i);              % targ param
    sc0 = (t-t0)/2; sc1 = (t1-t)/2;    % half-param-lengths to special interval
    tq = [(t+t0)/2 + sc0*az', (t1+t)/2 + sc1*az'];   % aux quad node params, row
    wq = [sc0*aw, sc1*aw];                       % aux quad wei, row
    if verb & ip==7 & k==5, figure(2); clf; plot(tq,0*tq,'.'); hold on;
      plot(t,1,'*'); plot(tpan, 2+0*tpan, '+'); end   % debug
    for sp=1:3               % loop over special sauce panels
      jp = ips{sp};          % src pan ind
      qq = inpan(tq,jp);     % aux inds in this src pan
      tj = pan(jp).t - 2*pi*(ip==1 & sp==1) + 2*pi*(ip==npan & sp==3);  % wrap
      tc = mean(tj);         % where center for the interpolation (stable)
      I = interpmat_1d(tq(qq)-tc, tj-tc);     % aux vals from src node vals
      A(i, (jp-1)*p + (1:p)) = (ker(t,tq(qq)) .* wq(qq)) * I;    % aux wei
    end
  end
end
if verb, figure(3); clf; imagesc(A); axis equal; colorbar; end

u = sin(s.t+0.7); norm(A*u-pi*u)  % try eigfunc, eigval = pi
u = 1 + 0*s.t; norm(A*u)  % try eigfunc, eigval = 0
u = sin(2*s.t); norm(A*u)  % try eigfunc, eigval = 0
                           %figure(4); plot(s.t, [u, A*u], '.-');
