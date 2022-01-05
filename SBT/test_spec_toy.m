% Examine spectrum of toy: scalar SBT operator (zz-part) on the unit circle
% The operator is (Ku)(t) := int_0^2pi (u(s)-u(t))/d(t,s) ds
% where d(t,s) = 2 sin |s-t|/2 is the distance btw angles s,t on unit circle.
%
% Barnett 12/23/21
clear

%ker = @(t,s) cos(mod(s-t,2*pi)/2);   % warm-up discont ker (half-cycle of cos)
% (this just tested convernence of point spectrum, weak decaying eigvals...)

% toy kernel:
ker = @(t,s) 0.5./sin(abs(s-t)/2);   % inv dist btw pts s,t on unit circ

p = 16;
npans = 40:20:100;   % convergence study
figure(1); clf; subplot(2,1,1);
for i=1:numel(npans); npan=npans(i)
  tpan = 2*pi*(0:npan)/npan;    % pan param breakpoints (first=0, last=2pi)
  pan = setup_pans(tpan,p);
  A = nyst_diagdiscont_sca(pan,tpan,ker);
  A = A - diag(sum(A,2));       % Nystrom for u(s) becoming u(s)-u(t) in K apply
  N = size(A,1);
  lam{i} = sort(real(eig(A)),'descend');
  plot(1:N, lam{i}, '.'); hold on; axis tight; drawnow; 
end
legnum(npans,3,'npan=');
xlabel('j'); ylabel('Re \lambda_j'); title('spec(A), convergence');
subplot(2,1,2);
Nt=numel(lam{end})/4; semilogy(1:Nt, abs(lam{end}(1:Nt)-lam{end-1}(1:Nt)),'+');
xlabel('j'); ylabel('est error in \lambda_j');
print -dpng figs/spec_toy_conv.png

% analytic study:
figure(2); clf;
nmax = 300;             % max Fourier mode to predict spectrum to
lampred = [0 kron(-4*cumsum(1./(1:2:2*nmax-1)), [1 1])];   % a zero then pairs of odd-harmonic sum eigvals
subplot(2,1,1); plot(1:(2*nmax+1),lampred,'.-');
hold on; semilogx(1:N, lam{end}, '.');
axis tight; xlabel('j'); ylabel('Re \lambda_j'); set(gca,'xlim',[1 2*nmax+1]);
legend('spec(A)', 'pred analytic');
title('check the spec(A) vs analytic formula');
subplot(2,1,2);
semilogy(1:(2*nmax+1), abs(lampred(:)-lam{end}(1:(2*nmax+1))),'+-');
title('differences'); axis tight;
print -dpng figs/spec_toy_form.png
% this shows the eigenvalues are precisely the analytic ones, good
