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

figure(2); clf; semilogx(1:N, lam{i}, '.'); hold on; axis tight;
xlabel('j'); ylabel('Re \lambda_j'); 
plot(1:N,-2*cumsum(1./(1:N)),'--');    % harmonic sum, guessed prefactor 2
legend('spec(A)', '-2 * harmonicsum');
title('fit the form of spec(A) negative growth');
print -dpng figs/spec_toy_form.png
