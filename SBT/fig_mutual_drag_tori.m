% read Dhairya data comparing SBT to convergent BIE for mutual grad of two
% tori.
% Barnett 3/22/23

clear

if 0   % ------------------------ VERSION 1
f=fopen('output/mutual_drag_tori_Dhairya_err.dat');   % paste from email
a = textscan(f,'%f %f %f %f %f','commentstyle','#');
eps=a{1}; delta=a{2};
FSBT=a{3}; FBIE=a{4};      % z-cmpts of F1 vector mutual drag
fclose(f);

relerr = abs((FSBT-FBIE)./FBIE);
err = abs(FSBT-FBIE);
ratio = delta./eps;
%model = 10 * eps.*log(1./eps) ./ ratio.^1.5;      % error model
model = 4*log(1./eps) .* eps.^2.5 ./ delta.^1.5;      % error model
                                                      % (powers ok, log not?)
figure;     % ............. plot relerr vs ratio, grouped by eps
eu = unique(eps)';  % has to be row vec for loop to work
rng(0);
for e=eu, j=find(eps==e);
  if ~isempty(j)
    c = 0.8*rand(1,3);
    loglog(ratio(j),err(j),'+-','color',c); hold on;
    plot(ratio(j),model(j),'o--','color',c); hold on;
    text(mean(ratio(j)), mean(err(j)), sprintf('\\epsilon = %g',eps(j(1))))
  end
end
xlabel('ratio \delta/\epsilon');
ylabel('err');
title('error in SBT F1_z (vs BIE ground truth)');
%legnum(eu,[],'\epsilon = ');
print -dpng figs/mutual_tori_err_draft.png
end  % ---------------------------

if 0    % ------------------- VERSION 2
f=fopen('output/mutual_drag_tori_twisted_grid.dat');  % no errs
a = textscan(f,'%f %f %f %f %f','commentstyle','#');
eps=a{1}; delta=a{2};
F1x=a{3}; F1y=a{4}; F1z=a{5};     % cmpts of F1 vector mutual drag
fclose(f);
epss = unique(eps);
dels = unique(delta);
ne = numel(epss); nd = numel(dels);
le = log10(epss);
figure; imagesc(le,log10(dels),reshape(-F1z,[nd ne])); colorbar;
xlabel('log_{10}\epsilon'); ylabel('log_{10}\delta'); axis xy tight;
v=axis; hold on;
plot(le,le,'k--'); text(-2,-2,"\delta=\epsilon");
plot(le,1+le,'k--'); text(-2,-1,"\delta=10\epsilon");
plot(le,2+le,'k--'); text(-2,0,"\delta=100\epsilon");
axis(v);
title('mutual drag twisted tori (-F1_z), SBT, grid samples');
print -dpng figs/mutual_tori_twisted_grid_F1z_SBT.png
end

% ----------------------- VERSION 3
f=fopen('output/mutual_drag_tori_twisted_grid_err.dat');  % Dhairya errs vs BIE
a = textscan(f,'%f %f %f %f %f %f %f %f %f','commentstyle','#');
eps=a{1}; delta=a{2};
F1x=a{3}; F1y=a{4}; F1z=a{5};     % cmpts of F1 vector mutual drag, SBT
B1x=a{6}; B1y=a{7}; B1z=a{8};     % cmpts of F1 vector mutual drag, BIE
fclose(f);
FSBT = [F1x,F1y,F1z]; FBIE = [B1x,B1y,B1z];
Fmag = sqrt(sum(FBIE.^2,2));       % accurate (BIE) ||F||
Frelerr = sqrt(sum((FBIE-FSBT).^2,2)) ./ Fmag;  % rel errs SBT vs BIE
epss = unique(eps);
dels = unique(delta);
ne = numel(epss); nd = numel(dels);
le = log10(epss);
figure; subplot(1,2,1);
imagesc(le,log10(dels),reshape(Fmag,[nd ne]));
caxis([0 25]); colorbar;
xlabel('log_{10}\epsilon'); ylabel('log_{10}\delta'); axis xy tight;
v=axis; hold on;
plot(le,le,'k--'); text(-2,-2,"\delta=\epsilon");
plot(le,1+le,'k--'); text(-2,-1,"\delta=10\epsilon");
plot(le,2+le,'k--'); text(-2,0,"\delta=100\epsilon");
axis(v);
title('mutual drag vector magnitude ||F||, twisted tori');
subplot(1,2,2);
%[C,h] = contourf(le,log10(dels),reshape(log10(Frelerr),[nd ne]), [-7:-3, -2:0.5:0]);
%clabel(C,h); caxis([-7 -0.24]);
imagesc(le,log10(dels),reshape(Frelerr,[nd ne])); set(gca,'colorscale','log'); caxis([1e-7 0.6]); % fails for contourf
colorbar;
xlabel('log_{10}\epsilon'); ylabel('log_{10}\delta'); axis xy tight;
hold on;
plot(le,le,'k--'); text(-2,-2,"\delta=\epsilon");
plot(le,1+le,'k--'); text(-2,-1,"\delta=10\epsilon");
plot(le,2+le,'k--'); text(-2,0,"\delta=100\epsilon");
axis(v);
title('mutual drag vector error ||F_{SBT}-F||/||F|| twisted tori');
%title('mutual drag vector error log_{10} ||F_{SBT}-F||/||F|| twisted tori');
set(gcf,'paperposition',[0 0 13 5])
%print -dpng figs/mutual_tori_twisted_grid_err.png
print -dpng figs/mutual_tori_twisted_grid_err2.png
