% read Dhairya data comparing SBT to convergent BIE for mutual grad of two
% tori.
% Barnett 3/22/23

close
f=fopen('output/mutual_drag_tori_Dhairya_err.dat');   % paste from email
a = textscan(f,'%f %f %f %f %f','commentstyle','#');
eps=a{1}; delta=a{2};
FSBT=a{3}; FBIE=a{4};      % z-cmpts of F1 vector mutual drag
fclose(f);

relerr = abs((FSBT-FBIE)./FBIE);
ratio = delta./eps;
figure;     % ............. plot relerr vs ratio, grouped by eps
eu = unique(eps)';  % has to be row vec for loop to work
for e=eu, j=find(eps==e);
  if ~isempty(j)
    loglog(ratio(j),relerr(j),'+-'); hold on;
  end
end
xlabel('ratio \delta/\epsilon');
ylabel('rel err');
legnum(eu,[],'\epsilon = ');



