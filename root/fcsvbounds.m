function [pts,res] = fcsvbounds(locs)

for i = 1:size(locs,1)
   tempm(i,1) = mean(locs{i,1}) ;
   temps(i,1) = std(locs{i,1});
   tempr(i,1) = rms(locs{i,1});
end

s = 20;
b = 1/s*ones(1,s);
a = 1;

%ftempm = filtfilt(b,a,tempm);
ftemps = filtfilt(b,a,temps);
ftempr = filtfilt(b,a,tempr);

[pts,res] = findchangepts(temps,'MaxNumChanges',1,'Statistic','rms');
[fpts,fres] = findchangepts(ftemps,'MaxNumChanges',1,'Statistic','rms');
[pts2,res2] = findchangepts(tempr,'MaxNumChanges',1,'Statistic','rms');
[fpts2,fres2] = findchangepts(ftempr,'MaxNumChanges',1,'Statistic','rms');

% figure
% %plot(tempm)
% hold on
% plot(temps)
% plot(tempr)
% %plot(ftempm)
% plot(ftemps)
% plot(ftempr)
% scatter(pts,temps(pts))
% scatter(pts2,tempr(pts2))
% scatter(fpts,ftemps(fpts))
% scatter(fpts2,ftempr(fpts2))