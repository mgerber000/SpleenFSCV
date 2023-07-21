%----------------------------------------------------------------------
function [ecglocs2,hr] = hrcalc2(x1,y1,y2,evestm,fs,thresh,f,plotidx)
% Set parameters
mpp = 400;
% thresh = 0.005;
% f = 1.8;
upperfilt = 2500;
lowerfilt = 50;
maxhr = 900;
mindist = 1/(maxhr/60);

%Filter ecg
dy1 = filter_data1(y1,fs,'butter',1,lowerfilt,upperfilt);
dy2 = filter_data1(y2,fs,'butter',1,lowerfilt,upperfilt);
ecg = dy1 - dy2;

%Find position of stim artifact
[~,art] = model_artifact(x1,ecg,evestm,fs);
[~,locs,~,p] = findpeaks(art);
artloc = (locs(p == max(p))-1)/fs;

%Find all potential peaks meeting minimum mpp threshold
[pks,locs,w,~] = findpeaks(ecg,x1,'MinPeakProminence',mpp); %'MinPeakDistance',mindist

%Eliminate peaks too close to stim artifacts
d = locs - (evestm+artloc);
d = abs(d);
mind = min(d,[],1,'omitnan');
inds = find(mind > thresh);
ecglocs = locs(inds);
ecgpks = pks(inds);
ecgw = w(inds);

%Apply minimum distance algorithm for ecg detection outside of stimulus window
idx = round((ecglocs - x1(1,1))*fs + 1);
idx = findPeaksSeparatedByMoreThanMinPeakDistance(dy1,x1,idx,mindist);
ecglocs2 = ecglocs(idx);
ecgw2 = ecgw(idx);

%Remove all peaks greater than the average widith of prestim peaks
meanecgw = mean(ecgw2(ecglocs2 < evestm(1,1)));
stdecgw = std(ecgw2(ecglocs2 < evestm(1,1)))*10;
idx = find((ecgw < meanecgw*f));
ecglocs = ecglocs(idx);
ecgpks = ecgpks(idx);

%Apply minimum distance algorithm on newly filtered peaks
idx = round((ecglocs - x1(1,1))*fs + 1);
idx = findPeaksSeparatedByMoreThanMinPeakDistance(dy1,x1,idx,mindist);
ecglocs2 = ecglocs(idx);
ecgpks2 = ecgpks(idx);

%Calculate hr
ibi = diff(ecglocs2);
hr = 60./[0 ibi];

%Determine stimstart/end
prestimlocs = find(ecglocs2 < evestm(1,1));
stimstart = ecglocs2(prestimlocs(1,end));
poststimlocs = find(ecglocs2 > evestm(end,1));
stimend = ecglocs2(poststimlocs(1,1));

% plot if called for
ax1 = subplot(2,9,plotidx(1));
hold off
plot(x1,ecg)
hold on
scatter(ecglocs2,ecgpks2,'o','g')
hold on
scatter(evestm+artloc,ones(size(evestm,1),1))
hold on
scatter(locs,pks,'*','r')
hold on
miny = ax1.YLim(1,1);
maxy = ax1.YLim(1,2);
patch([stimstart,stimstart,stimend,stimend],[miny,maxy,maxy,miny],[0.93,0.69,0.13],'EdgeColor','none','FaceAlpha',0.3)
stimparams = get_currentstimparams;
title(num2str(stimparams(plotidx(1))))

ax2 = subplot(2,9,plotidx(2));
plot(ecglocs2,smooth(hr,100))
ylim([200,800])
miny = ax2.YLim(1,1);
maxy = ax2.YLim(1,2);
patch([stimstart,stimstart,stimend,stimend],[miny,maxy,maxy,miny],[0.93,0.69,0.13],'EdgeColor','none','FaceAlpha',0.3)

linkaxes([ax1,ax2],'x');