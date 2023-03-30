function [dbeats,sigband,maxsig,tmax,tlowbound,thighbound,NEval,dtrace] = quantifyFSCV(reg,k,bands,timebounds,mxband)

cmap = reg.fcsvdata.cmap{k,1};
tmax = timebounds(find(cmap(mxband,timebounds) == max(cmap(mxband,timebounds))));
maxsig = cmap(mxband,tmax);

trace = cmap(mxband,timebounds);
poutliers = find(trace > prctile(trace,99.5));
noutliers = find(trace < prctile(trace,0.5));
ftrace = trace;
ftrace(poutliers) = NaN;
ftrace(noutliers) = NaN;
strace = smooth(ftrace,100);
maxi = max(strace);
tempmax = find(strace == max(strace));
tempend = find(~(diff(strace > mean(strace)*0.05) == 0));
validends = tempend(find(tempend > tempmax));
if isempty(validends)
    validends = size(trace,2);
end
thighbound = timebounds(validends(1,1));
tlowbound = timebounds(1,1);
timebounds = tlowbound:thighbound;
% trace = cmap(mxband,tlowbound:thighbound);
% NEval = trapz((tlowbound:thighbound)/10,trace);
% NEval = trapz((tlowbound:thighbound)/10,trace)/((thighbound-tlowbound)/10);


% [lowbound,highbound,~] = bounddetect(reg,k);
% timebounds = lowbound:highbound;

% for j = 1:size(timebounds,2)
%     wv(:,j) = cmap(:,timebounds(1,j))/max(cmap(:,timebounds(1,j)));
% end
% wv = mean(cmap(:,timebounds),2);
% load wvfilt.mat
% fwv = filtfilt(wvD,wv);
% [locs,~] = dynamicband(fwv);
fs = 1/10000;
% for j = 1
band = bands{1,1};
for i = 1:size(cmap,2)
    point1 = cmap(band(1,1),i);
    point2 = cmap(band(1,end),i);
    mxi = cmap(mxband,i);
%     m = (point2 - point1)/(band(1,end) - band(1,1));
%     b = point2 - m*band(1,end);
%     y = m * band + b;
%     dtrace(i) = fscv_intersect(point1,point2,band,mxband,mxi);
    dtrace(i) = fscv_peak2trough(point1,mxi) ;
end
% sigband(j) = sum(area);
sigband = NaN;
NEval = trapz(timebounds*0.1,dtrace(1,timebounds));
% end

% if reg.usephysio
%     hr = reg.fcsvdata.hr{k,1}(1:60000);
%     t = reg.fcsvdata.stm_ecg{k,1}(1:60000);
%     stmstart = find(diff(t < 0),1);
%     baselinehr = mean(hr(1:stmstart));
%     dbeats = trapz(t(stmstart+1:end),baselinehr - hr(stmstart+1:end))/60;
% else
%     dbeats = NaN;
% end

% reg.usephysio = true;
if reg.usephysio
    hr = reg.fcsvdata.hr{k,1}(1:60000);
    t = reg.fcsvdata.stm_ecg{k,1}(1:60000);
    stmstart = find(diff(t < 0),1);
    baselinehr = mean(hr(1:stmstart));
    poststmhr = mean(hr(stmstart+1:end));
    dbeats = (baselinehr - poststmhr)/baselinehr;
else
    dbeats = NaN;
end
