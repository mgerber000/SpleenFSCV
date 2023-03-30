function [lowbound,highbound,maxsig] = bounddetect(reg,k)

cmap = reg.fcsvdata.cmap{k,1};
load norms.mat

for i = 1:size(cmap,2)
    rho(i) = sum(cmap(:,i).*normsig);
end

% rho = filtfilt(D,rho);
poutliers = find(rho > prctile(rho,99.5));
noutliers = find(rho < prctile(rho,0.5));
rho(poutliers) = NaN;
rho(noutliers) = NaN;
rho = smooth(rho,100);

stm = reg.fcsvdata.stm_fcsv{k,1};

stmstart = find(diff(stm < 0),1)+1;
[maxsig,maxpos] = findpeaks(rho);
[~,minpos] = findpeaks(-rho);

maxsig = maxsig(maxpos < stmstart + floor(size(cmap,2)/2));
maxpos = maxpos(maxpos < stmstart + floor(size(cmap,2)/2));
% maxsig = maxsig(maxpos > stmstart);
% maxpos = maxpos(maxpos > stmstart);
maxpos = maxpos(maxsig == max(maxsig));
maxsig = maxsig(maxsig == max(maxsig));

% lowbound = find((diff(rho(1:maxpos-1) > 0) == 1))+1;
% if ~isempty(lowbound)
%     lowbound = lowbound(1,end);
%     if lowbound < stmstart
%        lowbound = stmstart;
%     end
% else
%     lowbound = find(diff(minpos > maxpos) == 1);
%     if ~isempty(lowbound)
%         lowbound = minpos(lowbound(1,1));
%     end
    lowbound = stmstart;
% end
highbound = find((diff(rho(maxpos+1:end) > 0) == -1));
if ~isempty(highbound)
    highbound = highbound(1,1)+1+maxpos;
else
%     highbound = find(diff(minpos > maxpos) == 1);
%     if ~isempty(highbound)
%         highbound = minpos(highbound(1,1)+1);
%     end
    highbound = size(cmap,2);
end
if highbound < lowbound
    highbound = lowbound + 1;
end