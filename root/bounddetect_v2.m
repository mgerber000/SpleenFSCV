function [lowbound,highbound,validpeaks,peakvals] = bounddetect_v2(cmap,stmstart,stmend)

load('norms.mat','normsig')

%Calculate dot product of raw voltammagrams with NE signal norm
rho = zeros(1,size(cmap,2));
for i = 1:size(cmap,2)
    rho(i) = sum(cmap(:,i).*normsig);
end

rho(rho > prctile(rho,99.5)) = NaN;
rho(rho < prctile(rho,0.5)) = NaN;
rho = smooth(rho,100);

%Find peaks that meet threshold prominence (determined from lowest injected
%NE concentration peaks)
[maxsig,maxpos] = findpeaks(rho,'MinPeakProminence',2200);
[~,minpos] = findpeaks(-rho,'MinPeakProminence',2200);

%Identify peaks and troughs
[ipoints,idx] = sort([minpos;maxpos]);
pt = [ones(length(minpos),1);ones(length(maxpos),1)*2];
pt = pt(idx);

%Find peaks that are greater than zero and occur between stim start 
%and 1 minute after stim end
validpeaks = maxpos((maxpos > stmstart) & (maxpos < stmend + 600) & maxsig > 0);

%A valid peak must be bounded by local minima or zero before and after.
%Exclude all peaks that do not satisfy this condition
if ~isempty(validpeaks)
    if find(ipoints == validpeaks(1)) - 1 < 1
        ipoints = [stmstart;ipoints];
        pt = [1;pt];
    end
    if ipoints(find(ipoints == validpeaks(1)) - 1) < stmstart
        ipoints(find(ipoints == validpeaks(1)) - 1) = stmstart;
    end
    lowbound = zeros(1,length(validpeaks))*NaN;
    highbound = zeros(1,length(validpeaks))*NaN;
    peakvals = zeros(1,length(validpeaks))*NaN;
    for i = 1:length(validpeaks)
        if length(i) ~= 1 && i == 1
            pre_t = (stmstart-1):validpeaks(i);
            post_t = validpeaks(i)+1:validpeaks(i+1);
            bound1 = pre_t(diff(rho(pre_t) > 0) == 1);
            bound2 = post_t((diff(rho(post_t) > 0) == -1));
        elseif length(i) ~= 1 && i == length(i)
            pre_t = validpeaks(i-1)+1:validpeaks(i);
            post_t = validpeaks(i)+1:length(rho);
            bound1 = pre_t(diff(rho(pre_t) > 0) == 1);
            bound2 = post_t((diff(rho(post_t) > 0) == -1));
        elseif length(i) == 1
            pre_t = (stmstart-1):validpeaks(i);
            post_t = validpeaks(i)+1:length(rho);
            bound1 = pre_t(diff(rho(pre_t) > 0) == 1);
            bound2 = post_t((diff(rho(post_t) > 0) == -1));
        else
            pre_t = validpeaks(i-1)+1:validpeaks(i);
            post_t = validpeaks(i)+1:validpeaks(i+1);
            bound1 = pre_t(diff(rho(pre_t) > 0) == 1);
            bound2 = post_t((diff(rho(post_t) > 0) == -1));
        end
        if ~isempty(bound1)
            lowbound(i) = bound1(end);
        elseif isempty(bound1) && ~isempty(ipoints(ipoints > pre_t(1) & ipoints < pre_t(end) & pt == 1))
            lpoints = ipoints(ipoints >= pre_t(1) & ipoints < pre_t(end) & pt == 1);
            lowbound(i) = lpoints(end);
        else
            lowbound(i) = NaN;
        end
        if ~isempty(bound2)
            highbound(i) = bound2(1);
        elseif isempty(bound2) && ~isempty(ipoints(ipoints > post_t(1) & ipoints < post_t(end) & pt == 1))
            hpoints = ipoints(ipoints > post_t(1) & ipoints <= post_t(end) & pt == 1);
            highbound(i) = hpoints(1);
        else
            highbound(i) = NaN;
        end
        if ~isnan(lowbound(i)) && ~isnan(highbound(i))
            peakvals(i) = maxsig(maxpos == validpeaks(i));
        else
            validpeaks(i) = NaN;
            peakvals(i) = NaN;
        end
    end
    lowbound = min(lowbound);
    highbound = max(highbound);
    if length(validpeaks) == 1
        if isnan(validpeaks)
        else
            validpeaks = validpeaks(peakvals == max(peakvals));
        end
    else
        validpeaks = validpeaks(peakvals == max(peakvals));
    end
    if length(peakvals) == 1
        if isnan(peakvals)
        else
            peakvals = peakvals(peakvals == max(peakvals));
        end
    else
        peakvals = peakvals(peakvals == max(peakvals));
    end
else
    lowbound = NaN;
    highbound = NaN;
    validpeaks = NaN;
    peakvals = NaN;
end