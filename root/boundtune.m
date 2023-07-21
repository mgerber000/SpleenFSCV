function [lowbound,highbound,maxpos] = boundtune(reg,k,lowbound,highbound,mxloc)

cmap = reg.fcsvdata.cmap{k,1};
timebound = lowbound:highbound;
sig = cmap(mxloc,timebound);
sigmax = find(sig == max(sig));

t1 = 1:timebound(sigmax)-1;
t2 = timebound(sigmax):size(cmap,2);
s1 = cmap(mxloc,t1);
fs1 = smooth(s1,10);
s2 = cmap(mxloc,t2);
fs2 = smooth(s2,10);
lowbound = find(diff(fs1 < max(sig)*0.1) == -1);
highbound = t2(find(diff(fs2 < max(sig)*0.1) == 1));

% lowbound = closest(max(sig)*0.1,cmap(mxloc,t1));
% highbound = t2(closest(max(sig)*0.1,cmap(mxloc,t2)));

if ~isempty(lowbound)
    lowbound = lowbound(end,1);
else
    lowbound = timebound(1,1);
end

if ~isempty(highbound)
    highbound = highbound(1,1);
else
    highbound = size(cmap,2);
end

maxpos = timebound(sigmax);