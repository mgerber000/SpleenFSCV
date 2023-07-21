function [bands,mxloc] = dynamicband(fwv)

% [pks,locs] = findpeaks(fwv);
% vpks = pks(find(locs > 350 & locs < 650));
% vlocs = locs(find(locs > 350 & locs < 650));
% maxpk = vlocs(find(vpks == max(vpks)));
% mxlocs = locs;
% 
% [~,locs] = findpeaks(-fwv);
% p2 = locs(diff(locs < maxpk) == -1);
% p1 = locs(find(diff(locs < maxpk) == -1) - 1);
% p3 = locs(find(diff(locs < maxpk) == -1) + 1);
% mxlocs = mxlocs(mxlocs > p1 & mxlocs < p2);
% bands{1,1} = p1:p2;
% bands{1,2} = p2:p3;

load NEID.mat

[~,locs] = findpeaks(fwv);
mxloc = locs((pm - pstd*4) < locs & (pm + pstd*4) > locs);

[~,locs] = findpeaks(-fwv);
p1 = locs(find(diff(locs < mxloc) == -1));
p2 = locs(find(diff(locs > mxloc) == 1) + 1);
p3 = locs(find(locs == p2) + 1);


if isempty(p1) || isempty(p2) || isempty(mxloc)
    bands{1,1} = NaN;
    bands{1,2} = NaN;
else
<<<<<<< HEAD
    bands = [p1,p2,p3];
=======
    bands{1,1} = p1:p2;
    bands{1,2} = p2:p3;
>>>>>>> 6bdb159a1ab162fedc9fd3f11bd7af14af5d5c3f
end

% minpeakprom = 0;
% df = 0.1;
% [~,locs] = findpeaks(-fwv);
% 
% while size(locs,1) ~= 5
%     d = (size(locs,1) - 5)*df;
%     if d < 0
%         df = df*0.1;
%         locs = slocs;
%         minpeakprom = smpp;
%     else
%         slocs = locs;
%         smpp = minpeakprom;
%         minpeakprom = minpeakprom + d;
%         [~,locs] = findpeaks(-fwv,'MinPeakProminence',minpeakprom);
%     end
% end
% 
% bands{1,1} = locs(2):locs(3);
% bands{1,2} = locs(3):locs(4);
% bands{1,3} = locs(4):locs(5);