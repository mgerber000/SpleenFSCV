function plotfscv(I,V,stm,cl,bounds,bands,Eo)

%Plot Controls
sr = 0.1; %sis 0.1s. Change if otherwise.
smpV = 50; % n. of voltage points for smoothing
smpT = 50; % n. of time points for smoothing

t = 0:sr:(size(I,2)*sr - sr);

Vi = 1:length(V); % voltage indices
Vt = 1000 - [50 473 898]; % voltage axis tick mark labels
V = flip(V,1);
% pi = reshape(I,1,numel(I));
% up_edge = prctile(pi,99.9);
% lo_edge = prctile(pi,0.1);

% figure
colormap hsv
imagesc(t,1:1000,flip(I,1),cl);
Vtl = num2cell(ndecp(V(Vt),2));
set(gca,'YTick',[50 473 898],'YTickLabel',Vtl)
% xlabel ('time (s)')
% ylabel ('voltage (V)')
hold on

if ~isempty(stm)
     mxy = max(get(gca,'YLim')); mny = min(get(gca,'YLim'));
     yp = [mny mny mxy mxy];
     stm = stm*sr - sr;
     xp = [stm(1),stm(2),stm(2),stm(1)];
     patch(xp,yp,[1 0 0],'FaceAlpha',0.3)
end

if ~isnan(bounds(1)) && ~isnan(bounds(2))
    bounds = bounds*sr-sr;
    scatter(bounds(1):bounds(2),(1000 - Eo)*ones(1,length(bounds(1):bounds(2))),3,'k','filled')
    scatter(bounds(1):bounds(2),(1000 - bands(1))*ones(1,length(bounds(1):bounds(2))),3,'k','filled')
    scatter(bounds(1):bounds(2),(1000 - bands(2))*ones(1,length(bounds(1):bounds(2))),3,'k','filled')

    scatter(bounds(1)*ones(1,length((1000-bands(2)):(1000-bands(1)))),(1000-bands(2)):(1000-bands(1)),3,'k','filled')
    scatter(bounds(2)*ones(1,length((1000-bands(2)):(1000-bands(1)))),(1000-bands(2)):(1000-bands(1)),3,'k','filled')
else
    scatter(bounds(1)*ones(1,length((1000-bands(2)):(1000-bands(1)))),(1000-bands(2)):(1000-bands(1)),3,'k','filled')
end
