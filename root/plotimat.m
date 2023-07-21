function plotimat(I,V,stm,bsl,seg)

%Plot Controls
sr = 0.1; %sis 0.1s. Change if otherwise.
smpV = 50; % n. of voltage points for smoothing
smpT = 50; % n. of time points for smoothing

t = 0:sr:size(I,2)*sr-sr;
t_bsl = [closest(bsl(1),t) closest(bsl(2),t)]; % baseline start/end indices
i_bsl = I(:,t_bsl(1):t_bsl(2)); % all current values during baseline
mi_bsl = mean(i_bsl,2,'omitnan'); % take mean I curve during baseline
for i = 1:size(I,2)
    I(:,i) = I(:,i) - mi_bsl;
end

if ~isfinite(seg(1))
    seg(1) = 0;
end
if ~isfinite(seg(2))
    seg(2) = t(end);
end

inds = (seg)/sr + 1;
nt = t(inds(1):inds(2));
nI = I(:,inds(1):inds(2));

Vi = [1:length(V)]; % voltage indices
Vt = [1 200 400 600 800 1000]; % voltage axis tick mark labels

col = makeColorMap([1 0 0],[0 0 1],size(nI,2));
pi = reshape(nI,1,numel(nI));
mx_edge = max(pi);
mn_edge = min(pi);
up_edge = prctile(pi,99.9);
lo_edge = prctile(pi,0.1);

figure
colormap hsv
imagesc(nt,Vi,nI,[lo_edge up_edge]);
colorbar
Vtl = num2cell(ndecp(V(Vt),2));
set(gca,'YDir','normal','YTick',[0 200 400 600 800 1000],'YTickLabel',Vtl)
xlabel ('time (s)')
ylabel ('voltage (V)')
hold on

if ~isempty(stm)
     mxy = max(get(gca,'YLim')); mny = min(get(gca,'YLim'));
     yp = [mny mny mxy mxy];
     for i = 1:length(stm)
         xp = [stm(i),stm(i)+20,stm(i)+20,stm(i)];
         patch(xp,yp,[1 0 0],'FaceAlpha',0.3)
     end
end