function [handles,ci] = createPlot(handles,I,V,sr,bsl,stm,smpT,segments,plottype,n,f)
%INPUTS
%handles - figure object for pushing plot display to correct window
%I - current values in mxnxp matrix where m = voltage index, n = segment size and p = segment number
%V - voltage index
%sr - sampling rate
%bsl - baseline current values in mxn matrix where m = voltage index, n = number of waveforms
%stm - stim location in 1x2 matrix representing start and end time
%smpT - smoothing coefficient over time
%segments - range of waveforms to be plotted
%plottype - choose whether to plot waveforms, colormap or neither
%n - number of waveforms to average
%f - filter object for filtering over time axis

%OUTPUTs
%ci - mxn matrix of current values where m = voltage index, n = number of waveforms
%% Baseline correction
mi_bsl = mean(bsl,2);
for i = 1:size(I,3)
    for j = 1:size(I,2)
        I(:,j,i) = I(:,j,i) - mi_bsl;
    end
end
%% Cutting and Squeezing Data to 2D matrix
if isempty(segments)
    startind = 1; endind = size(I,3);
elseif size(segments,2) == 1
    startind = segments; endind = segments;
elseif size(segments,2) == 2
    if segments(1) < 1 || segments(2) > size(I,3)*size(I,2)
        startind = 1; endind = size(I,3)*size(I,2);
    else
        startind = segments(1); endind = segments(2);
    end
end

ci = [];
for i = 1:size(I,3)
    ci = [ci , I(:,:,i)];
end
ci = ci(:,startind:endind);

% for i = 1:size(ci,1)
%     ci(i,:) = smooth(ci(i,:),smpT)';
% end
%% Select plot type
switch plottype
    case 'none'
        waveform = false; cmap = false; plotnum = 0;
    case 'waveform'
        waveform = true; cmap = false; plotnum = 1;
    case 'colormap'
        waveform = false; cmap = true; plotnum = 1;
    case 'both'
        waveform = true; cmap = true; plotnum = 2;
    otherwise
        waveform = true; cmap = true; plotnum = 2;
end
%% Filter over time axis
if ~isempty(f)
    for j = 1:size(cmap,1)
        cmap(j,:) = filtfilt(f,cmap(j,:));
    end
end
%% Average every n increments
if isempty(n)
    n = 1;
end

%Average every n waveforms
if n > 1
    avI = zeros(size(ci,1),floor(size(ci,2)/n));
    ind = 1;
    for j = 1:n:size(ci,2)-n
        avI(:,ind) = mean(ci(:,j:j+n-1),2);
        ind = ind+1;
    end
    avI = avI;
else
    avI = ci;
end
%% Create plot if applicable
%Define indices
if size(avI,2) > 1
    col = makeColorMap([1 0 0],[0 0 1],size(avI,2));
end
Vi = [1:length(V)]; % voltage indices
Vt = [1 200 400 600 800 1000]; % voltage axis tick mark labels

%Define edges of current plots
pi = reshape(avI,1,numel(avI));
mx_edge = max(pi);
mn_edge = min(pi);
up_edge = prctile(pi,99.9);
lo_edge = prctile(pi,0.1);

% figure

if waveform == true
    for j=1:size(avI,2)
        if size(avI,2) > 1
            handles.axes1 = plot(Vi,avI(:,j),'Color',col(j,:));
        else
            handles.axes1 = plot(Vi,avI(:,j),'-b');
        end
        hold on
    end
    xlabel ('voltage index')
    ylabel ('current')
    ylim([mn_edge mx_edge])
end

if cmap == true
    colormap hsv
    t = (startind-1)*sr:sr:(endind-1)*sr;
%     axes(handles.axes1);
    imagesc(t,Vi,ci,[lo_edge up_edge]);
    colorbar
    Vtl = num2cell(ndecp(V(Vt),2));
    set(gca,'YDir','normal','YTick',[0 200 400 600 800 1000],'YTickLabel',Vtl)
    xlabel ('time (s)')
    ylabel ('voltage (V)')

    mxy = max(get(gca,'YLim')); mny = min(get(gca,'YLim'));
%     xp_bsl = [bsl(2) bsl(3) bsl(3) bsl(2)];
    xp_stm = [(stm(1)-1)*sr (stm(2)-1)*sr (stm(2)-1)*sr (stm(1)-1)*sr];
    yp = [mny mny mxy mxy];
%     patch(xp_bsl,yp,[0.5 0.5 0.5],'FaceAlpha',0.2)
    patch(xp_stm,yp,[1 0 0],'FaceAlpha',0.3)
end