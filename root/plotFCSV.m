function [I,V,sr,bsl,stm,smpT,segments,plottype,averaging,trace] = plotFCSV

%Plot Controls
sr = 0.1; %sis 0.1s. Change if otherwise.
bsl = [1 1 4]; % define baseline [segment start-sec end-sec]
stm = [1 6 14]; % define stim train [segment start-sec end-sec]
smpV = 50; % n. of voltage points for smoothing
smpT = 50; % n. of time points for smoothing

segments = [1 400]; %List # of segments to plot. Leave empty for all. Example: For segments 1-4, enter [1 4]. If asked for more segments then available, plots all.
plottype = 'waveform'; %Plot waveforms, colormap or both. Enter 'waveform' for just waveform. Enter 'colormap' for just colormap. Enter 'both' for both. Default = both.
averaging = 1; %Average per n waveforms. Leave blank for no averaging.
trace = []; %Enter value for trace of current vs time at specified voltage. Leave empty for no trace. Value outside of voltage limit will be ignored.

%Get file name
[fname, fdir] = uigetfile('*.db3*');

%Get data file directory
fname(strfind(fname,'.db3'):end) = [];
tsv = dir([fdir fname '*.tsv']);

%Order list of files by date created.
for i = 1:size(tsv,1)
    list(i,:) = strsplit(tsv(i,1).date,{'-',' ',':'});
    list(i,2) = {monthconvert(char(list(i,2)))};
end
order = dateorder(list);

%Extract data
current{size(tsv,1),1} = [];
for i=1:length(tsv)
    copyfile ([fdir tsv(order(i)).name], [tsv(order(i)).name(1:end-4) '.txt']) % create text file
    data = readtable ([tsv(order(i)).name(1:end-4) '.txt']); % extract data
    delete ([tsv(order(i)).name(1:end-4) '.txt'])
    
    ii=table2array(data(:,2:end)); % current values
    for j=1:size(ii,2)
        ii(:,j) = smooth(ii(:,j),smpV);
    end
    current{i,1} = ii;
end

I = zeros(size(current{1,1},1),size(current{1,1},2),size(current,1));
for i = 1:size(current,1)
    I(:,:,i) = current{i,1};
end

V = table2array(data(:,1)); % voltage values (same for every file)

%Average baseline and normalize all data
t = 0:sr:size(I,2)*sr-sr;
t_bsl = [closest(bsl(2),t) closest(bsl(3),t)]; % baseline start/end indices
i_bsl = I(:,t_bsl(1):t_bsl(2),bsl(1)); % all current values during baseline
% for j=1:size(i_bsl,2)
%     i_bsl(:,j)=smooth(i_bsl(:,j),smpV);
% end

mi_bsl = mean(i_bsl,2); % take mean I curve during baseline
for i = 1:size(I,3)
    for j = 1:size(I,2)
        I(:,j,i) = I(:,j,i) - mi_bsl;
    end
end
handles = struct();
%Plot based on set parameters
createPlot(handles,I,V,sr,bsl,stm,smpT,segments,plottype,averaging,[])