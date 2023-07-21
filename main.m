clearvars -except datatable
%% Add function paths
p = pwd;
addpath(strcat(p,'\root'))
addpath(strcat(p,'\dependents'))
%% Compile and pre-process FSCV
fdir = uigetdir('F:\FSCV');
inputdir = strsplit(fdir,'/');
fname = genvarname(inputdir{1,end});
fcsvdata = extract_fcsv(fdir,fname);
V = fcsvdata.voltage{1,1};
select = multiselect(fcsvdata.name)';
%% Trim to stim triggers and baseline subtract
total_t = 1800;
stmstart = 101;
stmend = 200;

imat = cell(1,length(fcsvdata.name));
for i = 1:length(fcsvdata.name)
    temp = fcsvdata.current{i};
    for j = 1:size(temp,3)
        imat{i} = [imat{i},temp(:,:,j)];
    end
end

fscv = zeros(1000,total_t,length(select));
n = 1;
for i = select
    cmap = imat{i};
    if size(imat{i},2) < total_t
        cmap(:,end:total_t) = ones(size(cmap,1),total_t - size(cmap,2) + 1)*NaN;
    else
        cmap = cmap(:,1:total_t);
    end
    fscv(:,:,n) = cmap;
    bsl = mean(fscv(:,10:40,n),2,'omitnan');
    for j = 1:size(fscv,2)
        fscv(:,j,n) = fscv(:,j,n) - bsl;
    end
    n = n + 1;
end
%% Quantify FSCV data and output as table
for i = 1:size(fscv,3)
    [l(i),h(i),m(i),v(i)] = bounddetect_v2(fscv(:,:,i),stmstart,stmend);
end
midx = find(m == max(m));
wv = mean(fscv(:,l(midx):h(midx),midx),2,'omitnan');
load('wvfilt.mat','wvD')
fwv = filtfilt(wvD,wv);
[bands,Eo] = dynamicband(fwv);
for i = 1:size(fscv,3)
    [maxsig(i),tmax(i),Qo(i),io(i,:)] = quantifyFSCV(fscv(:,:,i),bands,[l(i),h(i)],Eo);
end
%% Plot Results
figure
pi = reshape(fscv,1,numel(fscv));
up_edge = prctile(pi,99.9);
lo_edge = prctile(pi,0.1);

q = 2;
p = round(size(fscv,3)/q);
while round(size(fscv,3)/q) ~= ceil(size(fscv,3)/q)
    q = q + 1;
    p = round(size(fscv,3)/q);
end
n = 1;
for i = 1:p
    for j = 1:q
        while n <= size(fscv,3)
            subplot(q,p,n)
            plotfscv(fscv(:,:,n),V,[stmstart stmend],[lo_edge,up_edge],[l(n),h(n)],bands,Eo); %set(gca,'Visible','off')
            xlabel(num2str(n))
            title(fcsvdata.name(select(n)))
            n = n + 1;
        end
    end
end