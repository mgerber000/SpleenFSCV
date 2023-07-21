idx = 5;
stimparams= [0,0,0,0]';
type = 'Sal inf';

bounds = [l(idx),h(idx)]';
ioM = maxsig(idx);
ioT = tmax(idx);

wv = mean(fscv(:,l(idx):h(idx),idx),2,'omitnan');
load('wvfilt.mat','wvD')
fwv = filtfilt(wvD,wv);

if ~exist('datatable','var')
    varNames = ["FileDir","StimName","StimParams",'StimBounds','Type','Bounds','fwv','Eo','io','Qo','io_max','t_max'];
    datatable = cell2table({fdir,[num2str(idx),fcsvdata.name{1}],stimparams,[stmstart;stmend],type,bounds,fwv,Eo,io(idx,:)',Qo(idx),ioM,ioT});
    datatable.Properties.VariableNames = varNames;
else
    T = cell2table({fdir,[num2str(idx),fcsvdata.name{1}],stimparams,[stmstart;stmend],type,bounds,fwv,Eo,io(idx,:)',Qo(idx),ioM,ioT});
    varNames = ["FileDir","StimName","StimParams",'StimBounds','Type','Bounds','fwv','Eo','io','Qo','io_max','t_max'];
    T.Properties.VariableNames = varNames;
    datatable = [datatable;T];
end

save('data_table.mat','datatable');
%fcsvdata.name(select(idx))
%[num2str(idx),fcsvdata.name{1}]