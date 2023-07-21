function [I,V] = readFCSV(fname,fdir,smpV)
    %Get data file directory
    fname(strfind(fname,'.db3'):end) = [];
    tsv = dir([fdir '/' fname '*.tsv']);

    %Order list of files by date created.
    for i = 1:size(tsv,1)
        list(i,:) = strsplit(tsv(i,1).date,{'-',' ',':'});
        list(i,2) = {monthconvert(char(list(i,2)))};
    end
    order = dateorder(list);

    %Extract data
    current{size(tsv,1),1} = [];
    for i=1:length(tsv)
        copyfile ([fdir '/' tsv(order(i)).name], [tsv(order(i)).name(1:end-4) '.txt']) % create text file
        data = readtable ([tsv(order(i)).name(1:end-4) '.txt']); % extract data
        delete ([tsv(order(i)).name(1:end-4) '.txt'])

        ii=table2array(data(:,2:end)); % current values
        for j=1:size(ii,2)
            ii(:,j) = smooth(ii(:,j),smpV);
        end
        current{i,1} = ii;
    end
    I = zeros(size(current{1,1},1),size(current{1,1},2)+5,size(current,1));
    for i = 1:size(current,1)
        I(:,:,i) = [current{i,1},zeros(size(current{i,1},1),5)*NaN];
    end

    V = table2array(data(:,1)); % voltage values (same for every file)
end