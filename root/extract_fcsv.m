function fcsv = extract_fcsv(fcsvdir,fname)
    %Define name of fcsvmat file
    fcsvmatfile = [fcsvdir '\' 'fcsv_' fname '.mat'];
    
    %Check if file already exists
    if exist(fcsvmatfile, 'file')
        disp('loading fcsv data...')
        load(fcsvmatfile,'fcsv') %If it does, load it.
    else
        disp('extracting fcsv data...')
        %If not, extract all data available.
        D = dir(fcsvdir);

        %Find and organize db3 files by date created
        names = {D.name}';
        db3ind = find(contains(names,'.db3'));
        for i = 1:size(db3ind,1)
            list(i,:) = strsplit(D(db3ind(i)).date,{'-',' ',':'});
            list(i,2) = {monthconvert(char(list(i,2)))};
        end
        order = dateorder(list);
        db3ind = db3ind(order);

        fcsv.current = {};
        fcsv.voltage = {};
        fcsv.name = {};
        fcsv.sr = 0.1; %Hard coded sampling rate
        smpV = 50; %Hard coded smoothing
        for i = 1:size(db3ind)
            [I,V] = readFCSV(D(db3ind(i,1)).name,fcsvdir,smpV);
            fcsv.voltage{i,1} = V;
            fcsv.current{i,1} = I;
            fcsv.name{i,1} = D(db3ind(i)).name;
        end
        disp('saving fcsv data...')
        save(fcsvmatfile,'fcsv')
    end
end