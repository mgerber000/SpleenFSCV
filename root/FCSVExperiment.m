classdef FCSVExperiment < dynamicprops
    properties
        inputdir
        physdata
        fcsvdata
        usephysio
    end
    
    properties (Constant)
        repsetregexp = '.*(?<!iso.*)(?<!until.*)(?<!x\s*)(\d+\.\d+)(?!x).*'
        repsetxregexp = '.*(?<!x\s?set\s?id\s*)(?<!iso.*)(?<!until.*)(?<!x\s*)(\d+\.\d+)(?!x).*'
    end
    
    methods (Access = public)
        function obj = FCSVExperiment
            %Load adi file dir
            usebr = false;
            fdir = uigetdir('F:\FSCV\NE-ELISA');
            files = dir(fdir);
            for i = 1:size(files,1)
                filenames{i} = files(i).name;
            end
            adifiles = strfind(filenames,'.adicht');
            for i = 1:size(adifiles,2)
                temp(i) = ~isempty(adifiles{1,i});
            end
            if ~isempty(find(temp,1))
                fname = filenames{find(temp,1)};
            else
                fname = {};
            end
            obj.inputdir = fdir;
            
            if isempty(fname)
                obj.physdata = [];
                obj.usephysio = false;
            else
                %Check if adimat file already exists and load physdata
                adifile = strcat(obj.inputdir,'/',fname);
                adimatfile = [adifile(1:end-6), 'mat']; 
                if exist(adimatfile, 'file')
                    obj.physdata = PhysioData(adimatfile, obj.repsetregexp, obj.repsetxregexp);
                    obj.physdata = hrcalc(obj.physdata);
                    if usebr
                        obj.physdata = brcalc(obj.physdata);
                    end
                elseif exist(strcat(obj.inputdir,'/',fname), 'file')
                    obj.physdata = PhysioData(adifile, obj.repsetregexp, obj.repsetxregexp);
                    obj.physdata = hrcalc(obj.physdata);
                    if usebr
                        obj.physdata = brcalc(obj.physdata);
                    end
                else
                    error('couldn''t find %s', adifile);
                end
                obj.usephysio = true;
            end
            
            inputdir = strsplit(obj.inputdir,'/');
            fname = genvarname(inputdir{1,end});
                
            %Check if fcsvmatfile exists already and load fcsv data
            obj.fcsvdata = extract_fcsv(obj.inputdir,fname);

            %Find stim triggers
            if obj.usephysio == true
                 obj = findtrigs(obj);
            end
            
            if obj.usephysio == true
                %Compile ECG and HR data
                ecgid = find(contains({obj.physdata.channel_meta.name},'ECG'));
                ecg_sr = obj.physdata.channel_meta(ecgid).dt(1,1);
                ecgchan = [];
                hrchan = [];
                for j = 1:size(obj.physdata.channel_meta(ecgid).n_samples,2)
                    ecgchan = [ecgchan ; obj.physdata.(sprintf('data__chan_%d_rec_%d',ecgid,j))];
                    locs = round(obj.physdata.thr{1,j}/ecg_sr);
                    locs = [locs size(obj.physdata.(sprintf('data__chan_%d_rec_%d',ecgid,j)),1)];
                    hr = [];
                    hr(1:locs(1,1)) = NaN;
                    for k = 1:size(locs,2)-1
                        hr(locs(1,k):locs(1,k+1)) = obj.physdata.hr{1,j}(1,k);
                    end
                    hrchan = [hrchan ; hr'];
                end
                if usebr
                    respid = find(contains({obj.physdata.channel_meta.name},'Nasal sensor'));
                    resp_sr = obj.physdata.channel_meta(respid).dt(1,1);
                    respchan = [];
                    brchan = [];
                    clear locs
                    for j = 1:size(obj.physdata.channel_meta(respid).n_samples,2)
                        respchan = [respchan ; obj.physdata.(sprintf('data__chan_%d_rec_%d',respid,j))];
                        locs = round(obj.physdata.tbr{j,1}/resp_sr);
                        locs = [locs size(obj.physdata.(sprintf('data__chan_%d_rec_%d',respid,j)),1)];
                        br = [];
                        br(1:locs(1,1)) = NaN;
                        for k = 1:size(locs,2)-1
                            if k == size(locs,2)-1
                                br(locs(1,k):locs(1,k+1)) = obj.physdata.br{j,1}(1,k);
                            else
                                br(locs(1,k):locs(1,k+1)) = obj.physdata.br{j,1}(1,k+1);
                            end
                        end
                        brchan = [brchan ; br'];
                    end
                    obj.physdata.respchan = respchan;
                    obj.physdata.brchan = brchan;
                end
                obj.physdata.ecgchan = ecgchan;
                obj.physdata.hrchan = hrchan;
            end
        end
        
        function obj = findtrigs(obj)
            ichan = find(cellfun(@(x) contains_(x, 'train', 'IgnoreCase', true), {obj.physdata.channel_meta.name}));            
            assert(length(ichan) == 1, 'did not uniquely identify the trigger channel');
            obj.physdata.trigid = ichan;
            obj.physdata.triggers = [];
            obj.physdata.trigchan = [];
            for irecord = 1:obj.physdata.file_meta.n_records
                trigchan = obj.physdata.(sprintf('data__chan_%d_rec_%d', ichan, irecord));
                obj.physdata.trigchan = [obj.physdata.trigchan; trigchan];
            end    
                
            inds1 = find(obj.physdata.trigchan(1:end-1, 1) < 2 & obj.physdata.trigchan(2:end, 1) >= 2);  % rising edges
            inds2 = find(obj.physdata.trigchan(1:end-1, 1) >= 2 & obj.physdata.trigchan(2:end, 1) < 2);  % falling edges

            if ~isempty(inds1) && ~isempty(inds2)
                if inds1(1) > inds2(1)
                    inds2 = inds2(2:end);
                end
                if inds1(end) > inds2(end)
                    inds1 = inds1(1:end-1);
                end
            elseif isempty(inds1) || isempty(inds2)
                inds1 = [];
                inds2 = [];
            end
            assert(length(inds1) == length(inds2));

            % assert(length(inds1) == length(inds2));
            if isempty(inds1) && isempty(inds2)
            else
                obj.physdata.triggers = [obj.physdata.triggers; [inds1,inds2]];
            end
        end
    end
end     