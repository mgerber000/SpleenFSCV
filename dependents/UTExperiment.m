classdef UTExperiment < handle
    % We will extrat each rhs segment separately and reference the
    % appropriate ADI record for comments
    % |--rhs 1--|   |---rhs 2---| |------rhs 3------| |--rhs 4--|
    %    |------ADI 1------|  |-ADI 2-|   |-----------ADI 3-----------|
    % multiple files can be contiguous. If t(1) == 0 it is restarted
    properties
        sr                  % sampling rate (Hz)
        neuraldata          % numsamples x numchannels array of neural data
        intanchannel        % definition of x numchannels of neural data
        trigchan            % numsamples x 1 array indicating pulse trains
        segtimes            % datenum from filename at start of segment
        trains              % array of PulseTrain objects
        physdata            % PhysioData object containing ADI data
        fcsv                % FCSV data 
        wparams             % WaveformParams object containing the Matlab script parameters
        inputdir            % directory name
        repsetlist          % list of unique repsets
        connectedelectrodes   % list of elec
        distance            % distance between cuffs
        feature_metric      % neural feature metric for plotting
        modeling            % modeling results
        results     % results table with mean values for physiological properties for pre, during and post stim
    end
    
    properties (Constant)
        repsetregexp = '.*(?<!iso.*)(?<!until.*)(?<!x\s*)(\d+\.\d+)(?!x).*'
        repsetxregexp = '.*(?<!x\s?set\s?id\s*)(?<!iso.*)(?<!until.*)(?<!x\s*)(\d+\.\d+)(?!x).*'
    end
    
    methods (Access = public)
        function obj = UTExperiment(experimentdir, connectedelectrodes, applynotch, filePrefix, useneural, usephysio)
            % load the neural data from the converted mat file
            obj.inputdir = experimentdir;
            obj.connectedelectrodes = connectedelectrodes;
            if ~exist('filePrefix','var'), filePrefix = []; end
            if ~exist('useneural','var') || isempty(useneural), useneural = true; end
            if ~exist('usephysio','var') || isempty(usephysio), usephysio = true; end
            
            if useneural
                % convert the rhs files if necessary
                Drhs = dir(sprintf('%s\\%s*.rhs', obj.inputdir, filePrefix));
                Dmat = dir(sprintf('%s\\%s*.mat', obj.inputdir, filePrefix));
                if isempty(Drhs) && ~isempty(Dmat)
                    % mat files only (intan files end in 'yymmdd_HHMMSS.mat')
                elseif ~isempty(Drhs) && isempty(Dmat)
                    % rhs files only (convert to mat)
                    fprintf('converting rhs files\n');
                    if ~exist('applynotch','var') || isempty(applynotch) || ~ismember(applynotch, {'Yes','No'})
                        applynotch = questdlg('Would you like to apply a 60 Hz notch filter?', 'Notch Filter', 'Yes', 'No', 'Yes');
                    end
                    read_intan_batch(obj.inputdir, applynotch, [], connectedelectrodes);
                    Dmat = dir(sprintf('%s\\%s*.mat', obj.inputdir, filePrefix));
                elseif ~isempty(Drhs) && ~isempty(Dmat)
                    % both mat and rhs files (verify each rhs file corresponds to a mat file)
                    % convert any rhs files that don't correspond to a mat file
                    filenames = [];
                    for irhs = 1:length(Drhs)
                        rhsfile = sprintf('%s\\%s', obj.inputdir, Drhs(irhs).name);
                        matfile_ = sprintf('%s\\%s', obj.inputdir, [Drhs(irhs).name(1:end-3) 'mat']);
                        if ~exist(matfile_, 'file')
                            filenames = [filenames {rhsfile}];  %#ok
                        end
                    end
                    if ~isempty(filenames)
                        % process any remaining files
                        if ~exist('applynotch','var') || isempty(applynotch) || ~ismember(applynotch, {'Yes','No'})
                            applynotch = questdlg('Would you like to apply a 60 Hz notch filter?', 'Notch Filter', 'Yes', 'No', 'Yes');
                        end
                        read_intan_batch(obj.inputdir, applynotch, filenames, connectedelectrodes);
                    end
                    Dmat = dir(sprintf('%s\\%s*.mat', obj.inputdir, filePrefix));
                end

                % filter out mat files that don't correspond to rhs files
                inds = cellfun(@(x) ~isempty(x), regexpi({Dmat.name}, '.*\d{6}_\d{6}.mat', 'match'));
                Dmat = Dmat(inds);
                % RHS files don't need to be present to do this
                % inds = ismember({Dmat.name}, cellfun(@(x) [x(1:end-3) 'mat'], {Drhs.name}, 'UniformOutput', false));
                % Dmat = Dmat(inds);

                % verify that the files are in chronological order
                t_ = arrayfun(@(x) datenum(x.name(end-16:end-4), 'yymmdd_HHMMSS'), Dmat);
                if ~issorted(t_)
                    fprintf('Files were not in chronological order. Rearranging.\n');
                    [t_, si] = sort(t_);
                    Dmat = Dmat(si);
                end

                if length(Dmat) == 0
                    fprintf('directory empty\n');
                    return
                end

                % load the neural data and the stim trigger channels
                iseg = 0;
                for ifile = 1:length(Dmat)
                    statusstr = sprintf('reading file %d of %d', ifile, length(Dmat));
                    fprintf('%s', statusstr);
                    filename = sprintf('%s\\%s', obj.inputdir, Dmat(ifile).name);
                    hmatfile = matfile(filename, 'Writable', false);

                    if isprop(hmatfile, 'y') && isprop(hmatfile, 'trig')
                        % This was used when I concatenated the data, but now I don't want to make assumptions about continuity
                        obj.neuraldata = hmatfile.y;
                        obj.trigchan = hmatfile.trig;
                    elseif isprop(hmatfile, 'board_dig_in_data') && isprop(hmatfile, 'amplifier_data')
                        temp = hmatfile.t;
                        t0 = temp(1);
                        if t0 == 0 || iseg == 0
                            iseg = iseg + 1;
                            obj.segtimes(iseg) = t_(ifile);
                        end

                        chans = 1:size(hmatfile.amplifier_data, 1);
                        % check ~connectedelectrodes and verify they're noise
                        % hmatfile.amplifier_data = hmatfile.amplifier_data(connectedelectrodes);

                        if length(obj.neuraldata) < iseg
                            obj.neuraldata{iseg} = hmatfile.amplifier_data';
                            obj.trigchan{iseg} = hmatfile.board_dig_in_data';
                            obj.intanchannel{iseg} = hmatfile.amplifier_channels;
                        else
                            obj.neuraldata{iseg} = [obj.neuraldata{iseg}; hmatfile.amplifier_data'];
                            obj.trigchan{iseg} = [obj.trigchan{iseg}; hmatfile.board_dig_in_data'];
%                             obj.intanchannel{iseg} = [obj.trigchan{iseg} hmatfile.amplifier_channels'];
                        end

                    end
                    obj.sr = 30e3;
                    fprintf(repmat('\b', [1 length(statusstr)]));
                end
            else
                obj.sr = 30e3;  % This is used to construct fake data so it should be set
            end
            
            if usephysio
                % load the physiological data from the LabChart file
                Dadi = dir(sprintf('%s\\%s*.adicht', obj.inputdir, filePrefix));
                
                obj.physdata = [];
                adit0 = [];
                for iadi = 1:length(Dadi)
                    if ~contains_(Dadi(iadi).name, 'pre', 'IgnoreCase', true)
                        adifile = sprintf('%s\\%s', obj.inputdir, Dadi(iadi).name);
                        adimatfile = sprintf('%s\\%s', obj.inputdir, [Dadi(iadi).name(1:end-6) 'mat']);
                        if exist(adimatfile, 'file')
                            obj.physdata = [obj.physdata PhysioData(adimatfile, obj.repsetregexp, obj.repsetxregexp)];
                        elseif exist(adifile, 'file')
                            obj.physdata = [obj.physdata PhysioData(adifile, obj.repsetregexp, obj.repsetxregexp)];
                        else
                            error('couldn''t find %s', adifile);
                        end
                        % if length(obj.physdata) > 1
                        % verify that adicht files are in chronological order
                        % recmeta = obj.physdata(end - 1).record_meta(end);
                        adit0 = [adit0 obj.physdata(end).record_meta(1).record_start];
                        % adit0 = [adit0 recmeta.record_start + (recmeta.n_ticks * recmeta.tick_dt)/86400];
                        % assert(recmeta.record_start + (recmeta.n_ticks * recmeta.tick_dt)/86400 < ...
                        %     obj.physdata(end).record_meta(1).record_start);
                        % end
                    end
                end
                
                if verLessThan('matlab', '9.1')
                    if ~issorted(adit0)
                        fprintf('adi files were not read in chronological order. Reordering them.\n');
                        [~, si] = sort(adit0);
                        obj.physdata = obj.physdata(si);
                    end
                else
                    if ~issorted(adit0, 'ascend')
                        fprintf('adi files were not read in chronological order. Reordering them.\n');
                        [~, si] = sort(adit0);
                        obj.physdata = obj.physdata(si);
                    end
                end
                
                % still need segtimes(1)
                useneuralsegtime = false;
                if useneuralsegtime
                    Dmat = dir(sprintf('%s\\%s*.rhs', obj.inputdir, filePrefix));
                    inds = cellfun(@(x) ~isempty(x), regexpi({Dmat.name}, '.*\d{6}_\d{6}.mat', 'match'));
                    Dmat = Dmat(inds);
                    if isempty(Dmat)
                        Drhs = dir(sprintf('%s\\%s*.rhs', obj.inputdir, filePrefix));
                        inds = cellfun(@(x) ~isempty(x), regexpi({Drhs.name}, '.*\d{6}_\d{6}.rhs', 'match'));
                        Drhs = Drhs(inds);
                        t_ = arrayfun(@(x) datenum(x.name(end-16:end-4), 'yymmdd_HHMMSS'), Drhs);
                    else
                        t_ = arrayfun(@(x) datenum(x.name(end-16:end-4), 'yymmdd_HHMMSS'), Dmat);
                    end
                    obj.segtimes = min(t_);
                else
                    % completely decouple from neural
                    obj.segtimes = arrayfun(@(x) obj.physdata(x).record_meta(1).record_start, 1:length(obj.physdata));
                end
            end
            
            % load the setID file
            hasparamsfile = true;
            Dcsv = dir(sprintf('%s\\*.csv', obj.inputdir));
            if isempty(Dcsv)
                Dxls = dir(sprintf('%s\\*.xls*', obj.inputdir));
                if isempty(Dxls)
                    hasparamsfile = false;
                else
                    assert(length(Dxls) == 1, 'can only be one xls* file');
                    wparamsfile = sprintf('%s\\%s', obj.inputdir, Dxls.name);
                end
            else
                assert(length(Dcsv) == 1, 'can only be one csv file');
                wparamsfile = sprintf('%s\\%s', obj.inputdir, Dcsv.name);
            end
            if hasparamsfile
                obj.wparams = WaveformParams(wparamsfile);
            end
                        
            % I don't want the constructor to fail from adjusting the segment times
            try
                % make sure the segtimes have been adjusted
                obj.adjust_segtimes();
            catch
                fprintf('adjust_segtimes() failed in the constructor\n');
            end
        end
        
        function extract_trains(obj, minT, T1, T2, useneural, usephysio, plotflag)
            % assert(isa(chanfcn, 'function_handle'));
            if ~exist('plotflag', 'var') || isempty(plotflag)
                plotflag = false;
            end
            
            if ~exist('minT', 'var') || isempty(minT)
                minT = 1; % Only consider trains that are longer than minT seconds
            end
            if ~exist('T1', 'var') || isempty(T1)
                T1 = 0.010; % Get T1 seconds before the start of the train
            end
            if ~exist('T2', 'var') || isempty(T2)
                T2 = 0.020; % Get T2 seconds after the end of the train
            end
            if ~exist('usephysio', 'var') || isempty(usephysio)
                usephysio = false;
            end
            
            obj.trains = [];
            obj.repsetlist = [];            
            if ~isempty(obj.physdata)
                % support array of PhysioData objects
                [setIDcomments, setidinds, tcomments, repsets_] = deal(cell(1, length(obj.physdata)));
                [releccomments, relecinds, tcommentsrelec, relecs] = deal(cell(1, length(obj.physdata)));
                [seleccomments, selecinds, tcommentsselec, selecs] = deal(cell(1, length(obj.physdata)));
                
                for iadi = 1:length(obj.physdata)
                    setIDcomments{iadi} = arrayfun(@(icomment) ~contains_(obj.physdata(iadi).comments(icomment).str, 'iso', 'IgnoreCase', true) && ...
                        ~isempty(regexpi(obj.physdata(iadi).comments(icomment).str, obj.repsetregexp, 'start')), 1:length(obj.physdata(iadi).comments));
                    setidinds{iadi} = find(setIDcomments{iadi});
                    tcomments{iadi} = arrayfun(@(icomment) (obj.physdata(iadi).comment_time(icomment) - obj.segtimes(1)) * 24 * 60 * 60, setidinds{iadi});
                    repsets_{iadi} = arrayfun(@(ind) obj.physdata(iadi).parse_comment_str(ind), setidinds{iadi}, 'UniformOutput', false);
                    
                    % TODO also look for R= and S=
                    releccomments{iadi} = arrayfun(@(icomment) contains_(obj.physdata(iadi).comments(icomment).str, 'record from', 'IgnoreCase', true) || ...
                        contains_(obj.physdata(iadi).comments(icomment).str, 'R=', 'IgnoreCase', true), 1:length(obj.physdata(iadi).comments));
                    relecinds{iadi} = find(releccomments{iadi});
                    tcommentsrelec{iadi} = arrayfun(@(icomment) (obj.physdata(iadi).comment_time(icomment) - obj.segtimes(1)) * 24 * 60 * 60, relecinds{iadi});
                    relecs{iadi} = arrayfun(@(ind) obj.physdata(iadi).parse_elec_comment_str(ind, 'R'), relecinds{iadi}, 'UniformOutput', false);
                    
                    seleccomments{iadi} = arrayfun(@(icomment) contains_(obj.physdata(iadi).comments(icomment).str, 'stim from', 'IgnoreCase', true) || ...
                        contains_(obj.physdata(iadi).comments(icomment).str, 'S=', 'IgnoreCase', true), 1:length(obj.physdata(iadi).comments));
                    selecinds{iadi} = find(seleccomments{iadi});
                    tcommentsselec{iadi} = arrayfun(@(icomment) (obj.physdata(iadi).comment_time(icomment) - obj.segtimes(1)) * 24 * 60 * 60, selecinds{iadi});
                    selecs{iadi} = arrayfun(@(ind) obj.physdata(iadi).parse_elec_comment_str(ind, 'S'), selecinds{iadi}, 'UniformOutput', false);
                    % relecs{iadi} = {obj.physdata(iadi).comments(releccomments{iadi}).str};
                    % TODO This is what I expected to see but the actual comments were based on the above
                    % releccomments{iadi} = arrayfun(@(icomment) ~isempty(regexpi(obj.physdata(iadi).comments(icomment).str, 'R\d{1, 2}(?<,\s?\d{1, 2})*;S\d{1, 2}(?<,\s?\d{1, 2})*', 'start')), 1:length(obj.physdata(iadi).comments));
                end
            else
                tcomments = cell(1, length(obj.physdata));
            end
            
            if useneural && ~isempty(obj.trigchan)
                % let extracting trains from neural data take precedence over extracting trains from physiological data
                dtseg = (obj.segtimes - obj.segtimes(1)) * 24*60*60;
                obj.extract_neural_trains(T1, T2, minT, dtseg, tcomments, setidinds, plotflag);
            elseif usephysio && ~isempty(obj.physdata)
                % neural data is not loaded, use physio data
                obj.trains = cell(size(obj.physdata));
                
                curpol = 0;
                obj.repsetlist = [];
                for iseg = 1:length(obj.physdata)
                    [obj.trains{iseg}, repsetlist, curpol] = obj.physdata(iseg).extract_trains(T1, T2, minT, tcomments{iseg}, setidinds{iseg}, obj.wparams, curpol, obj.segtimes(1), plotflag);
                    obj.repsetlist = [obj.repsetlist repsetlist];
                end
%                 obj.extractphysdata();
            end
            
            % ignore the initial assigments for the trains' repset and pol and reassign them here using nearest neighbor
            % assert(~all(cellfun(@(x) isempty(x), repsets_)), 'no rep sets');
            if ~isempty(obj.trains)
                if ~all(cellfun(@(x) isempty(x), repsets_))
                    obj.associate_repset_comments_with_trains(tcomments, repsets_, obj.wparams);
                end
                % associate the trains with the recording electrodes                
                if ~isempty([relecs{:}])
                    obj.associate_elecs_with_trains(tcommentsrelec, relecs, 'relecs');
                else
                    warning('Missing rec electrodes for a train')
                end
                if ~isempty([selecs{:}])
                    obj.associate_elecs_with_trains(tcommentsselec, selecs, 'selecs');
                else
                    warning('Missing stim electrodes for a train')
                end
            end
            
            % copy physiological data to the trains
            if ~isempty(obj.physdata)
                obj.store_phys_in_trains();
            end
        end
        
        function associate_elecs_with_trains(obj, tcommentselec, elecs, myfield)
            % both t0 and tcomments are cell arrays of segments containing arrays of times in seconds
            t0 = cellfun(@(x) [x.t0], obj.trains, 'UniformOutput', false);
            
            % flatten the segments
            allt0 = cell2mat(t0);
            alltcomments = cell2mat(tcommentselec);
            allelecs = [elecs{:}];
            
            % perform nearest neighbor interpolation to find which trains the repset comments should be associated with
            assoctrain = interp1(allt0, 1:length(allt0), alltcomments, 'nearest', 'extrap');
            % assert(isequal(assoctrain, unique(assoctrain)));  % could be multiple comments before a train 
            if ~isequal(assoctrain, unique(assoctrain))
                fprintf('~isequal(assoctrain, unique(assoctrain))\n');
                [assoctrain, ia] = unique(assoctrain, 'last');
                allelecs = allelecs(ia);
            end
            
            trainct = 0;    % keep track of the flattened train count for checking the membership in assoctrain
            elecct = 0;    % keep track of the flattened recording electrode count for indexing into allrepsets
            for iseg = 1:length(obj.trains)
                for itrain = 1:length(obj.trains{iseg})
                    trainct = trainct + 1;
                    if trainct < min(assoctrain)
                        % trains prior to any recording electrode
                        obj.trains{iseg}(itrain).(myfield) = [];                        
                    elseif ismember(trainct, assoctrain)
                        % recording electrode setting detected
                        elecct = elecct + 1;                        
                        obj.trains{iseg}(itrain).(myfield) = allelecs{elecct};                    
                    else
                        % continuation of previous recording electrodes setting
                        obj.trains{iseg}(itrain).(myfield) = allelecs{elecct};
                    end
                end
            end
            assert(elecct == length(allelecs));
            assert(trainct == length(allt0));
        end
        
        function associate_repset_comments_with_trains(obj, tcomments, repsets, wparams)
            % both t0 and tcomments are cell arrays of segments containing arrays of times in seconds
            % TODO handle case where a cell is empty
            t0 = cell(1, length(obj.trains));
            for itrain = 1:length(obj.trains)
                if isa(obj.trains{itrain}, 'PulseTrain')
                    t0{itrain} = [obj.trains{itrain}.t0];                
                end
            end
            % t0 = cellfun(@(x) [x.t0], obj.trains, 'UniformOutput', false);            
            
            % flatten the segments
            allt0 = cell2mat(t0);
            alltcomments = cell2mat(tcomments);
            allrepsets = [repsets{:}];
            
            % perform nearest neighbor interpolation to find which trains the repset comments should be associated with
            assoctrain = interp1(allt0, 1:length(allt0), alltcomments, 'nearest', 'extrap');
            if ~isequal(assoctrain, unique(assoctrain))
                fprintf('~isequal(assoctrain, unique(assoctrain))\n');
                [assoctrain, ia] = unique(assoctrain);
                allrepsets = allrepsets(ia);
            end
            
            trainct = 0;    % keep track of the flattened train count for checking the membership in assoctrain
            repsetct = 0;   % keep track of the flattened repset count for indexing into allrepsets
            for iseg = 1:length(obj.trains)
                for itrain = 1:length(obj.trains{iseg})
                    trainct = trainct + 1;
                    if trainct < min(assoctrain)
                        % trains prior to any repset
                        obj.trains{iseg}(itrain).repset = '0.0';
                        obj.trains{iseg}(itrain).ampgain = 1;
                        obj.trains{iseg}(itrain).pol = 0;  % unknown polarity
                    elseif ismember(trainct, assoctrain)
                        % repset detected
                        repsetct = repsetct + 1;
                        pol = 1;
                        obj.trains{iseg}(itrain).repset = allrepsets{repsetct};
                        obj.trains{iseg}(itrain).ampgain = 1;
                        obj.trains{iseg}(itrain).pol = pol;
                    elseif ~isempty(wparams) && wparams.altpol
                        % alternating polarity
                        pol = -pol;
                        obj.trains{iseg}(itrain).repset = allrepsets{repsetct};
                        obj.trains{iseg}(itrain).ampgain = 1;
                        obj.trains{iseg}(itrain).pol = pol;
                    else
                        % constant polarity
                        obj.trains{iseg}(itrain).repset = allrepsets{repsetct};
                        obj.trains{iseg}(itrain).ampgain = 1;
                        obj.trains{iseg}(itrain).pol = 1;
                    end
                end
            end
            assert(repsetct == length(allrepsets));
            assert(trainct == length(allt0));
        end
        
        function plot_trains(obj)
            % plot the trains
            
            % This is a hack for when the neural data isn't loaded
            if isempty(obj.sr)
                sr = 1 / obj.physdata(1).record_meta(1).tick_dt;
            else
                sr = obj.sr;
            end
            
            % subplotsz = ceil(sqrt(length(obj.trains)));
            figure;ax = gca;hold(ax, 'on');
            dtseg = (obj.segtimes - obj.segtimes(1)) * 24*60*60;
            for iseg = 1:length(obj.trains)
                
                % This is a hack for when the neural data isn't loaded
                if isempty(obj.neuraldata)                    
                    szneuraldata = 0;
                    dt1 = obj.physdata(iseg).record_meta(1).tick_dt;
                    for irecord = 1:obj.physdata(iseg).file_meta.n_records                        
                        dt = obj.physdata(iseg).record_meta(irecord).tick_dt;
                        szneuraldata = szneuraldata + round(dt/dt1 * obj.physdata(iseg).record_meta(irecord).n_ticks);
                    end
                else
                    szneuraldata = size(obj.neuraldata{iseg}, 1);
                end
            
                % form train interval vector
                if isa(obj.trains{iseg}, 'PulseTrain')
                    t0 = [obj.trains{iseg}.t0];  % t0 is already referenced to the start of the first segment
                    dt = arrayfun(@(x) length(x.data)/sr, obj.trains{iseg});
                    t{iseg} = [dtseg(iseg) reshape([t0;t0;t0+dt;t0+dt], 1, []) dtseg(iseg) + szneuraldata/sr];
                    y{iseg} = [0 reshape(repmat([0;1;1;0], [1 length(t0)]), 1, []) 0];
                else
                    t{iseg} = dtseg(iseg) + [0 szneuraldata/sr];
                    y{iseg} = [0 0];
                end
                
                % plot the trains
                plot(ax, t{iseg}/60, y{iseg});
                % superimpose the pulses
                arrayfun(@(x) x.plot_pulses(ax), obj.trains{iseg});
            end
            
            xlabel('Time (min)');
            ylim([0 1.1]);
            
            % add the comments
            % TODO support array of PhysioData objects
            for iadi = 1:length(obj.physdata)
                for icomment = 1:length(obj.physdata(iadi).comments)
                    tc = obj.physdata(iadi).comment_time(icomment);
                    x = (tc - obj.segtimes(1)) * 24*60;
                    text(ax, x, 0, ...
                        obj.physdata(iadi).comments(icomment).str, ...
                        'Interpreter', 'none', 'Rotation', 90);
                end
            end
        end
        
        function [y, ax] = plot_queried_pulses(obj, pulseinds, plottype, varargin)
            % pulseinds        - pulse indices to use for every train
            % plottype         - 'all', 'avg', 'none'
            % Name/Value pairs
            % repset           - cell array of repset strings to match
            % pol              - array of polarities to match -1, 1, or [-1 1]
            % trainnumberinset - array of train numbers (e.g. [1:5 15:30 50:80]) if you know the train numbers to use)
            % whereinsession   - [start time, stop time] (seconds relative to first segment)
            % relec            - arry of recording electrodes
            
            % query the train indices
            traininds = obj.query_trains(varargin{:});
            [y, ax] = obj.plot_by_train(pulseinds, plottype, traininds);
        end
        
        function set_colors(obj, repsets, pols, cols)
            assert(length(repsets) * length(pols) == size(cols, 1));
            assert(size(cols, 2) == 3);
            
            ct = 1;
            for irepset = 1:length(repsets)
                for ipol = 1:length(pols)
                    obj.set_color(repsets{irepset}, pols(ipol), cols(ct, :));
                    ct = ct + 1;
                end
            end
        end
        
        function set_color(obj, repset, pol, col)
            % set the color for each train
            for iseg = 1:length(obj.trains)
                for itrain = 1:length(obj.trains{iseg})
                    if obj.trains{iseg}(itrain).pol == pol && strcmp(obj.trains{iseg}(itrain).repset, repset)
                        obj.trains{iseg}(itrain).col = col;
                    end
                end
            end
        end
        
        function [y, ax] = plot_by_train(obj, chanfcn, pulseinds, plottype, traininds, suppressartifact, extendsuppression, outlierpct, plotflag, detrendartifact, T)
            if ~exist('suppressartifact', 'var') || isempty(suppressartifact)
                suppressartifact = false;
            end
            if ~exist('extendsuppression', 'var') || isempty(extendsuppression)
                extendsuppression = 0;  % microseconds
            end
            if ~exist('plotflag', 'var') || isempty(plotflag)
                plotflag = true;
            end

            if ~exist('detrendartifact', 'var') || isempty(detrendartifact)
                detrendartifact = false;
            end
            if ~exist('T', 'var') || isempty(T)                
                T = 0.005;                
            end

            % loop over all trains and plot the waveforms
            vis = containers.Map([true false], {'on', 'off'});
            h = figure('Visible', vis(plotflag));ax = gca;
            for iseg = 1:length(traininds)
                inds = traininds{iseg};
                for itrain = 1:length(inds)
                    
                    % I don't know how to assert that the channels are consistent because chanfcn is a function, not a set
                    % chanfcn(1:size(obj.neuraldata{iseg}, 2))
                    % obj.trains{iseg}.relec
                    
                    switch plottype
                        case 'all'
                            obj.trains{iseg}(inds(itrain)).plot_all_responses(ax, chanfcn, pulseinds, 0, extendsuppression, outlierpct, detrendartifact, T);
                            keepseg = iseg;
                            keeptrain = inds(itrain);
                        case 'avg'
                            avgonly = true;
                            obj.trains{iseg}(inds(itrain)).plot_avg_response(ax, chanfcn, pulseinds, avgonly, 0, extendsuppression, outlierpct, detrendartifact, T);
                            keepseg = iseg;
                            keeptrain = inds(itrain);
                        case 'avgstd'
                            avgonly = false;  % also show stds
                            obj.trains{iseg}(inds(itrain)).plot_avg_response(ax, chanfcn, pulseinds, avgonly, 0, extendsuppression, outlierpct, detrendartifact, T);
                            keepseg = iseg;
                            keeptrain = inds(itrain);
                        otherwise
                            % pass
                    end
                end
            end

            
            temp = get(findobj(ax, 'Type', 'Line'), 'YData');
            if iscell(temp)
                y = cell2mat(temp)';
            elseif isa(temp, 'double')
                y = temp;
            end
            % obj.add_legend(ax, ikey, pol, cols);
            
            if suppressartifact
                obj.trains{keepseg}(keeptrain).blank_artifact_with_patch(ax, extendsuppression*1e-3);
            end

            if ~plotflag
                close(h);
            end
            
            % This is separate code and could become inconsistent
            % [y, pols] = obj.stack_query(traininds, pulseinds);
        end
        
        function [traininds] = trains_after_comment(obj, iadi, icomment)
            traininds = cell(1, length(obj.segtimes));
            for iseg = 1:length(obj.segtimes)
                traintimes = obj.segtimes(iseg) + [obj.trains{iseg}(:).t0] / 86400;
                % don't use the halfwidth rule for this function 0 * ...
                halfwidth = 0 * arrayfun(@(trains) size(trains.data, 1) - trains.n1 - trains.n2, obj.trains{iseg}) / obj.sr / 2 / 86400;
                % support array of PhysioData objects
                if icomment < length(obj.physdata(iadi).comments)
                    traininds{iseg} = find(obj.physdata(iadi).comment_time(icomment) - halfwidth < traintimes & ...
                        obj.physdata(iadi).comment_time(icomment + 1) - halfwidth > traintimes);
                elseif iadi < length(obj.physdata)
                    traininds{iseg} = find(obj.physdata(iadi).comment_time(icomment) - halfwidth < traintimes & ...
                        obj.physdata(iadi + 1).comment_time(1) - halfwidth > traintimes);
                else
                    traininds{iseg} = find(obj.physdata(iadi).comment_time(icomment) - halfwidth < traintimes);
                end
            end
        end
        
        function traininds = query_trains(obj, varargin)
            % repset           - cell array of repset strings to match
            % pol              - array of polarities to match -1, 1, or [-1 1]
            % trainnumberinset - array of train numbers (e.g. [1:5 15:30 50:80]) if you know the train numbers to use)
            % whereinsession   - [start time, stop time] (seconds relative to first segment)
            
            % parse name/value pairs
            p = inputParser;  % default values of [] are handled below in the segment loop
            addParameter(p, 'repset', [], @(x) isempty(x) || iscellstr(x) || isstr(x));
            addParameter(p, 'pol', [], @(x) isempty(x) || all(ismember(x, [-1 1])));
            addParameter(p, 'trainnumberinset', [], @(x) isempty(x) || isnumeric(x));
            addParameter(p, 'whereinsession', [], @(x) isempty(x) || isnumeric(x));
            addParameter(p, 'relec', [], @(x) isempty(x) || isnumeric(x));
            % addParameter(p, 'pulsesintrain', [], @isnumeric);
            parse(p, varargin{:});
            
            % for each segment of neural data
            % y = [];
            traininds = cell(1, length(obj.trains));
            for iseg = 1:length(obj.trains)
                repset = p.Results.repset;
                pol = p.Results.pol;            
                trainnumberinset = p.Results.trainnumberinset;
                whereinsession = p.Results.whereinsession;            
                relec = p.Results.relec;
                
                trainrepsets = arrayfun(@(x) x.repset, obj.trains{iseg}, 'UniformOutput', false);
                assert(iscellstr(trainrepsets));  %#ok
                
                trainpols = arrayfun(@(x) x.pol, obj.trains{iseg});
                trainnums = 1:length(obj.trains{iseg});
                trainwis = [obj.trains{iseg}(1).tlocal(1) obj.trains{iseg}(end).tlocal(end)];
                trainrelec = {obj.trains{iseg}.relecs};
                
                % default to all
                if ~exist('repset', 'var') || isempty(repset)
                    repset = trainrepsets;
                end
                if ~exist('pol', 'var') || isempty(pol)
                    pol = trainpols;
                end
                if ~exist('trainnumberinset', 'var') || isempty(trainnumberinset)
                    trainnumberinset = trainnums;
                end
                if ~exist('whereinsession', 'var') || isempty(whereinsession)
                    whereinsession = trainwis;
                end
                if ~exist('relec', 'var') || isempty(relec)
                    relec = [];
                end
                
                traininds{iseg} = find(ismember(trainrepsets, repset) & ismember(trainpols, pol) & ismember(trainnums, trainnumberinset) & ...
                    [obj.trains{iseg}.t0] >= whereinsession(1) & [obj.trains{iseg}.t0] <= whereinsession(2) & ...
                    cellfun(@(x) all(ismember(relec, x)), trainrelec));
            end
        end % of function query_trains
        
        function [offsetDays] = adjust_segtimes(obj)
            % Intan data has timestamps with 1 second precision
            % Physiological data has timestamps with 1 millisecond precision
            % Interpolate the physiodata so that all trains align (compensate for drift)
            
            % TODO align trains even when there aren't the same number of trains
            
            % get the times of all train rising edges from the Intan data
            numsegs = length(obj.trigchan);
            tre = cell(1, numsegs);
            % figure;hold on;
            for iseg = 1:numsegs
                % get the start time (to 1 sec precision) from the Intan data
                t0 = obj.segtimes(iseg);
                % find rising edges in the train channel
                risingedges = find(obj.trigchan{iseg}(1:end-1, 1) == 0 & obj.trigchan{iseg}(2:end, 1) == 1);
                assert(nnz(obj.trigchan{iseg}(1:end-1, 1) == 0 & obj.trigchan{iseg}(2:end, 1) == 1) < ...
                    nnz(obj.trigchan{iseg}(1:end-1, 2) == 0 & obj.trigchan{iseg}(2:end, 2) == 1), ...
                    'assumed train channel came before pulse channel');
                tre{iseg} = t0 + risingedges/obj.sr/86400;
                % reduce_plot((t0 - obj.segtimes(1))*86400 + (0:size(obj.trigchan{iseg}, 1)-1)/obj.sr, obj.trigchan{iseg}(:, 1));
            end
            % hold off;
            try
                tre_ = cell2mat(tre')';
            catch
                tre_ = cell2mat(tre);
            end
            
            % get the times of all train rising edges from the physiological data
            numsegs2 = length(obj.physdata);
            tre2 = cell(1, numsegs2);
            % figure;hold on;
            for iseg2 = 1:numsegs2
                nrecords = obj.physdata(iseg2).file_meta.n_records;
                ichan = find(cellfun(@(x) contains_(x, 'train', 'IgnoreCase', true), {obj.physdata(iseg2).channel_meta.name}));
                assert(~isempty(ichan), 'couldn''t find train sync in physiological data records');
                for irecord = 1:nrecords
                    tre2{iseg2} = cell(1, nrecords);
                    t02 = obj.physdata(iseg2).record_meta(irecord).record_start;
                    datachanname = sprintf('data__chan_%d_rec_%d', ichan, irecord);
                    datachan = obj.physdata(iseg2).(datachanname);
                    sr2 = 1/obj.physdata(iseg2).channel_meta(ichan).dt(irecord);
                    % there is a ramp time (so use a single threshold, not a range)
                    risingedges2 = find(datachan(1:end-1) < 2 & datachan(2:end) >= 2);  % 3.3 V logic
                    tre2{iseg2}{irecord} = [t02 + risingedges2/sr2/86400];
                    % reduce_plot((t02 - obj.segtimes(1))*86400 + (0:length(datachan)-1)/sr2, datachan);
                end
            end
            % hold off;
            try
                tre2_ = cell2mat(cellfun(@(x) cell2mat(x), tre2, 'UniformOutput', false)')';
            catch
                % Subsequent calls used this
                tre2_ = cell2mat(cellfun(@(x) cell2mat(x), tre2, 'UniformOutput', false));
            end
            
            % determine the subsecond offset (otherwise mutliple segments won't be consistent to 1 ms)
            % figure;plot((tre_ - tre2_)*86400, '.');

            
            % numedges2 = cumsum(cell2mat(cellfun(@(x) cellfun(@(y) length(y), x), tre2, 'UniformOutput', false)));            
            % figure;plot((tre - tre{1}(1))*86400, (tre2 - tre)*86400, '.');  % tre is now a cell array            
            
            % align trains and exclude trains that don't have a match
            % assert(length(tre_) == length(tre2_), 'not a 1:1 mapping of edges');  % don't require a 1:1 mapping 
            D = pdist2(reshape(tre_, [], 1), reshape(tre2_, [], 1)) * 86400 < 2.0;  % 2 sec threshold
            assert(max(sum(D, 1)) <= 1 && max(sum(D, 2)) <= 1);
            inds1 = any(D, 2);
            inds2 = any(D, 1);
            % modify tre_ and tre and tre2_ and tre2
            ct = 1;

            for iseg = 1:length(tre)           
                inds = (0:length(tre{iseg})-1);
                tre{iseg} = tre{iseg}(inds1(ct + inds));
                ct = ct + length(inds);
            end
            tre_ = tre_(inds1);
            assert(sum(cellfun(@(x) length(x), tre)) == length(tre_));
            
            ct2 = 1;
            for iseg2 = 1:length(tre2)
                nrecords = obj.physdata(iseg2).file_meta.n_records;
                for irecord = 1:nrecords
                    inds = (0:length(tre2{iseg2}{irecord})-1);
                    tre2{iseg2}{irecord} = tre2{iseg2}{irecord}(inds2(ct2 + inds));
                    ct2 = ct2 + length(inds);
                end
            end
            tre2_ = tre2_(inds2);
            
            numedges = cumsum(cellfun(@(x) length(x), tre));
            % compute the subsecond offset for segtimes
            offsetDays = zeros(1, length(tre));
            for iseg = 1:length(tre)
                if iseg == 1
                    inds = 1:numedges(iseg);
                else
                    inds = numedges(iseg-1)+1:numedges(iseg);
                end
                offsets = reshape(tre{iseg}, [], 1) - reshape(tre2_(inds), [], 1);
                assert(range(offsets*86400) < 2.0, 'offsets span more than 2 sec');
                offsetDays(iseg) = offsets(1);
            end
            disp('');
            % adjust offset
            % assume no drift before first train and after last train
            obj.segtimes = obj.segtimes - offsetDays; % first train in each segment will align
            
            % interpolate physio data
            for ichan = 1:obj.physdata(iseg).file_meta.n_channels
                ct = 0;
                ct2 = 1;
                for iseg = 1:length(obj.physdata)
                    for irecord = 1:obj.physdata(iseg).file_meta.n_records
                        datachanname = sprintf('data__chan_%d_rec_%d', ichan, irecord);
                        data = obj.physdata(iseg).(datachanname);
                        dt = obj.physdata(iseg).channel_meta(ichan).dt(irecord);
                        sr = 1/dt;
                        t0 = obj.physdata(iseg).record_meta(irecord).record_start;
                        tre2b = tre2{iseg}{irecord};
                        inds = ct + (1:length(tre2b));
                        offset = offsetDays(ct2);  % subsequent calls to this function should result in offsetDays close to 0, so I need to remove the offset from tre and tre_
                        ct = ct + length(inds);
                        ct2 = ct2 + 1;
                        tref = obj.segtimes(1);
                        t = (t0-tref)*86400 + (0:length(data)-1)*dt;
                        ptime = reshape((tre2b-tref)*86400, [], 1);
                        ntime = reshape((tre_(inds)-offset-tref)*86400, [], 1);
                        t2 = interp1(ptime, ntime, t, 'linear', 'extrap');  % nonuniform neural times for the physiological data
                        obj.physdata(iseg).(datachanname) = interp1(t2, data, t, 'linear', 'extrap'); % interpolate to the presumed times
                    end
                end
            end
        end
        
        function [inds] = get_samples_from_time(obj, t, timeunits, datasource)
            % t         - length 2 array of a start time and stop time
            % timeunits - units for t from {'datenum', 'sec'} where datenum
            %             is an absolute time and sec is seconds relateive 
            %             to obj.segtimes(1)
            % datasource - source of data from {'intan' or 'adi'} where
            %              intan is the neural data and adi is the
            %              physiological data
            % inds is a structure array with the fields:
            %  seg      - segment number (could span multiple segments)
            %  record   - record number (could span multiple records)
            %             set to [] for 'intan'
            %  chaninds - start and stop indices for each channel as a numchans x 2 array
            
            switch timeunits
                case 'datenum'
                    % convert to seconds
                    t = (t - obj.segtimes(1)) * 86400;
                case 'sec'
                    % pass  % t = obj.segtimes(1) + t/86400;
                otherwise
                    fprintf('invalid time units\n');
                    return
            end
            
            switch datasource
                case 'intan'
                    inds = [];
                    for iseg = 1:length(obj.segtimes)
                        t0 = (obj.segtimes(iseg) - obj.segtimes(1)) * 86400;
                        tf = t0 + size(obj.neuraldata{iseg}, 1) / obj.sr;
                        
                        % Overlapping segment detected
                        if ~all(t < t0) && ~all(t > tf)
                            % append an inds structure
                            inds = [inds struct('seg', iseg, 'record', [], 'chaninds', zeros(1, 2))];
                            % compute the indices for each channel
                            % x = t0 + (0:size(obj.neuraldata{iseg}, 1)-1) / obj.sr;
                            %inds(end).chaninds(1, :) = interp1(x, 1:length(x), t, 'nearest', 'extrap');
                            if t(1) < t0 && t(2) >= t0
                                inds(end).chaninds(1, 1) = 1;
                            elseif t(1) >= t0 && t(1) <= tf
                                inds(end).chaninds(1, 1) = round((t(1) - t0) * obj.sr) + 1;
                            else
                                fprintf('unhandled exception\n');
                            end
                            if t(2) > tf && t(1) <= tf
                                inds(end).chaninds(1, 2) = size(obj.neuraldata{iseg}, 1);
                            elseif t(2) >= t0 && t(2) <= tf
                                inds(end).chaninds(1, 2) = round((t(2) - t0) * obj.sr) + 1;
                            else
                                fprintf('unhandled exception\n');
                            end
                        end
                    end
                case 'adi'
                    inds = [];
                    for iseg = 1:length(obj.physdata)
                        physdat = obj.physdata(iseg);                        
                        for irecord = 1:physdat.file_meta.n_records
                            recmet = physdat.record_meta(irecord);
                            t0rec = (recmet.record_start - obj.segtimes(1)) * 86400;
                            tfrec = t0rec + recmet.n_ticks * recmet.tick_dt;                            
                            
                            % Overlapping segment detected
                            if ~all(t < t0rec) && ~all(t > tfrec)
                                % append an inds structure
                                inds = [inds struct('seg', iseg, 'record', irecord, 'chaninds', zeros(physdat.file_meta.n_channels, 2))];
                                % compute the indices for each channel
                                for ichan = 1:physdat.file_meta.n_channels
                                    dt = physdat.channel_meta(ichan).dt(irecord);
                                    % x = t0rec + (0:physdat.channel_meta(ichan).n_samples(irecord)-1) * dt;
                                    % inds(end).chaninds(ichan, :) = interp1(x, 1:length(x), t, 'nearest', 'extrap');
                                    if t(1) < t0rec && t(2) >= t0rec
                                        inds(end).chaninds(ichan, 1) = 1;
                                    elseif t(1) >= t0rec && t(1) <= tfrec
                                        inds(end).chaninds(ichan, 1) = round((t(1) - t0rec) / dt) + 1;
                                    else
                                        fprintf('unhandled exception\n');
                                    end
                                    if t(2) > tfrec && t(1) <= tfrec
                                        channame = sprintf('data__chan_%d_rec_%d', ichan, irecord);
                                        inds(end).chaninds(ichan, 2) = length(obj.physdata(iseg).(channame));
                                    elseif t(2) >= t0rec && t(2) <= tfrec
                                        inds(end).chaninds(ichan, 2) = round((t(2) - t0rec) / dt) + 1;
                                    else
                                        fprintf('unhandled exception\n');
                                    end
                                end
                            end

                        end
                    end
                otherwise
                    fprintf('invalid time units\n');
                    return
            end
        end
        
        function save_all_data(obj, matfilename, varargin)
            % stitch segments and write a mat file for each neural and physiological data channel
            % matfilename   - a filename that will be modified with a suffix for each channel
            % MultipleFiles - (default=true) write a separate mat file for each channel
            
            p = inputParser;
            addRequired(p, 'matfilename', @(x) ischar(x));
            addParameter(p, 'MultipleFiles', true, @(x) islogical(x));
            parse(p, matfilename, varargin{:});
            
            if isempty(matfilename)
                % file not specified, select from GUI
                if p.Results.MultipleFiles
                    putfiletitle = 'Save mat files';
                else
                    putfiletitle = 'Save mat file';
                end
                myfilter = {'*.mat', 'mat files (*.mat)'; ...
                    '*.*', 'All files (*.*)'};
                [myfile, mypath] = uiputfile(myfilter, putfiletitle, sprintf('%s\\data.mat', obj.inputdir));
                if myfile == 0
                    return
                end
                
                matfilename = [mypath myfile];
            end
            
            if p.Results.MultipleFiles
                sr = obj.sr;  % Intan sampling rate
                nsegs = length(obj.neuraldata);
                
                % neural data
                for ichan = 1:size(obj.neuraldata{1}, 2)
                    y = [];
                    for iseg = 1:nsegs
                        % stitch together each segment
                        y = [y; obj.neuraldata{iseg}(:, ichan)];
                        if iseg < nsegs
                            % pad time between segments with NaN
                            dtdays = obj.segtimes(iseg+1) - (obj.segtimes(iseg) + size(obj.neuraldata{iseg}, 1)/obj.sr/86400);
                            dtsamp = round(dtdays * 86400 * obj.sr);
                            y = [y; NaN(dtsamp, 1)];
                        end
                    end
                    % modify the matfile name
                    newmatfile = sprintf('%s_chan%d.mat', matfilename(1:end-4), ichan);
                    t0 = 0;  % neural data will start at t = 0.0 s
                    fprintf('writing %s\n', newmatfile);
                    % save a neural data channel
                    save(newmatfile, 'sr', 'y', 't0', '-v7.3');
                end
                
                % trigger channels
                for ichan = 1:size(obj.trigchan{1}, 2)
                    y = [];
                    for iseg = 1:nsegs
                        % stitch together each segment
                        y = [y; obj.trigchan{iseg}(:, ichan)];
                        if iseg < nsegs
                            % pad time between segments with NaN
                            dtdays = obj.segtimes(iseg+1) - (obj.segtimes(iseg) + size(obj.trigchan{iseg}, 1)/obj.sr/86400);
                            dtsamp = round(dtdays * 86400 * obj.sr);
                            y = [y; NaN(dtsamp, 1)];
                        end
                    end
                    % modify the matfile name
                    newmatfile = sprintf('%s_trig%d.mat', matfilename(1:end-4), ichan);
                    fprintf('writing %s\n', newmatfile);
                    % save a trigger channel from the Intan
                    t0 = 0;
                    save(newmatfile, 'sr', 'y', 't0', '-v7.3');
                end
                
                % comments
                comments = [];
                tcomment = [];
                for iseg = 1:length(obj.physdata)
                    comments = [comments {obj.physdata(iseg).comments.str}];
                    tcomment = [tcomment arrayfun(@(ii) obj.physdata(iseg).comment_time(ii), 1:length(obj.physdata(iseg).comments))];
                end
                tcomment = (tcomment - obj.segtimes(1))*86400;
                newmatfile = sprintf('%s_comments.mat', matfilename(1:end-4));
                save(newmatfile, 'tcomment', 'comments', '-v7.3');
                
                % physiological data
                for ichan = 1:obj.physdata(iseg).file_meta.n_channels
                    y = [];
                    % get all sample rates
                    sr = [];
                    for iseg = 1:length(obj.physdata)
                        for irecord = obj.physdata(iseg).file_meta.n_records
                            sr = [sr 1/obj.physdata(iseg).channel_meta(ichan).dt(irecord)];
                        end
                    end
                    sr_ = sprintf('%d ', unique(sr));
                    fprintf('%s\n', sr_);
                    sr_ = max(sr);
                    
                    for iseg = 1:length(obj.physdata)
                        % name = obj.physdata(iseg).channel_meta(ichan).name;
                        
                        for irecord = obj.physdata(iseg).file_meta.n_records
                            % TODO test with multiple records
                            % time in seconds releative to the start of the neural data (neural data has 1 sec precision and physdata has ms precision)
                            % nsamples = obj.physdata(iseg).channel_meta(ichan).n_samples(irecord);
                            % t0 = (obj.physdata(iseg).record_meta(irecord).record_start - obj.segtimes(1)) * 86400;
                            sr = 1/obj.physdata(iseg).channel_meta(ichan).dt(irecord);
                            
                            datachanname = sprintf('data__chan_%d_rec_%d', ichan, irecord);
                            data = obj.physdata(iseg).(datachanname);
                            if sr ~= sr_
                                [n, d] = rat(sr_ / sr);
                                y = [y; reshape(resample(data, n, d), [], 1)];
                            else
                                y = [y; reshape(data, [], 1)];
                            end
                            
                            if irecord < obj.physdata(iseg).file_meta.n_records
                                % pad NaNs
                                dtdays = obj.physdata(iseg).record_meta(irecord+1).record_start - (obj.physdata(iseg).record_meta(irecord).record_start + length(data)/sr/86400);
                                dtsamp = round(dtdays * 86400 * sr_);
                                y = [y; NaN(dtsamp, 1)];
                            end
                        end
                        
                        if iseg < length(obj.physdata)
                            % pad NaNs
                            dtdays = obj.physdata(iseg+1).record_meta(1).record_start - (obj.physdata(iseg).record_meta(end).record_start + length(data)/sr/86400);
                            dtsamp = round(dtdays * 86400 * sr_);
                            y = [y; NaN(dtsamp, 1)];
                        end
                    end
                    
                    t0 = (obj.physdata(1).record_meta(1).record_start - obj.segtimes(1)) * 86400;
                    sr = sr_;
                    newmatfile = sprintf('%s_%s.mat', matfilename(1:end-4), obj.physdata(iseg).channel_meta(ichan).name);
                    fprintf('writing %s\n', newmatfile);
                    save(newmatfile, 't0', 'sr', 'y', '-v7.3');
                end
            else
                error('single mat file not yet implemented');
                % save(matfilename, '-v7.3');
            end
            
        end
        
        %% Functions for processing physiological signals
        % Subash Padmanaban 08/06/2018
        function obj = computeheartrate(obj, method, minpeakprom, minpeakdist, plotvar)
            clc, sprintf('Computing heart rate...')
            for jj=1:numel(obj.physdata)
                obj.physdata(jj).getheartrate(method, minpeakprom, minpeakdist, plotvar);
            end
        end
        
        function obj = computeresprate(obj, minpeakprom, smoothvar, plotvar)
            clc, sprintf('Computing respiratory rate...')
            for jj=1:numel(obj.physdata)
                obj.physdata(jj).getresprate(minpeakprom, smoothvar, plotvar);
            end
        end
        
        function obj = extractphysdata(obj, usephysio)
            % Extract physiological data corresponding to each train. Store
            % baseline and during stimulation values for raw ECG, heart
            % rate, raw Resp signal, breathing rate.
            
            ecgchan = obj.physdata.ecgid;
            respchan = obj.physdata.respid;
            
            % Loop through each record
            for iseg=1:numel(obj.trains)
                
                % Loop through each train within a record
                for itrain=1:numel(obj.trains{iseg})
                    clc, sprintf('Extracting physiological features for record %d, train %d', iseg, itrain)
                    ctr = obj.trains{iseg}(itrain); % current train
                    
                    % change the inds extraction to include -x seconds
                    % before and + x seconds after stimulation
                    buffertime = 5; % 5 seconds before and after stim extracted
                    obj.trains{iseg}(itrain).buffertime = buffertime;
                    if usephysio
                        stimdur = (ctr.data(end) - ctr.data(1))./ctr.sr;
                        obj.trains{iseg}(itrain).stimdur = stimdur;
                        inds = [ctr.data(1)/ctr.sr - buffertime ctr.data(end)/ctr.sr + buffertime];
                    else
                        stimdur = ctr.tlocal(end) - ctr.tlocal(1);
                        obj.trains{iseg}(itrain).stimdur = stimdur;
                        inds = get_samples_from_time(obj, [(ctr.t0 - buffertime) (ctr.t0 + stimdur + buffertime)], 'sec','adi');
                    end
                    ecg = obj.physdata(iseg).ecg{ctr.record};
                    resp = obj.physdata(iseg).resp{ctr.record};
                    ecgsr = round(1/obj.physdata(iseg).channel_meta(ecgchan).dt(ctr.record));
                    bsr = round(1/obj.physdata(iseg).channel_meta(respchan).dt(ctr.record));
                    if any(inds < 1) || any(inds > numel(ecg)/ecgsr)
                        warning('Stimulation train does not have enough samples to extract..')
                        continue;
                    end
                    ecginds = round(inds.*ecgsr);
                    respinds = round(inds.*bsr);
                    obj.trains{iseg}(itrain).ecg = ecg(ecginds(1):ecginds(2));
                    obj.trains{iseg}(itrain).resp = resp(respinds(1):respinds(2));
                    
                    % create time vector for raw ecg and resp
                    % base of [-buffertime, stimdur+buffertime]
                    %                     obj.trains{rec}(tr).tecg = linspace(-buffertime, stimdur + buffertime, numel(obj.trains{rec}(tr).ecg));
                    %                     obj.trains{rec}(tr).tresp = linspace(-buffertime, stimdur + buffertime, numel(obj.trains{rec}(tr).resp));
                    
                    % base of [stimstart- buffertime, stimstart +
                    % buffertime]
                    obj.trains{iseg}(itrain).tecg = linspace(ecginds(1)/ecgsr, ecginds(2)/ecgsr, numel(obj.trains{iseg}(itrain).ecg));
                    obj.trains{iseg}(itrain).tresp = linspace(respinds(1)/bsr, respinds(2)/bsr, numel(obj.trains{iseg}(itrain).resp));
                    
                    % extract heart rate and breathing rate
                    hr = obj.physdata(iseg).hr{ctr.record};
                    thr = obj.physdata(iseg).thr{ctr.record};
                    hpks = obj.physdata(iseg).hpks{ctr.record};
                    hinds = thr > ecginds(1)/ecgsr & thr < ecginds(2)/ecgsr;
                    subset_thr = thr(hinds);
                    obj.trains{iseg}(itrain).thr = subset_thr;
                    %                     obj.trains{rec}(tr).thr = subset_thr - subset_thr(1) - buffertime;
                    obj.trains{iseg}(itrain).hr = hr(hinds);
                    obj.trains{iseg}(itrain).hpks = hpks(hinds);
                    
                    
                    br = obj.physdata(iseg).br{ctr.record};
                    tbr = obj.physdata(iseg).tbr{ctr.record};
                    bpks = obj.physdata(iseg).bpks{ctr.record};
                    binds = tbr > respinds(1)/bsr & tbr < respinds(2)/bsr;
                    obj.trains{iseg}(itrain).bpks = bpks(binds);
                    subset_tbr = tbr(binds);
                    obj.trains{iseg}(itrain).tbr = subset_tbr;
                    %                     obj.trains{rec}(tr).tbr = subset_tbr - subset_tbr(1) - buffertime;
                    obj.trains{iseg}(itrain).br = br(binds);
                end
            end
        end
        
        function obj = compute_physio_mean(obj, usephysio)
            % Loop through each train object and compute the mean of
            % physiological parameters. Store this as a table with mean
            % values [pre, during] based on polarity. Only store up to 4
            % trains of a single polarity.
            
            ecgchan = obj.physdata.ecgid;
            respchan = obj.physdata.respid;
            % Loop through each segment
            
            mydict = containers.Map();            
            for iseg=1:numel(obj.trains)
                
                % Loop through each train within a record
                for itrain=1:numel(obj.trains{iseg})
                    clc, sprintf('Computing mean for record %d, train %d', iseg, itrain)
                    ctr = obj.trains{iseg}(itrain);
                    if isempty(ctr.thr) || isempty(ctr.tbr)
                        T = table([NaN; NaN], [NaN; NaN],[NaN; NaN]);
                        T.Properties.VariableNames = {'Pre', 'During', 'Post'};
                        T.Properties.RowNames = {'hr','br'};
                        ctr.results = T;
                        continue
                    end
                    
                    if usephysio
                        stimdur = (ctr.data(end) - ctr.data(1))./ctr.sr;
                        obj.trains{iseg}(itrain).stimdur = stimdur;
                        inds = [ctr.data(1)/ctr.sr - ctr.buffertime ctr.data(end)/ctr.sr + ctr.buffertime];
                    else
                        stimdur = ctr.tlocal(end) - ctr.tlocal(1);
                        obj.trains{iseg}(itrain).stimdur = stimdur;
                        inds = get_samples_from_time(obj, [(ctr.t0 - buffertime) (ctr.t0 + stimdur + buffertime)], 'sec','adi');
                    end
                    ecgsr = round(1/obj.physdata(iseg).channel_meta(ecgchan).dt(ctr.record));
                    bsr = round(1/obj.physdata(iseg).channel_meta(respchan).dt(ctr.record));
                    ecginds = round(inds.*ecgsr);
                    respinds = round(inds.*bsr);
                    
                    % TODO should we use >= to cover every point? Or is it better to leave out points on the boundary?
                    indsprehr = ctr.thr < ctr.thr(1) + ctr.buffertime;
                    indsprebr = ctr.tbr < ctr.tbr(1) + ctr.buffertime;
                    indsduringhr = ctr.thr > ctr.thr(1) + ctr.buffertime & ctr.thr < ctr.thr(end) - ctr.buffertime;
                    indsduringbr = ctr.tbr > ctr.tbr(1) + ctr.buffertime & ctr.tbr < ctr.tbr(end) - ctr.buffertime;
                    indsposthr = ctr.thr > ctr.thr(end) - ctr.buffertime;
                    indspostbr = ctr.tbr > ctr.tbr(end) - ctr.buffertime;
                    
                    if ismember(ctr.repset, mydict.keys())                        
                        mydict(ctr.repset) = mydict(ctr.repset) + 1;
                    else
                        mydict(ctr.repset) = 1;
                    end
                    
                    
                    if stimdur > 2
                        T = table([nanmean(ctr.hr(indsprehr)); nanmean(ctr.br(indsprebr))],...
                            60.*[sum(indsduringhr)/stimdur; sum(indsduringbr)/stimdur],...
                            [nanmean(ctr.hr(indsposthr)); nanmean(ctr.br(indspostbr))]);
                        
                        % write results 2 structure
                        durparams = 60.*[sum(indsduringhr)/stimdur; sum(indsduringbr)/stimdur];
                        numbrpoints = sum(indsduringbr);
                        if numbrpoints == 0, numbrpoints = 1; end
                        numhrpoints = sum(indsduringhr);
                        if numhrpoints == 0, numhrpoints = 1; end
                        ctr.results2 = struct('HR_pre', ctr.hr(indsprehr), 'BR_pre', ctr.br(indsprebr), ...
                                          'HR_during', durparams(1)*ones(1,numhrpoints), 'BR_during', durparams(2)*ones(1,numbrpoints), ...
                                          'HR_post', ctr.hr(indsposthr), 'BR_post', ctr.br(indspostbr), ...
                                          'setIDcount', mydict(ctr.repset), ...
                                          'polarity', ctr.pol, ...
                                          'setid', ctr.repset);

                    elseif stimdur < 2 && stimdur > 0.5 % Greater than 500ms and 2 seconds
%                          T = table([60.*sum(indsprehr)/stimdur; nanmean(ctr.br(indsprebr))], ...
%                         [60.*sum(indsduringhr)/stimdur; nanmean(ctr.br(indsduringbr))], ...
%                         [60.*sum(indspostbr)/stimdur; nanmean(ctr.br(indspostbr))]);
                        postbr_points = find(indspostbr);
                        T = table([nanmean(ctr.hr(indsprehr)); nanmean(ctr.br(indsprebr))],...
                            [nanmean(ctr.hr(indsduringhr)); nanmean([ctr.br(indsduringbr) ctr.br(postbr_points(1))])],...
                            [nanmean(ctr.hr(indsposthr)); nanmean(ctr.br(indspostbr))]); 
                        
                        ctr.results2 = struct('HR_pre', ctr.hr(indsprehr), 'BR_pre', ctr.br(indsprebr), ...
                                          'HR_during', ctr.hr(indsduringhr), 'BR_during', [ctr.br(indsduringbr) ctr.br(postbr_points(1))], ...
                                          'HR_post', ctr.hr(indsposthr), 'BR_post', ctr.br(indspostbr), ...
                                          'setIDcount', mydict(ctr.repset), ...
                                          'polarity', ctr.pol, ...
                                          'setid', ctr.repset);

                    else % stim duration shorter than 500 ms
%                         T = table([nanmean(ctr.hr(indsprehr)); nanmean(ctr.br(indsprebr))], ...
%                         [nanmean(ctr.hr(indsduringhr)); nanmean(ctr.br(indsduringbr))], ...
%                         [nanmean(ctr.hr(indsposthr)); nanmean(ctr.br(indspostbr))]);
                    posthr_points = find(indsposthr);
                    postbr_points = find(indspostbr);
                        T = table([nanmean(ctr.hr(indsprehr)); nanmean(ctr.br(indsprebr))],...
                            [nanmean([ctr.hr(indsduringhr) ctr.hr(posthr_points(1))]); nanmean([ctr.br(indsduringbr) ctr.br(postbr_points(1))])],...
                            [nanmean(ctr.hr(indsposthr)); nanmean(ctr.br(indspostbr))]); 
                        
                        ctr.results2 = struct('HR_pre', ctr.hr(indsprehr), 'BR_pre', ctr.br(indsprebr), ...
                                          'HR_during', [ctr.hr(indsduringhr) ctr.hr(posthr_points(1))], 'BR_during', [ctr.br(indsduringbr) ctr.br(postbr_points(1))], ...
                                          'HR_post', ctr.hr(indsposthr), 'BR_post', ctr.br(indspostbr), ...
                                          'setIDcount', mydict(ctr.repset), ...
                                          'polarity', ctr.pol, ...
                                          'setid', ctr.repset);
                    end
                    
                    
                    T.Properties.VariableNames = {'Pre', 'During', 'Post'};
                    T.Properties.RowNames = {'hr','br'};
                    ctr.results = T;
                end
                
            end
        end
        
        function obj = writeresultstable(obj, phase, attr, excelind)
            % Aggregate results of the means pre and during stimulation
            % 'phase' - cell array of combination of 'pre','during','post'
            % values to be extracted
            % 'attr' - specify physiological properties to be extracted as
            % a cell array 'hr','br'
            
            if ~exist('phase','var') || isempty(phase), phase = {'pre','during'}; end
            if ~exist('attr','var') || isempty(attr), attr = {'hr','br'}; end
            allattr = {'HR','BR'}; allphase = {'Pre','During','Post'};
            [~, cmdout] = system('whoami');
            cmdout = strtrim(cmdout);
            switch cmdout
                case 'tnpsoff\tnpla'
                    valphase = {allphase{cell2mat(cellfun(@find, cellfun(@(x,y) strcmpi(x, phase), allphase , 'un',0), 'un',0))}};
                    valattr = {allattr{cell2mat(cellfun(@find, cellfun(@(x) strcmpi(x, attr), allattr, 'un',0), 'un', 0))}};
                    
                otherwise
                    valphase = {allphase{cell2mat(cellfun(@(x) contains_(x, phase, 'IgnoreCase',true), allphase, 'un',0))}};
                    valattr = {allattr{cell2mat(cellfun(@(x) contains_(x, attr, 'IgnoreCase',true), allattr, 'un',0))}};
            end
            
            % Find all trains in all the records
            alltrains = horzcat(obj.trains{:});
            repsets = {alltrains(:).repset};
            uniquereps = setdiff(unique(repsets), '0.0');
            trainidxs = zeros(1,numel(repsets));
            polarities = [alltrains(:).pol]; % polarity of each stim train
            
            % Get indices of trains of the same repset
            for jj=1:numel(repsets)
                if jj == 1
                    trainidxs(jj) = 1;
                else
                    if strcmp(repsets{1,jj}, repsets{1,jj-1})
                        trainidxs(jj) = trainidxs(jj-1) + 1;
                    else
                        trainidxs(jj) = 1;
                    end
                end
            end
            [~, filename] = fileparts(obj.inputdir);
            
            % Create template for table
            T = table();
            %             T = table('Size',[numel(uniquereps) 15], 'VariableTypes', {'double','double','double','double','double','double','double','double',...
            %                 'double','double','double','double','double','double','double'});
            %             T.Properties.VariableNames = {'RatID','RepID','SetID','HR_Pre_pos','HR_During_pos','HR_Post_pos','HR_Pre_neg','HR_During_neg','HR_Post_neg',...
            %                 'BR_Pre_pos','BR_During_pos','BR_Post_pos', 'BR_Pre_neg','BR_During_neg','BR_Post_neg'};
            
            % Run through each set-ID and populate the results table
            T_ = cell(1, numel(uniquereps)); % for combinepolarities
            Ts = cell(1, numel(uniquereps)*2); % for separate polarities
            
            % Choose mode for extracting results
            usefirsttrain = false; % returns only first train for each set-ID
            combinepolarities = false; % One line for each set-id and sequentially report each train
            separatepols = true; % one line for each set-id and polarity. report only trains of particular set-id and polarity
            polcounter = 0;
                    
            for sid = 1:numel(uniquereps)
                
                % ASSIGN RAT - ID. Fix this later.
                T(uniquereps{1,sid},'RatID') = {filename};
                
                % Assign rep ID and set ID
                parserep = strsplit(uniquereps{1,sid}, '.');
                T(uniquereps{1,sid},'RepID') = {str2double(parserep{1})};
                T(uniquereps{1,sid},'SetID') = {str2double(parserep{2})};
                
                % Find trains that have the current set ID
                valididxs = ismember(repsets, uniquereps{1,sid}); % valid indices for a given set-ID
                curpols = polarities(valididxs); % polarities of current trains
                polmap = containers.Map([1 -1],{'Pos'  'Neg'});
                pols = arrayfun(@(x) polmap(x), curpols, 'UniformOutput', false);
                subsettrains = alltrains(valididxs);
                
                % Initial write results function which returned only the
                % first train of each stimulation
                checkpols = unique(curpols);
                if usefirsttrain
                    for p = 1:numel(checkpols)
                        cpol = checkpols(p);
                        subsetpols = [subsettrains.pol];
                        polidx = find(subsetpols == cpol, 1, 'first');
                        valtrain = subsettrains(polidx);
                        switch cpol, case 1, polvar = 'pos'; case -1, polvar = 'neg'; end
                        for jj=1:numel(valattr)
                            for kk=1:numel(valphase)
                                eval(sprintf('T(uniquereps{1,sid}, ''%s_%s_%s'') = valtrain.results(''%s'',''%s'');', ...
                                    valattr{1,jj}, valphase{1,kk}, polvar, attr{1,jj}, valphase{1,kk}));
                            end
                        end
                    end
                end
                
                if combinepolarities
                    temp1 = cell2mat(arrayfun(@(x) table2array(x.results(:, valphase)), subsettrains, 'UniformOutput', false));
                    temp2 = reshape(temp1, 2, 2, []);    % attributes x phase x trains
                    temp3 = permute(temp2, [2 1 3]);     % phase x attributes x trains
                    alltables = reshape(temp3, [], 1)';  % flatten the array
                    
                    varnames = arrayfun(@(itrain) {sprintf('Pre_HR_%s_%d', pols{itrain}, itrain) ...
                                     sprintf('Dur_HR_%s_%d', pols{itrain}, itrain) ...
                                     sprintf('Pre_BR_%s_%d', pols{itrain}, itrain) ...
                                     sprintf('Dur_BR_%s_%d', pols{itrain}, itrain)}, 1:length(pols), 'UniformOutput', false);
                    varnames = [varnames{:}];
                    
                    repid = str2double(parserep{1});
                    setid = str2double(parserep{2});
                    repset = uniquereps(sid);
                    T_{sid} = array2table([0 repid setid alltables], 'VariableNames', [{'RatID', 'RepID' 'SetID'} varnames]);
                    T_{sid}.Properties.RowNames = repset;
                    T_{sid}(repset, 'RatID') = {excelind};                    
                end
                
                if separatepols
                    availablepols = unique(curpols);
                    for p=1:numel(availablepols)
                        polcounter = polcounter + 1;
                        polx = availablepols(p); % current polarity to extract
                        resultspolsidx = [subsettrains.pol] == polx;
                        temp1 = arrayfun(@(x) table2array(x.results(:, valphase)), subsettrains(resultspolsidx), 'UniformOutput', false);
                        temp2 = cellfun(@transpose, temp1, 'un',0);
                        temp3 = cellfun(@(x) x(:), temp2, 'un',0);
                        alltables = vertcat(temp3{:})'; % flatten the array
                        varnames = arrayfun(@(itrain) {sprintf('Pre_HR_%d', itrain) ...
                            sprintf('Dur_HR_%d', itrain) ...
                            sprintf('Pre_BR_%d', itrain) ...
                            sprintf('Dur_BR_%d', itrain)}, 1:sum(resultspolsidx), 'UniformOutput', false);
                        varnames = [varnames{:}];
                        
                        repid = str2double(parserep{1});
                        setid = str2double(parserep{2});
                        repset = uniquereps(sid);
                        Ts{polcounter} = array2table([0 repid setid polx alltables], 'VariableNames', [{'RatID', 'RepID' 'SetID', 'Polarity'} varnames]);
                        Ts{polcounter}.Properties.RowNames = repset;
                        Ts{polcounter}(repset, 'RatID') = {excelind};
                    end
                    
                end
            end
            
            % write final results table to 'hut' object
            if usefirsttrain
                T.Properties.RowNames = uniquereps;
                obj.results = T;
            elseif combinepolarities
                obj.results = T_;
            elseif separatepols
                obj.results = Ts;
            end
        end
        
        function save_results2_file(obj)
            getreps = cellfun(@(y) arrayfun(@(x) x.repset, y, 'un',0), obj.trains, 'UniformOutput', false);
            allrepsets = horzcat(getreps{:});
            validsetids = setdiff(unique(allrepsets), '0.0');
            valididxs = ismember(allrepsets, validsetids);
            alltrains = horzcat(obj.trains{:});
            results_ = arrayfun(@(y) arrayfun(@(x) x.results2, y, 'un',0), alltrains(valididxs), 'UniformOutput', false);
            counter = 0;
            if iscell(results_)
                for jj=1:numel(results_)
                    if ~isempty(results_{jj}{1})
                        counter = counter + 1;
                        results{counter} = results_{jj}{1};
                    end
                end
                
            else
                results = results_;
            end
            results = cell2mat(results);
            [~, name] = fileparts(obj.inputdir);
            filename = sprintf('%s\\%s_results.mat', obj.inputdir, name);
            save(filename, 'results');
        end
        
        function obj = saveplotdata(obj)
            % Store low-level plotting functions as a separate master
            % structure of cell arrays.
            
            
        end
        
        function rawplot(obj, setid, polarity, trainnum)
            % Plot low-level figure for a given set-ID and polarity.
            % trainnum - ordinal number of train in a given sequence of
            % trains of the same set - IDs
            
            if ~exist('trainnum','var') || isempty(trainnum), trainnum = 1; end
            
            % Find all trains in all the records
            alltrains = horzcat(obj.trains{:});
            repsets = {alltrains(:).repset};
            polarities = [alltrains(:).pol]; % polarity of each stim train
            for jj=1:numel(setid)
                for kk =1:numel(polarity)
                    loopidx = find(ismember(repsets, setid{jj}) & ismember(polarities, polarity(kk)));
                    if max(trainnum) <= numel(loopidx)
                        trainidxs{jj,kk} = loopidx(trainnum);
                    else
                        warning('Train number does not exist for set-ID %s. Using highest train number for this set ID instead', setid{jj});
                        trainidxs{jj,kk} = loopidx;
                    end
                end
            end
            
            curtrainidx = [trainidxs{:}];
            subsettrains = {alltrains(curtrainidx)};
            
            for jj=1:numel(subsettrains{1})
                subsettrains{1}(jj).plot_physio();
                suptitle(sprintf('Set - ID: %s, polarity = %d',subsettrains{1}(jj).repset, subsettrains{1}(jj).pol));
            end
        end
        
        function save_params(obj)
            % save parameters used for generating plots and analyses
            ecgmindists = arrayfun(@(iseg) [obj.physdata(iseg).ecg_params.mindist], 1:length(obj.physdata), 'UniformOutput', false);
            breathmindists = arrayfun(@(iseg) [obj.physdata(iseg).breath_params.mindist], 1:length(obj.physdata), 'UniformOutput', false);
            breathminpeakprom = arrayfun(@(iseg) [obj.physdata(iseg).breath_params.minpeakprom], 1:length(obj.physdata), 'UniformOutput', false);
            breathsmoothvar = arrayfun(@(iseg) [obj.physdata(iseg).breath_params.smoothvar], 1:length(obj.physdata), 'UniformOutput', false);
            existflag = any(arrayfun(@isempty, [ecgmindists breathmindists breathminpeakprom breathsmoothvar]));
            zeroflag = any(arrayfun(@(x) isequal(x, 0), [ecgmindists breathmindists breathminpeakprom breathsmoothvar]));
            if existflag || zeroflag
                warning('One or more variables not assigned during call to store');
                return;
            else                
                [~, filename] = fileparts(obj.physdata(1).inputfile);
                filesavename = sprintf('%s\\%s', obj.inputdir, strcat(filename, '_params.mat'));
                ecg_params = {obj.physdata(:).ecg_params};
                breath_params = {obj.physdata(:).breath_params};
                result = obj.results;
                save(filesavename, 'ecg_params', 'breath_params', 'result');                
            end
        end
        
        function load_params(obj, filesavename)
            if ~exist('filesavename', 'var') || ~exist(filesavename, 'file')
                % load save parameters            
                [~, filename] = fileparts(obj.physdata(1).inputfile);
                filesavename = sprintf('%s\\%s', obj.inputdir, strcat(filename, '_params.mat'));
            end
            
            if ~exist(filesavename, 'file')
                warning('File with saved parameters not found.');
                return;
            else
                r1 = load(filesavename);
                for iseg = 1:length(r1.ecg_params)
                    if isempty([r1.ecg_params{iseg}.mindist]) || isempty([r1.breath_params{iseg}.mindist]) || isempty([r1.breath_params{iseg}.minpeakprom])...
                            || isempty([r1.breath_params{iseg}.smoothvar])
                        warning('One or more variables is empty during call to load');
                        return;
                    else
                        disp('Loading existing saved values for physiological signal analysis..')
                        for irecord = 1:obj.physdata(iseg).file_meta.n_records
                            obj.physdata(iseg).ecg_params(irecord).mindist = r1.ecg_params{iseg}(irecord).mindist;
                            obj.physdata(iseg).ecg_params(irecord).minpeakprom = r1.ecg_params{iseg}(irecord).minpeakprom;
                            obj.physdata(iseg).breath_params(irecord).mindist = r1.breath_params{iseg}(irecord).mindist;
                            obj.physdata(iseg).breath_params(irecord).minpeakprom = r1.breath_params{iseg}(irecord).minpeakprom;
                            obj.physdata(iseg).breath_params(irecord).smoothvar = r1.breath_params{iseg}(irecord).smoothvar;
                        end
                    end
                end
            end
            
        end
        
        function save_results_to_excel(obj, overwriteflag, capcolumns)
            % Converts cell array of tables to a combined table and
            % rewrites 'results' property of UTExperiment class with this
            % modified table. 
            % Save results table to excel
%             if ~exist('capcolumns','var') || isempty(capcolumns), capcolumns = true; end
%             if capcolumns
%                 if ~isnumeric(capcolumns)
%                     [~, fname1] = fileparts(obj.physdata(1).inputfile);
%                     filename = sprintf('%s\\%s', obj.inputdir, strcat(fname1, '_results.xlsx'));
%                     
%                     if exist(filename, 'file') == 2
%                         disp('Loading existing results file')
%                         [~, ~, raw] = xlsread(filename, 'Sheet',1);
%                         
%                     end
%                 end
%             end
            
            if ~exist('overwriteflag','var') || isempty(overwriteflag), overwriteflag = false; end
            
            % check if results excel sheet already exists
            [excelpath, excelsuffix] = fileparts(obj.physdata(1).inputfile);
            checkvar = [excelsuffix, '_results.xlsx'];
            checkflag = exist(sprintf('%s\\%s', excelpath, checkvar), 'file');
            
            if checkflag
                disp('Results file already exists\n')
                
                if overwriteflag
                    disp('Deleting old excel sheet and writing new file')
                    delete(sprintf('%s\\%s', excelpath, checkvar));
                else
                    disp('File already exists. Set ''overwriteflag'' to ''true''')
                    return;
                end
            end
            
            % write table to excel
            if isempty(obj.results)
                warning('Results table empty. Exiting this function..');
                return;
            else                
                [~, fname1] = fileparts(obj.physdata(1).inputfile);
                filename = sprintf('%s\\%s', obj.inputdir, strcat(fname1, '_results.xlsx'));
                if istable(obj.results)
                    writetable(obj.results, filename, 'Sheet',1);                
                elseif iscell(obj.results)
%                     currow = 1;
%                     for itable = 1:length(obj.results)
%                         numcols = size(obj.results{itable}, 2);
%                         myrange = sprintf('A%d:%s%d', currow, xlscol(numcols), currow+1);
%                         writetable(obj.results{itable}, filename, 'Sheet', 1, 'Range', myrange);
%                         currow = currow + 3;
%                     end
                    
                    % get size of rows
                    nonemptyresults = find(~cellfun(@isempty, obj.results));
                    rowsizes = cellfun(@width, obj.results(nonemptyresults));
                    [maxsize, maxind] = max(rowsizes);
                    
                    % header info for final table
                    maxheader = obj.results{maxind}.Properties.VariableNames;
                    
                    % replicate smaller length rows with NaN's
                    for r=1:numel(nonemptyresults)
                        
                        newtablevals(r,:) = [table2array(obj.results{nonemptyresults(r)}) NaN(1, maxsize - width(obj.results{nonemptyresults(r)}))];
                            
                    end
                    
                    newtable = array2table(newtablevals);
                    newtable.Properties.VariableNames = maxheader;
                    
                    % reset results property of hut object with modified
                    % table
                    obj.results = newtable;
                    writetable(obj.results, filename, 'Sheet',1); 
                    
                else
                    error('unexpected results format\n');
                end
            end
        end
        %% Functions for UT data modeling and analyses
        % Subash  Padmanaban 07/16/2018
        
        function obj = stackalltrains(obj, chanfcn)
            % Stack all neural responses to start at stimulus triggered
            % points
            numperiods = numel(obj.trains);
            for jj=1:numperiods
                numtrains = numel(obj.trains{jj});
                
                for kk = 1:numtrains
                    clc, sprintf('Stacking pulses for record %d, train %d',jj, kk)
                    obj.trains{jj}(kk).stack_ep(chanfcn);
                end
            end
        end
        
        function obj = getneuralfeatures(obj, metric)
            % This function takes the stimulus triggered neural response of
            % the vagus nerve signal and computes neural features that
            % correspond to alpha, beta and gamma neural activity.
            % metric - 'auc': area under the curve (default)
            if ~exist('metric','var') || isempty(metric), metric = 'auc'; end
            cvelocities = [0.001 3; 3 35; 35 80].*1000; % Conduction velocities
            
            switch metric
                case 'auc'
                    obj.feature_metric = 'auc';
                    numperiods = numel(obj.trains);
                    
                    % For each time period
                    for jj=1:numperiods
                        numtrains = numel(obj.trains{jj});
                        
                        % For each train
                        for kk=1:numtrains
                            if ~strcmp(obj.trains{jj}(kk).repset, '0.0')
                                rectifiedsig = abs(obj.trains{jj}(kk).evokedpotentials);
                                tvec = linspace(0, size(rectifiedsig,1)/obj.sr, size(rectifiedsig,1));
                                timelims = cell2mat(arrayfun(@(x) obj.distance./x, cvelocities, 'un',0));
                                neural_idxs = arrayfun(@(x,y) tvec > x & tvec < y, timelims(:,2),timelims(:,1), 'un',0);
                                
                                % For each pulse
                                obj.trains{jj}(kk).nfeats = [];
                                for mm=1:size(rectifiedsig,2)
                                    clc, sprintf('Calculating neural features for record %d, train %d, pulse %d',jj,kk,mm)
                                    obj.trains{jj}(kk).nfeats(:,mm) = cell2mat(cellfun(@(x) trapz(tvec(x), ...
                                        rectifiedsig(x, mm)), neural_idxs, 'un',0));
                                end
                            end
                        end
                    end
            end
        end
        
        
        function obj = plotneuralfeatures(obj, setid, trainnum, polarity, plotvar)
            % Plot extracted neural features
            %             if ~exist('record','var') || isempty(record), record = 1:numel(obj.trains); end
            %             repsets = {obj.trains{record}(:).repset}; % repset corresponding to each stim train
            alltrains = horzcat(obj.trains{:});
            repsets = {alltrains(:).repset};
            uniquereps = setdiff(unique(repsets), '0.0');
            trainidxs = zeros(1,numel(repsets));
            polarities = [alltrains(:).pol]; % polarity of each stim train
            
            % Get indices of trains of the same repset
            for jj=1:numel(repsets)
                if jj == 1
                    trainidxs(jj) = 1;
                else
                    if strcmp(repsets{1,jj}, repsets{1,jj-1})
                        trainidxs(jj) = trainidxs(jj-1) + 1;
                    else
                        trainidxs(jj) = 1;
                    end
                end
            end
            
            % True if variable exists
            setidflag = exist('setid','var') && ~isempty(setid);
            trainnumflag = exist('trainnum','var') && ~isempty(trainnum);
            polarityflag = exist('polarity','var') && ~isempty(polarity);
            
            % If user does not provide setid, trainnum and polarities, set
            % defaults
            if ~setidflag, c_setids = uniquereps; else, c_setids = setid; end
            if ~trainnumflag, c_trains = trainidxs; else, c_trains = trainnum; end
            if ~polarityflag, c_pol = [1 -1]; else, c_pol = polarity; end
            
            valididxs = ismember(repsets, c_setids) & ismember(trainidxs, c_trains) & ismember(polarities, c_pol);
            
            % Subset of trains to be plotted as requested by the user
            subsettrains = {alltrains(valididxs)};
            
            % Create colors for plotting
            trainlegs = arrayfun(@num2str, c_trains(valididxs), 'un',0);
            pollegs = cellfun(@num2str, {subsettrains{1}(:).pol}, 'un',0);
            setlegs = {subsettrains{1}(:). repset};
            plot_legends = cellfun(@(x,y,z) [x, ',', y, ',', z], setlegs, trainlegs, pollegs, 'un',0);
            plotrepsets = {subsettrains{1}(:). repset};
            plotcols = hsv(numel(plot_legends));
            
            %             figure;
            
            switch plotvar
                case 'individual'
                    for jj=1:numel(plot_legends)
                        plotdata = [subsettrains{1}(jj).nfeats];
                        scatter3(plotdata(1,:), plotdata(2,:), plotdata(3,:), 500, plotcols(jj,:), 'Marker','.');
                        hold on;
                    end
                case 'group'
                    plotidxs = find(ismember(plotrepsets, plot_legends{1,jj}(1:3)));
                    plotdata = [subsettrains{1}(plotidxs).nfeats];
                    scatter3(plotdata(1,:), plotdata(2,:), plotdata(3,:), 500, plotcols(jj,:), 'Marker','.');
                    hold on;
            end
            
            xlabel('Gamma (0 - 3 m/s)')
            ylabel('A - delta (3 - 35 m/s)')
            zlabel('A - beta (35 - 80 m/s)')
            legend(plot_legends)
            title(sprintf('Pulse triggered neural response, %s',obj.feature_metric))
            set(gca, 'FontSize',20)
            
        end
        
        function obj = compareclusters(obj, clustervals, polarity, specialval)
            % This function computes the inter and intra class distance
            % between two clusters and computes the cluster similarity
            
            
            % get data for both the clusters
            alltrains = horzcat(obj.trains{:});
            repsets = {alltrains(:).repset};
            polarities = [alltrains(:).pol];
            if ~exist('specialval','var') || isempty(specialval), specialval = 'none'; end
            switch specialval
                case 'none'
                    for jj=1:numel(polarity)
                        dataidx1 = find(ismember(repsets, clustervals{1}) & ismember(polarities, polarity(jj)));
                        data1 = [alltrains(dataidx1).nfeats];
                        dataidx2 = find(ismember(repsets, clustervals{2}) & ismember(polarities, polarity(jj)));
                        data2 = [alltrains(dataidx2).nfeats];
                        
                        % compute Fisher's linear discriminant statistic
                        % compute covariance matrices for the data
                        S(jj) = UTExperiment.ldastat(data1, data2);
                    end
                    obj.modeling.indcompare.cluster1 = clustervals{1};
                    obj.modeling.indcompare.cluster2 = clustervals{2};
                    obj.modeling.indcompare.polarities = polarity;
                    obj.modeling.indcompare.Fisherstat = S;
                    
                case 'samepairs'
                    % Find same pairs of stimulation
                    uniquereps = setdiff(unique(repsets), '0.0');
                    oneidxs = find(cell2mat(cellfun(@(x) startsWith(x, '1.'), uniquereps, 'un',0)));
                    startswith1 = {uniquereps{oneidxs}};
                    twoidxs = find(cell2mat(cellfun(@(x) startsWith(x, '2.'), uniquereps, 'un',0)));
                    startswith2 = {uniquereps{twoidxs}};
                    for jj=1:numel(startswith1)
                        temp1 = strsplit(startswith1{jj}, '.');
                        setvals(jj) = str2double(temp1{2});
                    end
                    
                    for jj=1:numel(startswith2)
                        temp2 = strsplit(startswith2{jj}, '.');
                        setvals2(jj) = str2double(temp2{2});
                    end
                    
                    pairvals = sort(setvals(ismember(setvals, setvals2)));
                    for jj=1:numel(pairvals)
                        clusterpairs{1,jj} = ['1.' num2str(pairvals(jj))];
                        clusterpairs{2,jj} = ['2.' num2str(pairvals(jj))];
                        dataidx1 = find(ismember(repsets, clusterpairs{1,jj}));
                        data1 = [alltrains(dataidx1).nfeats];
                        dataidx2 = find(ismember(repsets, clusterpairs{2, jj}));
                        data2 = [alltrains(dataidx2).nfeats];
                        
                        % compute Fisher's linear discriminant statistic
                        % compute covariance matrices for the data
                        sigma1 = cov(data1'); sigma2 = cov(data2');
                        sigma_avg = (sigma1 + sigma2)./2;
                        
                        % decision criteria
                        mu1 = mean(data1, 2);
                        mu2 = mean(data2, 2);
                        w_ = inv(sigma_avg)*(mu1 - mu2);
                        
                        % Fisher's statistic
                        numerator = (dot(w_, (mu1 - mu2)))^2; % between-cluster variance
                        denominator = w_'*(sigma1 + sigma2)*w_; % within-cluster variance
                        S(jj) = numerator/denominator;
                    end
                    T = table({clusterpairs{1,:}}', {clusterpairs{2,:}}', S', 'VariableNames',{'SetID1','SetID2','LDAStatistic'});
                    obj.modeling.samepairs = T;
                    
                case 'samepairspol'
                    % Find same pairs of stimulation
                    uniquereps = setdiff(unique(repsets), '0.0');
                    oneidxs = find(cell2mat(cellfun(@(x) startsWith(x, '1.'), uniquereps, 'un',0)));
                    startswith1 = {uniquereps{oneidxs}};
                    twoidxs = find(cell2mat(cellfun(@(x) startsWith(x, '2.'), uniquereps, 'un',0)));
                    startswith2 = {uniquereps{twoidxs}};
                    for jj=1:numel(startswith1)
                        temp1 = strsplit(startswith1{jj}, '.');
                        setvals(jj) = str2double(temp1{2});
                    end
                    
                    for jj=1:numel(startswith2)
                        temp2 = strsplit(startswith2{jj}, '.');
                        setvals2(jj) = str2double(temp2{2});
                    end
                    
                    pairvals = sort(setvals(ismember(setvals, setvals2)));
                    counter = 0;
                    for jj=1:numel(pairvals)
                        counter = counter +1;
                        clusterpairs{1,counter} = ['1.' num2str(pairvals(jj))];
                        clusterpairs{2,counter} = ['2.' num2str(pairvals(jj))];
                        
                        % Compare clusters for +ve polarity
                        polidx(counter) = 1;
                        dataidx1 = find(ismember(repsets, clusterpairs{1,counter}) & ismember(polarities, 1));
                        data1 = [alltrains(dataidx1).nfeats];
                        dataidx2 = find(ismember(repsets, clusterpairs{2, counter}) & ismember(polarities, 1));
                        data2 = [alltrains(dataidx2).nfeats];
                        
                        % compute Fisher's linear discriminant statistic
                        % compute covariance matrices for the data
                        S(counter) = UTExperiment.ldastat(data1, data2);
                        
                        % Compare clusters for -ve polarity
                        counter = counter +1;
                        clusterpairs{1,counter} = ['1.' num2str(pairvals(jj))];
                        clusterpairs{2,counter} = ['2.' num2str(pairvals(jj))];
                        polidx(counter) = -1;
                        dataidx1 = find(ismember(repsets, clusterpairs{1,counter}) & ismember(polarities, -1));
                        data1 = [alltrains(dataidx1).nfeats];
                        dataidx2 = find(ismember(repsets, clusterpairs{2, counter}) & ismember(polarities, -1));
                        data2 = [alltrains(dataidx2).nfeats];
                        
                        % Fisher's LDA statistics
                        S(counter) = UTExperiment.ldastat(data1, data2);
                    end
                    T = table({clusterpairs{1,:}}', {clusterpairs{2,:}}', polidx', S', 'VariableNames',{'SetID1','SetID2','Polarity','LDAStatistic'});
                    obj.modeling.samepairspol = T;
            end
        end
        
        function obj = stimfeatures(obj)
            % Extract input features for modeling electrical
            % stimulation parameters to evoked neural potentials
            
            numperiods = numel(obj.trains);
            % For each time period
            for jj=1:numperiods
                numtrains = numel(obj.trains{1});
                
                % For each train
                for kk=1:numtrains
                    curtrain = obj.trains{jj}(kk);
                    
                    % only process trains that have a specific stimulation
                    % parameter mapped
                    if ~strcmp(curtrain.repset, '0.0')
                        clc, sprintf('Extracting stimulation parameters for period %d, train %d', jj,kk)
                        paramset = obj.wparams.params(curtrain.repset);
                        numpulses = size(curtrain.nfeats, 2);
                        pwfeature = paramset.pulsewidth_us.*ones(1, numpulses);
                        freqfeature = paramset.freq_Hz.*ones(1, numpulses);
                        pulsenum = 1:numpulses;
                        threshfeat = paramset.amp.*ones(1,numpulses);
                        polfeat = obj.trains{jj}(kk).pol.*ones(1,numpulses);
                        obj.trains{jj}(kk).stimfeats = [pwfeature; freqfeature; pulsenum; threshfeat; polfeat];
                    end
                end
            end
            
        end
        
        function obj = gettraintest(obj, biosignal, splitvar, testvar)
            % Use this function to split neural and stimulation parameter
            % features into training and testing sets for modeling.
            % splitvar: 'holdout' (MUST be called with testing indices input 'testvar')
            % '1vs2': splits all stimulation trials with 1.x into training
            % and 2.x into testing
            % '+vs-' : splits stimulation trials into +ve and -ve
            
            % get data
            alltrains = horzcat(obj.trains{:});
            repsets = {alltrains(:).repset};
            uniquereps = setdiff(unique(repsets), '0.0');
            if ~exist('biosignal','var') || isempty(biosignal), ...
                    error('Call this function by specifying whether you want to use neural or physiological signal as input'); end
            
            switch splitvar
                case 'holdout'
                    if ~exist('testvar','var') || isempty(testvar), disp('When calling a holdout split, testvar should be specified'); return; end
                    if sum(ismember(uniquereps, testvar)) ~= numel(unique(testvar)), ...
                            disp('Some testing repsets may not be available for this session'); return; end
                    trainingsets = setdiff(uniquereps, testvar);
                    traininds = ismember(repsets, trainingsets);
                    testinds = ismember(repsets, testvar);
                    
                    switch biosignal
                        case 'neural'
                            obj.modeling.trainX = [obj.trains{:}(traininds).stimfeats];
                            obj.modeling.trainY = [obj.trains{:}(traininds).nfeats];
                            obj.modeling.testX = [obj.trains{:}(testinds).stimfeats];
                            obj.modeling.testY = [obj.trains{:}(testinds).nfeats];
                            obj.modeling.trainsets = trainingsets;
                            
                        case 'physio'
                            alltrains = horzcat(obj.trains{:});
                            subsettrains = {alltrains(traininds)};
                            features = [subsettrains{:}.physfeats];
                            obj.modeling.trainX = [features(:).phys];
                            obj.modeling.trainY = [features(:).stim];
                            
                            subsettrains = {alltrains(testinds)};
                            features = [subsettrains{:}.physfeats];
                            obj.modeling.testX = [features(:).phys];
                            obj.modeling.testY = [features(:).stim];
                            obj.modeling.trainsets = trainingsets;
                    end
                    obj.modeling.testsets = testvar;
                    obj.modeling.method = splitvar;
                    
                case '1vs2'
                    % In this context, testvar is used as a logical to
                    % assign repetitions 2.x as the test set
                    if ~exist('testvar','var') || isempty(testvar), testvar = true; end
                    oneinds = cell2mat(cellfun(@(x) contains(x, '1.'), repsets, 'un',0));
                    twoinds = cell2mat(cellfun(@(x) contains(x, '2.'), repsets, 'un',0));
                    if testvar
                        testinds = twoinds; traininds = oneinds;
                    else
                        testinds = oneinds; traininds = twoinds;
                    end
                    
                    switch biosignal
                        case 'neural'
                            obj.modeling.trainX = [obj.trains{:}(traininds).stimfeats];
                            obj.modeling.trainY = [obj.trains{:}(traininds).nfeats];
                            obj.modeling.testX = [obj.trains{:}(testinds).stimfeats];
                            obj.modeling.testY = [obj.trains{:}(testinds).nfeats];
                            
                        case 'physio'
                            alltrains = horzcat(obj.trains{:});
                            subsettrains = {alltrains(traininds)};
                            features = [subsettrains{:}.physfeats];
                            obj.modeling.trainY = [features(:).phys];
                            obj.modeling.trainX = [features(:).stim];
                            
                            subsettrains = {alltrains(testinds)};
                            features = [subsettrains{:}.physfeats];
                            obj.modeling.testY = [features(:).phys];
                            obj.modeling.testX = [features(:).stim];
                            
                    end
                    obj.modeling.trainsets = unique(repsets(traininds));
                    obj.modeling.testsets = unique(repsets(testinds));
                    obj.modeling.method = splitvar;
                    
                case '+vs-'
                    
                    
            end
        end
        
        function plotmodelfeatures(obj, plotset, stimparam, neuralparam)
            if ~exist('plotset','var') || isempty(plotset), plotset = 'train'; end
            if ~exist('stimparam','var') || isempty(stimparam), stimparam = 'amp'; end
            if ~exist('neuralparam','var') || isempty(neuralparam), neuralparam = 'ABeta'; end
            
            allsets = {'train'; 'test'};
            allx = {obj.modeling.trainX; obj.modeling.testX};
            ally = {obj.modeling.trainY; obj.modeling.testY};
            plotinds = cell2mat(cellfun(@(x) contains(x, plotset), allsets, 'un',0));
            plotx = horzcat(allx{plotinds}); ploty = horzcat(ally{plotinds});
            
            allparams = {'pw'; 'freq'; 'pnum'; 'amp'};
            paraminds = find(cell2mat(cellfun(@(x) contains(x, stimparam), allparams, 'un',0)));
            
            allneuralp = {'ADelta'; 'ABeta'; 'Gamma'};
            neuralinds = find(cell2mat(cellfun(@(x) contains(x, neuralparam), allneuralp, 'un',0)));
            
            for nfeat=1:numel(neuralinds)
                figure;
                
                switch numel(paraminds)
                    case 1
                        plot(plotx(paraminds,:), ploty(neuralinds(nfeat),:), '.')
                        xlabel(allparams{paraminds})
                        ylabel(allneuralp{neuralinds(nfeat)})
                        title(sprintf('AUC of %s as a function of %s', allneuralp{neuralinds(nfeat)}, allparams{paraminds}));
                        axis('tight');
                        set(gca, 'FontSize',20)
                        
                    case 2
                        plot3(plotx(paraminds(1),:), plotx(paraminds(2),:), ploty(neuralinds(nfeat),:), '.')
                        xlabel(allparams{paraminds(1)})
                        ylabel(allparams{paraminds(2)})
                        zlabel(allneuralp{neuralinds(nfeat)})
                        title(sprintf('AUC of %s as a function of %s and %s', allneuralp{neuralinds(nfeat)}, allparams{paraminds(1)}, allparams{paraminds(2)}));
                        axis('tight');
                        set(gca, 'FontSize',20)
                        
                end
            end
        end
        
        function obj = extract_physio_features(obj)
            % extract physiological and stimulation features based on the
            % number of stimulation pulses delivered
            
            traincount = 0; % variable to count the count of stim train delivered to the animal
            % For each time period
            for jj=1:numel(obj.trains)
                numtrains = numel(obj.trains{jj});
                
                % For each train
                for kk=1:numtrains
                    curtrain = obj.trains{jj}(kk);
                    
                    % only process trains that have a specific stimulation
                    % parameter mapped
                    if ~strcmp(curtrain.repset, '0.0')
                        traincount = traincount + 1;
                        clc, sprintf('Extracting physiological and stimulation features for period %d, train %d',jj, kk)
                        
                        % initialize physresults structure
                        features = struct('phys',[], 'stim',[]); 
                        
                        % Obtain stimulation parameters for set-ID
                        stimparams = curtrain.params; 
                        
                        % Get pulse times
                        pulsedur = 1/stimparams.freq_Hz;
                        pulsetimes = curtrain.tecg(1) + curtrain.buffertime: pulsedur: curtrain.tecg(end) - curtrain.buffertime;
                        
                        % Interpolate heart rate and breathing rate values
                        % for the extracted pulse times
                        
                        int_hr = interp1(curtrain.thr, curtrain.hr, pulsetimes);
                        
                        int_br = interp1(curtrain.tbr, curtrain.br, pulsetimes);
                        
                        numpulses = numel(pulsetimes);
                        features.phys = [int_hr; int_br; ...
                            table2array(curtrain.results('hr','Pre'))*ones(1,numpulses);...
                            table2array(curtrain.results('br','Pre'))*ones(1,numpulses)];
                        
                        % extract stimulation parameters
                        features.stim = [stimparams.amp*ones(1,numpulses); stimparams.freq_Hz*ones(1,numpulses);...
                            stimparams.duration_s*ones(1,numpulses); 1:numpulses; traincount*ones(1, numpulses); ...
                            curtrain.pol*ones(1,numpulses); stimparams.pulsewidth_us*ones(1,numpulses)];
                        curtrain.physfeats = features;
                    end
                end
            end
        end
        
        function add_derivative_channel(obj, name, fc)
            if ~exist('fc', 'var') || isempty(fc)
                fc = 50;
            end
            
            for iseg = 1:length(obj.physdata)
                [dx, dt, ind] = obj.physdata(iseg).derivative_chan(name, fc); 
                if isempty(dx)
                    continue
                end
                units = cellfun(@(x) [x '/sec'], obj.physdata(iseg).channel_meta(ind).units, 'UniformOutput', false);
                obj.physdata(iseg).add_channel(dx, dt, sprintf('d%s/dt', name), units);
            end
        end
        
        function resample_channel(obj, name, sr)
            for iseg = 1:length(obj.physdata)
                obj.physdata(iseg).resample_chan(name, sr);
            end
        end
        
        function plot_chan(obj, name, ax)
            data = [];
            t = [];
            t0 = obj.physdata(1).record_meta(1).record_start;
            for iseg = 1:length(obj.physdata)
                ichan = find(ismember({obj.physdata(iseg).channel_meta.name}, name));
                assert(length(ichan) == 1);
                
                for irecord = 1:obj.physdata(iseg).file_meta.n_records
                    datastr = sprintf('data__chan_%d_rec_%d', ichan, irecord);
                    data = [data; obj.physdata(iseg).(datastr)];
                    nsamps = obj.physdata(iseg).channel_meta(ichan).n_samples(irecord);
                    dt = obj.physdata(iseg).channel_meta(ichan).dt(irecord);
                    recordStart = obj.physdata(iseg).record_meta(irecord).record_start;
                    t = [t; (0:nsamps-1)' * dt + (recordStart - t0) * 86400];
                end
            end
            
            if ~exist('ax', 'var') || ~isvalid(ax)
                figure;ax = gca;
            end
            
            holdState = ishold(ax);
            if ~holdState
                hold(ax, 'on');
            end
            axes(ax);
            reduce_plot(t, data);
            if ~holdState  % only turn hold off if it was off to begin with
                hold(ax, 'off');
            end
        end
        
        function [comments, ia, trains, B] = associate_comments_with_trains(obj)
            % associate each train with a set of comments
            % each group of comments that precedes a set of trains is associated
            
            % comment times (sec)
            tc = [];
            comments = [];
            for iseg = 1:length(obj.physdata)
                tc = [tc (arrayfun(@(x) obj.physdata(iseg).comment_time(x), 1:length(obj.physdata(iseg).comments)) - obj.segtimes(1)) * 86400];
                comments = [comments {obj.physdata(iseg).comments.str}];
            end
            
            % train times (sec)
            tt = [];
            trains = [];
            for iseg = 1:length(obj.trains)
                tt = [tt obj.trains{iseg}.t0];
                trains = [trains obj.trains{iseg}];  % create a linear array of PulseTrain objects (handle class)
            end
            
            % compute the distance matrix
            if verLessThan('matlab', '9.1')
                D = bsxfun(@minus, reshape(tt, 1, []), reshape(tc, [], 1));
            else
                D = reshape(tt, 1, []) - reshape(tc, [], 1);
            end
            D(D < -5) = NaN;
            % D = pdist2(reshape(tc, [], 1), reshape(tt, [], 1));
            [~, mi] = min(D, [], 2);
            [B, ia] = unique(mi, 'last');            
            mi = [mi;length(trains)+1];
            
            % clear the comments property
            for itrain = 1:length(trains)
                trains(itrain).comments = {};
            end
            
            % associte the comments
            for icomment = 1:length(mi)-1
                ind = find(mi > mi(icomment), 1, 'first');
                for itrain = mi(icomment):mi(ind)-1
                    trains(itrain).comments = [trains(itrain).comments comments(icomment)];
                end
            end
        end
        
        function store_phys_in_trains(obj, trains)
            % copy the physiological data into each train
            if ~exist('trains', 'var')
                trains = [obj.trains{:}];
            end
            
            % loop over each train (iseg and itrain)
            for itrain = 1:length(trains)
                % store a handle to the current train
                curtrain = trains(itrain);

                % compute the indices into the physiological data for the current train
                T = (size(curtrain.data, 1) - curtrain.n1 - curtrain.n2) / curtrain.sr;
                adiInds = obj.get_samples_from_time(curtrain.t0 + curtrain.n1 / curtrain.sr + [0 T], 'sec', 'adi');

                % initialize the PulseTrain properties that store physiological data
                numchans = size(adiInds(1).chaninds, 1);                    
                [curtrain.physioNames, curtrain.physioData] = deal(cell(1, numchans));
                curtrain.physioSr = zeros(1, numchans);

                if length(adiInds) ~= 1
                    fprintf('Expected a single data segment for assigning physiological data to trains and got %d\n', length(adiInds));
                end

                % loop over the physiological data segments that contain the data for the train
                for iadiseg = 1:length(adiInds)
                    % get the ADI segment and record
                    seg = adiInds(iadiseg).seg;
                    rec = adiInds(iadiseg).record;
                    % loop over each data channel
                    for ichan = 1:size(adiInds(iadiseg).chaninds, 1)                            
                        % assign the channel name and sampling rate only once
                        if iadiseg == 1
                            curtrain.physioNames{ichan} = obj.physdata(seg).channel_meta(ichan).name;
                            curtrain.physioSr(ichan) = 1 / obj.physdata(seg).channel_meta(ichan).dt(rec);
                        end
                        % assign the data (and concatenate if necessary)
                        chanstr = sprintf('data__chan_%d_rec_%d', ichan, rec);
                        if isnan(adiInds(iadiseg).chaninds(ichan, 1)) || isnan(adiInds(iadiseg).chaninds(ichan, 2))
                            curtrain.physioData{ichan} = [curtrain.physioData{ichan}; {[]}];
                        else
                            inds = adiInds(iadiseg).chaninds(ichan, 1):adiInds(iadiseg).chaninds(ichan, 2);
                            curtrain.physioData{ichan} = [curtrain.physioData{ichan}; obj.physdata(seg).(chanstr)(inds)];
                        end
                    end
                end                    
            end
        end
        
        function [mytable, results] = create_pressure_table(obj, subjectid_, Twindow)
            % create a table with comment labels along the rows and pressure metrics in the columns
            % also return an array of results structures            
            if ~exist('Twindow', 'var') || isempty(Twindow)
                Twindow = 4.9;
            end            
            
            % these are the channel I will operate on
            chanNames = {'PAP', 'RVP', 'dRVP/dt', 'RAP', 'SAP', 'LVP', 'dLVP/dt', 'LAP'};
            
            % I computed the results for getting the row labels in associate_comments_with_trains so I returned them 
            obj.associate_comments_with_trains();  % [comments, ia, trains, B] = 
            % B = [B; length(trains)+1];
            
            % create an empty table
            [comment, pulsewidth_us, freq_Hz, numpulses, duration_s, amp, iti_s, wait_min, ...
                SubjectID, RepID, SetID, Polarity, HR, RR, ...
                sPAP_right, dPAP_right, mPAP_right, ...
                sRVP_right, dRVP_right, mRVP_right, ...
                dRVPdtmax_right, dRVPdtmin_right, ...
                sRAP_right, dRAP_right, mRAP_right, ...
                sSAP_left, dSAP_left, mSAP_left, ...
                sLVP_left, dLVP_left, mLVP_left, ...
                dLVPdtmax_left, dLVPdtmin_left, ...
                sLAP_left, dLAP_left, mLAP_left] = deal([]);
            mytable = table(comment, pulsewidth_us, freq_Hz, numpulses, duration_s, amp, iti_s, wait_min, ...
                SubjectID, RepID, SetID, Polarity, HR, RR, ...
                sPAP_right, dPAP_right, mPAP_right, ...
                sRVP_right, dRVP_right, mRVP_right, ...
                dRVPdtmax_right, dRVPdtmin_right, ...
                sRAP_right, dRAP_right, mRAP_right, ...
                sSAP_left, dSAP_left, mSAP_left, ...
                sLVP_left, dLVP_left, mLVP_left, ...
                dLVPdtmax_left, dLVPdtmin_left, ...
                sLAP_left, dLAP_left, mLAP_left);
            
            SubjectID = subjectid_;
            
            % initialize the figure to show the results of findpeaks
            h = figure;
            ax = cell(1, length(chanNames));
            for ichan = 1:length(chanNames)
                ax{ichan} = subplot(2, 4, ichan);
            end
            
            % create PulseTrain objects that capture the baseline periods            
            [preTrains, postTrains] = obj.create_pre_post_trains(Twindow);
            originalTrains = [obj.trains{:}];
            
            % sort the baseline trains with the original trains
            allTrains = reshape([preTrains; originalTrains; postTrains], 1, []);            
            % allTrains = obj.sort_trains([baselineTrains postTrains], [obj.trains{:}]);
            
            ntrains = length([obj.trains{:}]);  % number of original trains
            obj.extractphysdata(allTrains, 0);  % stores data to the trains
            
            % initialize results
            pdp = {'pre', 'during', 'post'};
            temp = struct('setid', [], 'polarity', [], 'setIDcount', [], 'pulsewidth_us', [], 'freq_Hz', [], 'numpulses', [], 'duration_s', [], 'amp', [], 'iti_s', [], 'wait_min', []);
            excludeList = {'comment', 'pulsewidth_us', 'freq_Hz', 'numpulses', 'duration_s', 'amp', 'iti_s', 'wait_min', 'SubjectID', 'RepID', 'SetID', 'Polarity', ...
                'mPAP_right', 'mRVP_right', 'mRAP_right', 'mSAP_left', 'mLVP_left', 'mLAP_left'};
            for ipdp = 1:length(pdp)
                for ivarname = 1:length(mytable.Properties.VariableNames)
                    if ~ismember(mytable.Properties.VariableNames{ivarname}, excludeList)
                        fn = sprintf('%s_%s', mytable.Properties.VariableNames{ivarname}, pdp{ipdp});
                        temp.(fn) = [];
                    end
                end
            end            
            results = repmat(temp, [1 ntrains]);             
            setIDcount = containers.Map();
            trainct = 1;
            
            % loop over each row in the table
            for itrain = 1:length(allTrains)
                curtrain = allTrains(itrain);
                rowstr = curtrain.comments{end};
                
                % --- fill in the initial part of the results structure ---
                hasValidSetID = ~isempty(obj.wparams) && ismember(curtrain.repset, obj.wparams.params.keys());
                if hasValidSetID
                    results(trainct).setid = curtrain.repset;
                    results(trainct).polarity = curtrain.pol;
                    wp = obj.wparams.params(curtrain.repset);
                    fnwp = fieldnames(wp);
                    for ifnwp = 1:length(fnwp)
                        results(trainct).(fnwp{ifnwp}) = wp.(fnwp{ifnwp});
                    end                
                end                                
                if isequal(curtrain.repset, 'pre')
                    curphase = 'pre';
                elseif isequal(curtrain.repset, 'post')
                    curphase = 'post';
                else
                    curphase = 'during';
                    if ismember(curtrain.repset, setIDcount.keys())
                        setIDcount(curtrain.repset) = setIDcount(curtrain.repset) + 1;
                        results(trainct).setIDcount = setIDcount(curtrain.repset);                        
                    else
                        setIDcount(curtrain.repset) = 1;
                        results(trainct).setIDcount = setIDcount(curtrain.repset);
                    end
                end
                results(trainct).(sprintf('HR_%s', curphase)) = curtrain.hr;
                results(trainct).(sprintf('RR_%s', curphase)) = curtrain.br;  % respiratory rate (breathing rate)
                % ---------------------------------------------------------

                % get the waveform parameters
                if ~isempty(obj.wparams) && ismember(curtrain.repset, obj.wparams.params.keys())
                    params = obj.wparams.params(curtrain.repset);
                    temp = struct2cell(params);
                    [pulsewidth_us, freq_Hz, numpulses, duration_s, amp, iti_s, wait_min] = deal(temp{:});
                else
                    [pulsewidth_us, freq_Hz, numpulses, duration_s, amp, iti_s, wait_min] = deal(NaN);
                end
                
                % get the set id
                if contains_(curtrain.repset, '.')
                    temp = strsplit(curtrain.repset, '.');
                    temp = cellfun(@str2double, temp, 'UniformOutput', false);
                    [RepID, SetID] = temp{:};
                else
                    [RepID, SetID] = deal(NaN);
                end
                % get the polarity
                Polarity = curtrain.pol;
                
                % compute the metrics 
                myMetrics = [pulsewidth_us, freq_Hz, numpulses, duration_s, amp, iti_s, wait_min, ...
                    SubjectID, RepID, SetID, Polarity, mean(curtrain.hr), mean(curtrain.br)];
                for ichan = 1:length(chanNames)
                    % get the physiological data from the train
                    chanInd = find(ismember(curtrain.physioNames, chanNames{ichan}));                        
                    % assert(length(chanInd) == 1, '%s isn''t a channel\n', chanNames{ichan});
                    
                    if isempty(chanInd)
                        % There may only be one of either LAP or RAP but not both
                        pks_d = [];
                        pks_s = [];
                        data = [];                        
                        title(ax{ichan}, chanNames{ichan});
                    else
                        data = curtrain.physioData{chanInd};
                        sr_ = curtrain.physioSr(chanInd);
                        
                        % find the peaks
                        [pks_s, locs_s] = findpeaks(data, 'MinPeakDistance', .25 * sr_, ...
                            'MinPeakProminence', 0.5*range(data));
                        [pks_d, locs_d] = findpeaks(-data, 'MinPeakDistance', .25 * sr_, ...
                            'MinPeakProminence', 0.5*range(data));
                        pks_d = -pks_d;
                        
                        % handle when no peaks are found
                        hasPeaks = true;
                        if isempty(pks_s)
                            pks_s = NaN;                            
                            hasPeaks = false;
                        end
                        if isempty(pks_d)
                            pks_d = NaN;                            
                            hasPeaks = false;
                        end
                        
                        if hasPeaks
                            % plot the result to verify the findpeaks parameters were set appropriately
                            t = (0:length(data)-1) / sr_;
                            plot(ax{ichan}, t, data, '-b', t(locs_s), pks_s, 'r.', t(locs_d), pks_d, 'g.');
                            title(ax{ichan}, chanNames{ichan});
                        end
                    end

                    % assign the metrics
                    if ismember(chanNames{ichan}, {'PAP', 'RVP', 'SAP', 'LVP', 'RAP', 'LAP'})
                        sp = mean(pks_s);
                        dp = mean(pks_d);
                        mp = mean(data);  % TODO why use this heuristic (1/3)*sp + (2/3)*dp when I can compute it exactly?
                        myMetrics = [myMetrics sp dp mp];
                    elseif ismember(chanNames{ichan}, {'dRVP/dt', 'dLVP/dt'})
                        maxdpdt = max(pks_s);
                        mindpdt = min(pks_d);
                        myMetrics = [myMetrics maxdpdt mindpdt];
                    else
                        error('unhandled channel name %s\n', chanNames{ichan});
                    end
                    
                    switch chanNames{ichan}
                        case 'PAP'
                            results(trainct).(sprintf('sPAP_right_%s', curphase)) = reshape(pks_s, 1, []);
                            results(trainct).(sprintf('dPAP_right_%s', curphase)) = reshape(pks_d, 1, []);
                        case 'RVP'
                            results(trainct).(sprintf('sRVP_right_%s', curphase)) = reshape(pks_s, 1, []);
                            results(trainct).(sprintf('dRVP_right_%s', curphase)) = reshape(pks_d, 1, []);
                        case 'dRVP/dt'
                            results(trainct).(sprintf('dRVPdtmax_right_%s', curphase)) = reshape(pks_s, 1, []);
                            results(trainct).(sprintf('dRVPdtmin_right_%s', curphase)) = reshape(pks_d, 1, []);
                        case 'RAP'
                            results(trainct).(sprintf('sRAP_right_%s', curphase)) = reshape(pks_s, 1, []);
                            results(trainct).(sprintf('dRAP_right_%s', curphase)) = reshape(pks_d, 1, []);
                        case 'SAP'
                            results(trainct).(sprintf('sSAP_left_%s', curphase)) = reshape(pks_s, 1, []);
                            results(trainct).(sprintf('dSAP_left_%s', curphase)) = reshape(pks_d, 1, []);
                        case 'LVP'
                            results(trainct).(sprintf('sLVP_left_%s', curphase)) = reshape(pks_s, 1, []);
                            results(trainct).(sprintf('dLVP_left_%s', curphase)) = reshape(pks_d, 1, []);
                        case 'dLVP/dt'
                            results(trainct).(sprintf('dLVPdtmax_left_%s', curphase)) = reshape(pks_s, 1, []);
                            results(trainct).(sprintf('dLVPdtmin_left_%s', curphase)) = reshape(pks_d, 1, []);
                        case 'LAP'
                            results(trainct).(sprintf('sLAP_left_%s', curphase)) = reshape(pks_s, 1, []);
                            results(trainct).(sprintf('dLAP_left_%s', curphase)) = reshape(pks_d, 1, []);
                    end
                    
                end
                
                if isequal(curtrain.repset, 'post')
                    trainct = trainct + 1;
                end
                    
                drawnow;
                % add a row to the table
                temp = arrayfun(@(x) myMetrics(x), 1:length(myMetrics), 'UniformOutput', false);                    
                mytable = [mytable; [{rowstr} temp]];
            end        
        end
        
        function [preTrains, postTrains] = create_pre_post_trains(obj, Twindow)
            % reference the times to the first segtime to convert to seconds
            t0_ = obj.segtimes(1);
            % get all trains
            trains = [obj.trains{:}];
            % get the start time of each train
            trainTimes = [trains.t0];
            assert(issorted(trainTimes), 'create_pre_post_trains: trains are not chronological');
            % get the duration of each train
            trainDur = arrayfun(@(x) (size(x.data, 1) - x.n1 - x.n2) / x.sr, trains);
            
            % initialize the output
            preTrains = [];
            postTrains = [];
            % loop over each physiological data record (files and records)
            for iseg = 1:length(obj.physdata)
                for irecord = 1:obj.physdata(iseg).file_meta.n_records
                    % get the start time and duration of the current record
                    t0rec = (obj.physdata(iseg).record_meta(irecord).record_start - t0_) * 86400;
                    Trec = obj.physdata(iseg).record_meta(irecord).tick_dt * obj.physdata(iseg).record_meta(irecord).n_ticks;
                    % determine which trains are contained in the current record
                    inds = find(trainTimes >= t0rec & trainTimes < (t0rec + Trec));
                    % handle the first and last trains separately
                    for itrain = 1:length(inds)
                        data = (0:round(Twindow*obj.sr)-1)';
                        t0 = trainTimes(inds(itrain)) - Twindow;
                        t0B = trainTimes(inds(itrain)) + trainDur(inds(itrain));
                        trigchan_ = zeros(size(data));

                        pt = PulseTrain(data, obj.sr, t0, 0, 0, trigchan_, NaN, NaN, 'pre', false, [], []);
                        pt.comments = {'pre'};
                        preTrains = [preTrains pt];


                        pt2 = PulseTrain(data, obj.sr, t0B, 0, 0, trigchan_, NaN, NaN, 'post', false, [], []);
                        pt2.comments = {'post'};
                        postTrains = [postTrains pt2];
                    end
                end
            end
            
            obj.store_phys_in_trains([preTrains postTrains]);
        end
        
        function list_comments(obj)
            commentList = [];
            for iseg = 1:length(obj.physdata)
                commentList = [commentList arrayfun(@(icomment) sprintf('Segment %d, Record %d, Index %d, Comment %s', ...
                    iseg, obj.physdata(iseg).comments(icomment).record, ...
                    icomment, obj.physdata(iseg).comments(icomment).str), ...
                    1:length(obj.physdata(iseg).comments), 'UniformOutput', false)];                
            end
            disp(char(commentList));
        end
        
        function foundTrains = find_train_by_comment(obj, myregexp)
            % look for comments that contain all of the following substrings            
            foundTrains = [];
            % loop over all neural segments
            for iseg = 1:length(obj.trains)
                % loop over trains
                for itrain = 1:length(obj.trains{iseg})
                    % loop over comments associated with the train
                    for icomment = 1:length(obj.trains{iseg}(itrain).comments)
                        if ~isempty(regexpi(obj.trains{iseg}(itrain).comments{icomment}, myregexp, 'match'))
                        % if all(cellfun(@(x) contains_(obj.trains{iseg}(itrain).comments{icomment}, x), mystr))
                            foundTrains = [foundTrains [iseg;itrain]];
                        end                        
                    end
                end
            end
            % should find a single consecutive group of trains
            % assert(all(diff(foundTrains(2, :))==1));
        end
        
        function plot_neural(obj, ts, units)            
            neuralSamplesStruct = obj.get_samples_from_time(ts, units, 'intan');            
            
            % plot the neural data (example uses neuralSamplesStruct)
            neuraldata = [];
            t = [];
            for istruct = 1:length(neuralSamplesStruct)  % will be a scalar struct unless the data spans multiple segments
                iseg = neuralSamplesStruct(istruct).seg;
                inds = neuralSamplesStruct.chaninds(1):neuralSamplesStruct.chaninds(2);
                neuraldata = [neuraldata; obj.neuraldata{iseg}(inds, :)];
                t0 = (obj.segtimes(iseg) - obj.segtimes(1)) * 86400;
                t = [t; t0 + (inds - 1) / obj.sr];
            end
            figure;plot(t, neuraldata);
        end
        
        function plot_physio(obj, ts, units)            
            physioSamplesStruct = obj.get_samples_from_time(ts, units, 'adi');            
            
            % plot the physiological data
            channames = {obj.physdata(1).channel_meta.name};
            figure;
            for ichan = 1:obj.physdata(1).file_meta.n_channels  % assume all adi segments contain the same channels
                ax = subplot(length(channames), 1, ichan);
                data = [];
                t = [];
                for istruct = 1:length(physioSamplesStruct)  % will be a scalar struct unless the data spans multiple records/segments
                    iseg = physioSamplesStruct(istruct).seg;
                    irecord = physioSamplesStruct(istruct).record;
                    dt = obj.physdata(iseg).channel_meta(ichan).dt(irecord);
                    inds = physioSamplesStruct.chaninds(ichan, 1):physioSamplesStruct.chaninds(ichan, 2);
                    channame = sprintf('data__chan_%d_rec_%d', ichan, irecord);
                    data = [data; obj.physdata(iseg).(channame)(inds)];
                    t0 = (obj.physdata(iseg).record_meta(irecord).record_start - obj.segtimes(1)) * 86400;
                    t = [t; t0 + (inds - 1) * dt];
                end
                plot(ax, t, data);
                title(ax, channames{ichan});
            end
        end
        
    end % of public methods
    
    methods (Access = private)
        function [y, pols] = stack_query(obj, traininds, pulseinds)
            %
            y = [];
            pols = [];
            for iseg = 1:length(traininds)
                % default to all pulses
                if ~exist('pulseinds', 'var') || isempty(pulseinds)
                    pulseinds = 1:max(arrayfun(@(x) length(x.pulses), obj.trains{iseg}(traininds{iseg})));
                end
                
                % cell arrays {[numsamples x numpulses] x numtrains}
                baseline = arrayfun(@(x) x.baseline(:, ismember(1:length(x.pulses), pulseinds)), obj.trains{iseg}(traininds{iseg}), 'UniformOutput', false);
                evokedpotentials = arrayfun(@(x) x.evokedpotentials(:, ismember(1:length(x.pulses), pulseinds)), obj.trains{iseg}(traininds{iseg}), 'UniformOutput', false);
                pols_ = arrayfun(@(x) repmat(x.pol, [1 nnz(ismember(1:length(x.pulses), pulseinds))]), obj.trains{iseg}(traininds{iseg}), 'UniformOutput', false);
                
                if verLessThan('matlab', '9.1')
                    y = [y bsxfun(@minus, [cell2mat(baseline); cell2mat(evokedpotentials)], median(cell2mat(baseline), 1))];
                else
                    % Matlab 2016b and later implements implicit expansion
                    y = [y [cell2mat(baseline); cell2mat(evokedpotentials)] - median(cell2mat(baseline), 1)];
                end
                pols = [pols cell2mat(pols_)];
            end
        end
        
        function extract_neural_trains(obj, T1, T2, minT, dtseg, tcomments, setidinds, plotflag)
            [repsets] = deal(cell(1, length(obj.trigchan)));
            for iseg = 1:length(obj.trigchan)
                % They might have forgot to store the trigger channel
                if isempty(obj.trigchan{iseg})
                    continue
                end

                inds1 = find(obj.trigchan{iseg}(1:end-1, 1) == 0 & obj.trigchan{iseg}(2:end, 1) == 1);  % rising edges
                inds2 = find(obj.trigchan{iseg}(1:end-1, 1) == 1 & obj.trigchan{iseg}(2:end, 1) == 0);  % falling edges
                assert(length(inds1) == length(inds2));
                if isempty(inds1) && isempty(inds2)
                    continue
                end

                n1 = round(T1*obj.sr);
                n2 = round(T2*obj.sr);

                % I handled these conditions with padding
                % assert(inds1(1) - n1 > 0);
                % assert(inds2(end) + n2 <= size(obj.neuraldata{iseg}, 1));

                obj.trains{iseg} = [];
                trainct = 0; % TODO should I assume that the neural data won't start and stop during a train?
                curind = [0 0];  % current setID comment
                for itrain = 1:length(inds1)
                    if (inds2(itrain) - inds1(itrain)) / obj.sr > minT
                        % trainct = trainct + 1;

                        if inds1(itrain) <= n1
                            n0 = 1;                            
                        else
                            n0 = (inds1(itrain) - n1);                            
                        end
                        pad1 = max(0, n1 - inds1(itrain) + 1);
                        if inds2(itrain) + n2 > size(obj.neuraldata{iseg}, 1)
                            nf = size(obj.neuraldata{iseg}, 1);
                        else
                            nf = (inds2(itrain) + n2);
                        end
                        pad2 = max(0, inds2(itrain) + n2 - size(obj.neuraldata{iseg}, 1));

                        % data = chanfcn(obj.neuraldata{iseg}(n0:nf, :));
                        numchans = size(obj.neuraldata{iseg}, 2);
                        data = [NaN(pad1, numchans); obj.neuraldata{iseg}(n0:nf, :); NaN(pad2, numchans)];
                        t0 = dtseg(iseg) + n0 / obj.sr;  % Reference the time to the start of the first segment

                        % TODO map comment to repset and keep track of ampgain and pol
                        for iadi = length(tcomments):-1:1
                            halfwidth = (inds2(itrain) - inds1(itrain)) / obj.sr / 2;
                            ind = find(t0 > tcomments{iadi} - halfwidth, 1, 'last');
                            if ~isempty(ind)
                                break
                            end
                        end
                        % ind and iadi

                        % I left this alone, but the values are mostly overwritten by associate_repset_comments_with_trains
                        isemptytc = all(cellfun(@(tc) isempty(tc), tcomments));
                        if isempty(ind)
                            if isemptytc
                                warning('current train preceeds all setID comments\n');
                            end

                            repset = '0.0';
                            ampgain = 1;
                            pol = 0;  % unknown polarity

                            params = struct('pulsewidth_us',NaN,'freq_Hz',NaN,'numpulses',NaN,'duration_s',NaN,'iti_s',NaN,'wait_min',NaN,'distbwcuffs',NaN);

                        elseif isequal([ind iadi], curind)
                            if ~isempty(obj.wparams) && obj.wparams.altpol
                                pol = -pol;
                            end
                            trainct = trainct + 1;
%                             params = obj.wparams.params(repset);
                        else
                            if ~isempty(obj.wparams) && mod(trainct, obj.wparams.trainsperset) ~= 0
                                warning('Unexpected number of trains');
                            end

                            curind = [ind iadi];
                            assert(~isempty(obj.physdata) && ~isemptytc && ~isempty(ind));
                            % support array of PhysioData objects
                            [repset, ampgain] = obj.physdata(iadi).parse_comment_str(setidinds{iadi}(ind));
                            pol = 1;  % TODO Is positive polarity always first?
                            trainct = 1;
                            repsets{iseg} = [repsets{iseg} {repset}];

                            if str2double(repset) < 100 && ~isequal(repset, '0.0')
                                params = obj.wparams.params(repset);
                                params.distbwcuffs = obj.wparams.distbwcuffs;
                            end
                        end

                        % Get the PhysioData record for the PulseTrain constructor
                        % TODO Why do we need the record by not the segment?
                        % if ~isempty(obj.physdata)
                        %     ind_ = obj.get_samples_from_time([t0 t0], 'sec', 'adi');
                        %     record = ind_.record;
                        % else
                        %     record = [];
                        % end
                        
                        if size(obj.trigchan{iseg}, 2) == 1
                            obj.trains{iseg} = [obj.trains{iseg} PulseTrain(data, obj.sr, t0, n1, n2, [], pol, ampgain, repset, plotflag, params)];
                        else                            
                            obj.trains{iseg} = [obj.trains{iseg} PulseTrain(data, obj.sr, t0, n1, n2, obj.trigchan{iseg}(n0:nf, 2), pol, ampgain, repset, plotflag, params)];
                        end

                        % I don't think obj.repsetlist needs to be divided into segments
                        if ~isempty(repset) && (isempty(obj.repsetlist) || ~ismember({repset}, obj.repsetlist))
                            obj.repsetlist = [obj.repsetlist {repset}];
                        end

                    end
                end
            end
        end
    end % of private methods
    
    methods (Static)
        function add_legend(ax, repsets, pols, cols)
            % plot NaNs on the current axes to get a subset for the legend
            legstr = cell(1, length(repsets) * length(pols));
            polmap = containers.Map([-1 1], {'-' '+'});
            ct = 1;
            for irs = 1:length(repsets)
                for ipol = 1:length(pols)
                    legstr{ct} = [repsets{irs} ' ' polmap(pols(ipol))];
                    ct = ct + 1;
                end
            end
            
            hold(ax, 'on');
            if verLessThan('matlab', '9.3')
                s_ = arrayfun(@(icol) plot(ax, NaN, NaN, 'Color', cols(icol, :)), 1:size(cols, 1), 'UniformOutput', false);
                s = [];
                for ii = 1:length(s_)
                    % cell2mat doesn't work either
                    s = [s s_{ii}];
                end
            else
                % not sure what version of Matlab supports arrayfun Line array (2017b 9.3?)
                s = arrayfun(@(icol) plot(ax, NaN, NaN, 'Color', cols(icol, :)), 1:size(cols, 1));
            end
            legend(s, legstr);
            hold(ax, 'off');
        end
        
        % Compute Fisher's linear discriminant
        function S = ldastat(data1, data2)
            % compute Fisher's linear discriminant statistic
            % compute covariance matrices for the data
            sigma1 = cov(data1'); sigma2 = cov(data2');
            sigma_avg = (sigma1 + sigma2)./2;
            
            % decision criteria
            mu1 = mean(data1, 2);
            mu2 = mean(data2, 2);
            w_ = inv(sigma_avg)*(mu1 - mu2);
            
            % Fisher's statistic
            numerator = (dot(w_, (mu1 - mu2)))^2; % between-cluster variance
            denominator = w_'*(sigma1 + sigma2)*w_; % within-cluster variance
            S = numerator/denominator;
        end
        
        function [allTrains] = sort_trains(baselineTrains, trains)
           allTrains = [baselineTrains trains];
           [~, si] = sort([allTrains.t0]);
           allTrains = allTrains(si);
        end
        
        function write_results_table_for_pressure(mytable, filename)
            % format an Excel spreadsheet as Subash has done based on the data in the table mytable
            % I'm intentionally not going to concatenate rows because this would become very long
            assert(mod(size(mytable, 1), 3) == 0, 'mytable must be [pre during post] * N');
            
            temp = cellfun(@(varname) cellfun(@(phase) sprintf('%s_%s', phase, varname), {'Pre' 'Dur'}, 'UniformOutput', false), mytable.Properties.VariableNames(13:end), 'UniformOutput', false);
            temp = [temp{:}];            
            header = [mytable.Properties.VariableNames(9:12) temp];
            
            mysheet = header;
            for itrain = 1:3:size(mytable, 1)
                data = table2cell(mytable(itrain:itrain+2, :));
                [repid, setid, pol] = data{2, 10:12};
                myrow = [data(2, 9:12) reshape(data(1:2, 13:end), 1, [])];
                mysheet = [mysheet; myrow];
            end
            
            xlswrite(filename, mysheet);
        end
        
    end
end % of class UTExperiment
