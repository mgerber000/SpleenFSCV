classdef PhysioData < dynamicprops
    properties
        % comments    % comments structure
        inputfile
        % record_meta  % record meta data
        % channel_meta
        % stimchans  % cell array of stim channel records
        % toffset  % vector of length record_meta
        ecg           % raw ecg signal
        ecgchan
        respchan
        hr            % heart rate
        hrchan
        brchan
        thr           % time vector for heart rate
        ecgid         % ECG channel ID
        resp          % respiratory signal
        br            % breathing rate
        tbr           % time vector for breathing rate
        respid      % channel ID for respiratory signal
        hpks        % ecg peak values
        bpks        % breathing peak values
        ecg_params % store parameters used for ECG signal analysis
        breath_params % store parameters used for breathing signal
        trigid        %MCS train channel ID
        trigchan      %MCS train raw signal 
        triggers      %Rising and falling edges of signal
        repsetregexp
        repsetxregexp
    end
    
    methods
        function obj = PhysioData(aidchtfile, repsetregexp, repsetxregexp)
            obj.inputfile = aidchtfile;
            obj.repsetregexp = repsetregexp;
            obj.repsetxregexp = repsetxregexp;
            
            [filepath, name, ext] = fileparts(aidchtfile);
            switch ext
                case '.mat'
                    obj.load_mat_file(aidchtfile);
                case '.adicht'
                    matfile_ = [filepath filesep name '.mat'];
                    fprintf('converting %s to %s\n', aidchtfile, matfile_);
                    adi.convert(aidchtfile, 'save path', matfile_);
                    obj.load_mat_file(matfile_);
                    % case '.xlsx'
                    %     obj.toffset = 2151.6;  % hardcoded for intan427
                    %     [~, txt] = xlsread(aidchtfile);
                    %     obj.comments = txt(2:end, 4);
                    %     timestrings = txt(2:end, 3);
                    %     obj.commenttimes = cellfun(@(x) timestr2sec(x), timestrings) - obj.toffset;
                otherwise
                    error('%s not a valid aidcht file', aidchtfile);
            end
            
            assert(isequal([obj.record_meta.record_start], [obj.record_meta.data_start]));
        end
        
        function t = comment_time(obj, inds)
            t = zeros(1, length(inds));
            for ct = 1:length(inds)
                ind = inds(ct);
                % get the datenum for the comment
                record = obj.comments(ind).record;
                assert(isequal(obj.comments(ind).tick_dt, obj.record_meta(record).tick_dt));
                dt = obj.record_meta(record).tick_dt;
                t(ct) = obj.record_meta(record).record_start + obj.comments(ind).tick_position * dt / 86400;                
            end
        end
        
        function [repset, ampgain] = parse_comment_str(obj, ind)
            tokens = regexpi(obj.comments(ind).str, obj.repsetregexp, 'tokens');
            tokens2 = regexpi(obj.comments(ind).str, obj.repsetxregexp, 'tokens');
            
            if ~isempty(tokens) && ~isempty(tokens2)
                repset = tokens{1}{1};
                ampgain = [];  % TODO haven't extracted this token yet
            elseif ~isempty(tokens) && isempty(tokens2)
                repset = '0.0';
                ampgain = [];
            else
                % 'Begin set IDs' for example
                repset = '';
                ampgain = [];
            end
        end
        

        function [elec] = parse_elec_comment_str(obj, ind, mysymbol)
            if contains_(obj.comments(ind).str, 'record from', 'IgnoreCase', true)
                tokens = regexpi(obj.comments(ind).str, 'R?(\d+(?:,\s*\d+)*)(?:;S)?(\d+(?:,\s*\d+)*)?$', 'tokens');           
            else
                tokens = regexpi(obj.comments(ind).str, ['(?:' mysymbol '=?)(\d+(?:,\s*\d+)*)'], 'tokens');            
            end
            % tokens = regexpi(obj.comments(ind).str, '(?:R=?)(\d+(?:,\s*\d+)*)(?:;(?:S=?)?)(\d+(?:,\s*\d+)*)?$', 'tokens');
            % tokens = regexpi(obj.comments(ind).str, 'R?(\d+(?:,\s*\d+)*)(?:;S)?(\d+(?:,\s*\d+)*)?$', 'tokens');
            if ~isempty(tokens)
                elec = str2num(tokens{1}{1});
            else
                % 'Begin set IDs' for example
                elec = [];
            end
        end
        
        function obj = getheartrate(obj, method, minpeakprom, mindist, plotvar,neg)
            % Identify ECG channel, detect QRS complexes and compute heart
            % rate
            
            % Identify ECG channel
            chan_names = {obj.channel_meta(:).name};
            ecgidx = cell2mat(cellfun(@(x) contains_(x, 'ECG'), chan_names, 'un',0));
            ecgchan = find(ecgidx);
            obj.ecgid = ecgchan;
            numrecords = obj.file_meta.n_records;
            [obj.hr, obj.thr] = deal([]);
            obj.ecg_params = repmat(struct('mindist', 0, 'minpeakprom', 0), [1 numrecords]);
%            hmsg = MessageUpdater();
            for rec = 1:numrecords
                % Get sampling rate of ECG channel
                chan_dt = [obj.channel_meta(ecgchan).dt(rec)];
                if isnan(chan_dt)
                    continue
                end
                ecg_sr = round(1/chan_dt);
                
                % Assign ECG channel
                if isequal(obj.inputfile, 'E:\Data\UT - PH new\2018-06-26-u-1\exp-2018-06-26-u-1-1.mat') && rec == 2
                    obj.ecg{rec, 1} = -double(obj.(sprintf('data__chan_%d_rec_%d', ecgchan, rec)));
                else
                    obj.ecg{rec, 1} = double(obj.(sprintf('data__chan_%d_rec_%d', ecgchan, rec)));
                end
                % eval(sprintf('obj.ecg{%d,1} = double(obj.data__chan_%d_rec_%d);', rec, ecgchan, rec));
                %             eval(sprintf('obj.data__chan_%d_rec_1 = [];', ecgchan));
                
                %High pass filter signal to remove movement artifacts
                [b,a] = butter(2, 0.1/(ecg_sr/2), 'high');
%                 [A,B,C,D] = butter(2,[0.1 10]/(ecg_sr/2),'bandpass');
%                 [sos,g] = ss2sos(A,B,C,D);
                filt_ecg = filtfilt(b, a, obj.ecg{rec});
                if neg == true
                    filt_ecg = -filt_ecg;
                end
                
                % Detect peaks
                if ~exist('mindist', 'var') || isempty(mindist)
                    maxhr = 800; % 550 beats per min
                    mindist = 1/(maxhr/60);
                end
                t = linspace(0, numel(obj.ecg{rec})/ecg_sr, numel(obj.ecg{rec}));
                
                [pks, locs] = findpeaks(filt_ecg, t, 'MinPeakDistance', mindist, 'MinPeakProminence', minpeakprom);
                if isempty(locs)
                    continue
                end
                
                switch method
                    case 'template'
                        % Find template of ECG - QRS complexes
%                         temp_locs = find(ismember(t,locs(locs < 10))); % Get location of ECG peaks less than 10 seconds
%                         numsamples = 100;
%                         temp_locs = temp_locs(temp_locs > numsamples);
%                         gen_inds = arrayfun(@(x) x-numsamples:x+numsamples, temp_locs, 'un',0);
%                         template_inds = vertcat(gen_inds{:});
%                         conv_ecg = mean(filt_ecg(template_inds));
%                         cpdt = conv(filt_ecg, conv_ecg);
                end
                ibi = diff(locs); % inter-beat-interval
                obj.hr{rec} = 60./[0 ibi]; % heart rate
                diffhr = diff([0 obj.hr{rec}]);
                
                % Correct spuriously detected peaks
                diffthresh = 150;
                counter = 0;
                
                while sum(abs(diffhr > diffthresh)) > 0 && counter < 10000
                    counter = counter + 1;
                    msg = sprintf('Correcting heart rate for %d points..', sum(abs(diffhr > 35)));
 %                   hmsg.update_message(msg);                    
                    excinds = find(diffhr > diffthresh | diffhr < -diffthresh);
                    numpoints = min([10, excinds(1)-1]);
                    %numpoints = cell2mat(arrayfun(@(x) min([10, x]), excinds, 'un',0));
%                     obj.hr{rec}(excinds(1)) = nanmean(obj.hr{rec}(excinds(1)- numpoints: excinds(1) - 1));
                     diffhr = diff([0 obj.hr{rec}]);
                    excinds = find(diffhr > diffthresh | diffhr < -diffthresh);
                end
                
                diffhr = diff([0 obj.hr{rec}]);
                
                % print warning message about uncorrected spurious peaks
                if sum(abs(diffhr > diffthresh)) > 0
                    warning('%d spurious detections in ECG uncorrected...', sum(abs(diffhr > diffthresh)));
                end
                
                % Peak detection figure
                if plotvar
                    figure;
                    ax1 = subplot(311);
                    reduce_plot(t, filt_ecg);
                    ylim([-1,1])
                    hold on;
                    reduce_plot(locs, pks, 'r*');
                    ax2 = subplot(312);
                    plot(locs, obj.hr{rec});
                    ax3 = subplot(313);
                    plot(locs, diffhr)
                    linkaxes([ax1, ax2, ax3], 'x');
                end
                
                obj.thr{rec} = locs;
                obj.hpks{rec} = obj.ecg{rec}(ismember(t, locs));
                obj.ecg_params(rec) = struct('mindist',mindist, 'minpeakprom',minpeakprom);
            end
            %clear hmsg
        end
        
        function obj = getresprate(obj, minpeakprom, smoothvar, plotvar, neg)
            % Identify resp channel, detect inhalation peaks and compute
            % respiratory rate
            
            % Identify resp channel
            if ~exist('minpeakprom','var') || isempty(minpeakprom), minpeakprom = 0.05; end
            chan_names = {obj.channel_meta(:).name};
            respidx = cell2mat(cellfun(@(x) contains_(x, 'Resp', 'IgnoreCase',true) || contains_(x, 'Nasal', 'IgnoreCase',true), chan_names, 'un',0));
            respchan = find(respidx);
            obj.respid = respchan;
            numrecords = obj.file_meta.n_records;
            obj.breath_params = repmat(struct('mindist', 0, 'minpeakprom', 0, 'smoothvar', 0), [1 numrecords]);
            
            for rec = 1:numrecords
                % Get sampling rate of ECG channel
                chan_dt = [obj.channel_meta(respchan).dt];
                resp_sr = round(1/chan_dt(rec));
                
                % Assign ECG channel
                eval(sprintf('obj.resp{%d,1} = double(obj.data__chan_%d_rec_%d);', rec, respchan, rec));
                
                % High pass filter signal to remove movement artifacts
                [b,a] = butter(2, 0.1/(resp_sr/2), 'high');
                filt_resp = filtfilt(b, a, obj.resp{rec,1});
                filt_resp_smooth = smooth(filt_resp, smoothvar); 
                if neg == true
                    filt_resp_smooth = -filt_resp_smooth;
                end
                
                % Detect peaks
                maxbr = 120; % 120 breaths per min
                mindist = 1/(maxbr/60);
                t = linspace(0, numel(obj.resp{rec,1})/resp_sr, numel(obj.resp{rec,1}));
                [pks, locs] = findpeaks(filt_resp_smooth, t, 'MinPeakDistance', mindist, 'MinPeakProminence', minpeakprom);
                if isempty(locs)
                    disp('No respiratory peaks were detected with the given peak prominence and smoothing values')
                    continue
                end
                ibi = diff(locs); % inter-beat-interval
                obj.br{rec,1} = 60./[0 ibi]; % heart rate
                diffbr = diff([0 obj.br{rec,1}]);
                % Peak detection figure
                if plotvar
                    figure;
                    ax1 = subplot(311);
                    reduce_plot(t, filt_resp);
                    hold on;
                    reduce_plot(t, filt_resp_smooth,'c');
                    reduce_plot(locs, pks, 'r*');
                    ax2 = subplot(312);
                    plot(locs, obj.br{rec,1});
                    ax3 = subplot(313);
                    plot(locs, diffbr)
                    linkaxes([ax1, ax2, ax3], 'x');
                end
                
                excinds = obj.br{rec,1} < 0 | obj.br{rec,1} > 120;
                obj.br{rec,1}(excinds) = nanmean(obj.br{rec,1}(~excinds));
                obj.tbr{rec,1} = locs;
                obj.bpks{rec,1} = obj.resp{rec,1}(ismember(t, locs));
                obj.breath_params(rec) = struct('minpeakprom',minpeakprom, 'mindist',mindist, 'smoothvar',smoothvar);
            end
        end  
        
        function [trains, repsetlist, pol] = extract_trains(obj, T1, T2, minT, tcomments, setidinds, wparams, pol, segtime0, plotflag)         
            ichan = find(cellfun(@(x) contains_(x, 'train', 'IgnoreCase', true), {obj.channel_meta.name}));            
            assert(length(ichan) == 1, 'did not uniquely identify the trigger channel');
            
            [trains, repsets] = deal(cell(1, obj.file_meta.n_records));
            repsetlist = [];
            for irecord = 1:obj.file_meta.n_records
                trigchan = obj.(sprintf('data__chan_%d_rec_%d', ichan, irecord));
                
                inds1 = find(trigchan(1:end-1, 1) < 2 & trigchan(2:end, 1) >= 2);  % rising edges
                inds2 = find(trigchan(1:end-1, 1) >= 2 & trigchan(2:end, 1) < 2);  % falling edges
                
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
                    continue
                end
                
                sr = round(1/obj.channel_meta(ichan).dt(irecord));
                nsamples = obj.channel_meta(ichan).n_samples(irecord);
                
                n1 = round(T1*sr);
                n2 = round(T2*sr);
                
                % I handled these conditions with padding
                % assert(inds1(1) - n1 > 0);
                % assert(inds2(end) + n2 <= size(obj.neuraldata{iseg}, 1));
                
                % obj.trains{irecord} = [];
                trainct = 0; % TODO should I assume that the neural data won't start and stop during a train?
                curind = 0;  % current setID comment
                for itrain = 1:length(inds1)
                    if (inds2(itrain) - inds1(itrain)) / sr > minT
                        % trainct = trainct + 1;
                        
                        if inds1(itrain) <= n1
                            n0 = 1;
                        else
                            n0 = (inds1(itrain) - n1);
                        end
                        % pad1 = max(0, n1 - inds1(itrain) + 1);
                        
                        
                        if inds2(itrain) + n2 > nsamples
                            nf = nsamples;
                        else
                            nf = (inds2(itrain) + n2);
                        end
                        % pad2 = max(0, inds2(itrain) + n2 - nsamples);
%                         
%                         % data = chanfcn(obj.neuraldata{iseg}(n0:nf, :));
%                         numchans = size(obj.neuraldata{iseg}, 2);
%                         data = [NaN(pad1, numchans); obj.neuraldata{iseg}(n0:nf, :); NaN(pad2, numchans)];
%                         t0 = dtseg(iseg) + n0 / obj.sr;  % Reference the time to the start of the first segment

                        segtime0 = obj.record_meta(irecord).record_start;
                        t0 = (obj.record_meta(irecord).record_start - segtime0) * 86400 + inds1(itrain) / sr ;
                        
                        % TODO the loop is outside of this function so I may need to pass in the polarity state
                        % TODO map comment to repset and keep track of ampgain and pol
                        % for iadi = length(tcomments):-1:1
                            halfwidth = (inds2(itrain) - inds1(itrain)) / sr / 2;
                            % ind = find(t0 > tcomments{iadi} - halfwidth, 1, 'last');
                            ind = find(t0 > tcomments - halfwidth, 1, 'last');
                            % if ~isempty(ind)
                            %     break
                            % end
                        % end
                        % ind and iadi
                        
                        % I left this alone, but the values are mostly overwritten by associate_repset_comments_with_trains
                        % isemptytc = all(cellfun(@(tc) isempty(tc), tcomments));
                        params = struct('pulsewidth_us',NaN,'freq_Hz',NaN,'numpulses',NaN,'duration_s',NaN,'iti_s',NaN,'wait_min',NaN,'distbwcuffs',NaN);
                        isemptytc = isempty(tcomments);
                        if isempty(ind)
                            if isemptytc
                                warning('current train preceeds all setID comments\n');
                            end
                            repset = '0.0';
                            ampgain = 1;
                            pol = 0;  % unknown polarity                            
                        elseif isequal(ind, curind)
                            if ~isempty(wparams) && wparams.altpol
                                pol = -pol;
                            end
                            trainct = trainct + 1;                            
                            %                             params = obj.wparams.params(repset);
                        else
                            if ~isempty(wparams) && mod(trainct, wparams.trainsperset) ~= 0
                                warning('Unexpected number of trains');
                            end
                            
                            curind = ind;
                            assert(~isemptytc && ~isempty(ind));
                            % support array of PhysioData objects
                            [repset, ampgain] = obj.parse_comment_str(setidinds(ind));
                            pol = 1;  % TODO Is positive polarity always first?
                            trainct = 1;
                            repsets{irecord} = [repsets{irecord} {repset}];
                            if str2double(repset) < 100 && ismember(repset, wparams.params.keys())                                
                                params = wparams.params(repset);
                                params.distbwcuffs = wparams.distbwcuffs;
                            end
                        end  
                        trains{irecord} = [trains{irecord} PulseTrain((n0:nf)', sr, t0, n1, n2, trigchan(n0:nf, 1), pol, ampgain, repset, plotflag, params, irecord)];

                        % I don't think obj.repsetlist needs to be divided into segments
                        if ~isempty(repset) && (isempty(repsetlist) || ~ismember({repset}, repsetlist))
                            repsetlist = [repsetlist {repset}];
                        end
                        
                    end
                end                
            end
            
            % concatenate trains from the records
            trains_ = [];
            for irecord = 1:length(trains)
                trains_ = [trains_ trains{irecord}];
            end
            trains = trains_;
        end
        
        function [dx, dt, ind] = derivative_chan(obj, name, fc)
            % construct a channel containing the derivative of the pressure            
            
            ind = find(ismember({obj.channel_meta.name}, name));
            if length(ind) ~= 1
                dx = [];
                dt = [];
                fprintf('%s not found\n', name);
                return
            end
            
            % smoothing parameter
            if ~exist('fc', 'var') || isempty(fc)
                fc = 50;
            end            
            
            % get the pressure channel                        
            dx = cell(1, obj.file_meta.n_records);
            dt = zeros(1, obj.file_meta.n_records);
            for irecord = 1:obj.file_meta.n_records
                p = double(obj.(sprintf('data__chan_%d_rec_%d', ind, irecord)));
                dt_ = obj.channel_meta(ind).dt(irecord);
                
                % Pressure should be at 1 kHz, but I observed it at 40 kHz
                if dt_ < 1e-4
                    fprintf('dt = %f. Resample the channel first.\n', dt_);
                    return
                    % p = resample(p, 1, round(.001/dt_));
                    % dt_ = dt_ * round(.001/dt_);
                end
                dt(irecord) = dt_;
                
                % piecewise polynomial smoothing spline interpolation
                t = reshape((0:length(p)-1)*dt_, [], 1);                
                % pp = csaps(t, p, sp);
                
                % Labchart uses something more similar to this
                sr = 1/dt(irecord);
                n = round(0.0384 * sr);
                b = fir1(n, fc/(sr/2), 'low', kaiser(n+1, 6));
                
                % take the derivative
                % dx{irecord} = gradient(ppval(pp, t), dt(irecord));
                dx{irecord} = gradient(filtfilt(b, 1, p), dt(irecord));
            end
        end
        
        function add_channel(obj, data, dt, name, units)
            nrecs = obj.file_meta.n_records;
            assert(iscell(data) && length(data) == nrecs);
            assert(length(dt) == nrecs);
            assert(iscell(units) && length(units) == nrecs);
            
            if ismember(name, {obj.channel_meta.name})
                fprintf('physiological data alread contains a channel named %s\n', name); 
                return
            end
            
            ichan = obj.file_meta.n_channels + 1;
            
            old_fm = obj.file_meta;
            old_cm = obj.channel_meta;
            try
                obj.file_meta.n_channels = ichan;

                chanmeta = struct('id', ichan, 'name', name, 'units', [], 'n_samples', cellfun(@(x) length(x), data), 'dt', dt);
                chanmeta.units = units;  % cannot assign a cell array using struct()

                obj.channel_meta = [obj.channel_meta chanmeta];
                for irecord = 1:nrecs
                    chanstr = sprintf('data__chan_%d_rec_%d', ichan, irecord);
                    obj.addprop(chanstr);
                    obj.(chanstr) = single(data{irecord});
                end
            catch ex
               disp(getReport(ex));
                % cleanup from an exception
               obj.file_meta = old_fm;
               obj.channel_meta = old_cm;
               if isprop(obj, chanstr)
                   rmprops(obj, chanstr);
               end
            end
        end
        
        function resample_chan(obj, name, sr)
            % dt changes across records
            ichan = find(ismember({obj.channel_meta.name}, name));
            % assert(length(ichan) == 1);
            if isempty(ichan)
                fprintf('did not find channel %s\n', name);
                return
            end
            
            cm = obj.channel_meta(ichan);
            oldData = cell(1, obj.file_meta.n_records);
            
            try
                for irecord = 1:obj.file_meta.n_records
                    [p, q] = rat(cm.dt(irecord) * sr, .01);
                    
                    datastr = sprintf('data__chan_%d_rec_%d', ichan, irecord);
                    oldData{irecord} = obj.(datastr);
                    
                    if (p == 1 && q == 1) || (p == 0 && q == 0) || isempty(obj.(datastr))
                        continue
                    end
                    
                    if isa(obj.(datastr), 'single')
                        obj.(datastr) = single(resample(double(obj.(datastr)), p, q));
                    else
                        obj.(datastr) = resample(obj.(datastr), p, q);
                    end
                    obj.channel_meta(ichan).n_samples(irecord) = length(obj.(datastr));
                    obj.channel_meta(ichan).dt(irecord) = q/p * cm.dt(irecord);
                end
            catch ex
                disp(getReport(ex));
                
                obj.channel_meta(ichan) = cm;
                for irecord = 1:obj.file_meta.n_records
                    datastr = sprintf('data__chan_%d_rec_%d', ichan, irecord);
                    obj.(datastr) = oldData{irecord};
                end
            end
        end        
    end
    
    methods (Access = private)
        function load_mat_file(obj, matfile_)
            m = matfile(matfile_, 'Writable', false);
            if ~isprop(m, 'comments')
                fprintf('comments not found in aidcht mat file\n');
            end
            % obj.comments = m.comments;
            % obj.record_meta = m.record_meta;
            % obj.channel_meta = m.channel_meta;
            % also could be called 'MCS train sync' and 'MCS pulse sync'
            % containstrigger = arrayfun(@(cm) contains_(cm.name, 'trigger'), obj.channel_meta);
            % triggerchan = find(containstrigger);
            % if isempty(triggerchan)
            %     obj.stimchans = 'none';
            % else
            %     obj.stimchans = arrayfun(@(recnum) m.(sprintf('data__chan_%d_rec_%d', triggerchan, recnum)), 1:length(obj.record_meta), 'UniformOutput', false);
            % end
            myprops = properties(m);
            % datachans = reshape(find(cellfun(@(x) ~isempty(x), regexp(myprops, 'data__chan_\d+_rec_\d+', 'match'))), 1, []);
            for ichan = 1:length(myprops)
                if strcmp(myprops{ichan}, 'Properties')
                    continue
                end
                if ~isprop(obj, myprops{ichan})
                    obj.addprop(myprops{ichan});
                end
                obj.(myprops{ichan}) = m.(myprops{ichan});
            end
            if ~isprop(obj, 'comments')
                obj.addprop('comments');
            end
            assert(isprop(obj, 'record_meta'));
            assert(isprop(obj, 'channel_meta'));
        end
    end
end