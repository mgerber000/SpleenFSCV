classdef WaveformParams < handle
    properties
        datestr
        expid         %comment
        threshparams  % struct
        trainsperset  % includes alternating polarity if true
        altpol        % boolean
        setreps
        vnsside
        stimcuff
        reccuff
        distbwcuffs   % mm
        params     % containers.Map that takes 'r.s' and returns a struct
        timeest
    end
    
    methods
        function obj = WaveformParams(filename)
            [folderDir, ~, ext] = fileparts(filename);
            switch ext
                case {'.xls', '.xlsx'}
                    [numeric, ~, raw] = xlsread(filename);
                    obj.datestr = raw{1, 2};
                    obj.expid = raw{2, 2};
                    obj.threshparams = struct('pulsewidth_us', numeric(1, 2), 'freq_Hz', numeric(1, 3), 'numpulses', numeric(1, 4), 'duration_s', numeric(1, 5));
                    obj.trainsperset = numeric(3, 2);
                    obj.altpol = contains_(raw{7, 1}, 'with alt-pol');
                    obj.setreps = numeric(5, 2);
                    obj.vnsside = raw{10, 2};
                    obj.stimcuff = raw{11, 2};
                    obj.reccuff = raw{12, 2};
                    inds = 12:size(numeric, 1)-2;
                    keyset = arrayfun(@(repid, setid) [num2str(repid) '.' num2str(setid)], numeric(inds, 1), numeric(inds, 2), 'UniformOutput', false);
                    obj.params = containers.Map(keyset, ...
                        num2cell(arrayfun(@(ii) struct('pulsewidth_us', numeric(ii, 3), ...
                        'freq_Hz', numeric(ii, 4), 'numpulses', numeric(ii, 5), ...
                        'duration_s', numeric(ii, 6), 'amp', numeric(ii,7), 'iti_s', numeric(ii, 8), ...
                        'wait_min', numeric(ii, 9)), inds)));
                    obj.timeest = numeric(end, 2);
                case '.csv'
                    % implement csv reader
                    raw = WaveformParams.read_waveformparams_csv(filename, ',');
                    obj.datestr = raw{1, 2};
                    obj.expid = raw{2, 2};
                    obj.threshparams = struct('pulsewidth_us', str2double(raw{5, 2}), ...
                        'freq_Hz', str2double(raw{1, 3}), ...
                        'numpulses', str2double(raw{5, 4}), ...
                        'duration_s', str2double(raw{5, 5}));
                    obj.trainsperset = str2double(raw{7, 2});
                    obj.altpol = contains_(raw{7, 1}, 'with alt-pol');
                    obj.setreps = str2double(raw{9, 2});
                    obj.vnsside = raw{10, 2};
                    obj.stimcuff = raw{11, 2};
                    obj.reccuff = raw{12, 2};
                    ind1 = find(cellfun(@(x) isequal(x, 'Rep ID'), raw(:, 1)));
                    ind = ind1 + find(cellfun(@(x) isempty(x), raw((ind1+1):end, 1)), 1, 'first') - 1;
                    inds = (ind1+1):ind;
                    temp = cellfun(@(x) str2double(x), raw(inds, :));
                    keyset = arrayfun(@(repid, setid) [num2str(repid) '.' num2str(setid)], temp(:, 1), temp(:, 2), 'UniformOutput', false);
%                     obj.params = containers.Map(keyset, ...
%                         num2cell(arrayfun(@(ii) struct('pulsewidth_us', temp(ii, 3), ...
%                         'freq_Hz', temp(ii, 4), 'numpulses', temp(ii, 5), ...
%                         'duration_s', temp(ii, 6), 'amp', temp(ii,7), 'iti_s', temp(ii, 8), ...
%                         'wait_min', temp(ii, 9)), 1:length(temp))));
                    obj.params = containers.Map(keyset, ...
                        num2cell(arrayfun(@(ii) struct('pulsewidth_us', temp(ii, 3), ...
                        'freq_Hz', temp(ii, 4), 'numpulses', temp(ii, 5), ...
                        'duration_s', temp(ii, 6), 'amp', temp(ii,7), 'iti_s', temp(ii, 8), ...
                        'wait_min', temp(ii, 9)), 1:size(temp,1))));
                    obj.timeest = str2double(raw{inds(end)+2, 2});
            end
            
            try
                [~,foldername,~] = fileparts(folderDir);
                obj.distbwcuffs = getDistbwCuffs(foldername);
            catch
                obj.distbwcuffs = [];
            end
        end
    end
    
    methods (Static, Access = private)
        function raw = read_waveformparams_csv(filename, delim)
            fid = fopen(filename, 'r');
            hcleanup = onCleanup(@() fclose(fid));
            raw = [];
            while ~feof(fid)
                myline = fgetl(fid);
                splitline = strsplit(myline, [delim '(?=(?:[^"]*"[^"]*")*[^"]*$)'], 'DelimiterType', 'RegularExpression', 'CollapseDelimiters', false);
%                 splitline = strsplit(myline, delim, 'CollapseDelimiters', false);
                if ~isempty(raw)
                    assert(length(splitline) == size(raw, 2), ...
                        '%s and %s in csv file are inconsistent', ...
                        strjoin(raw(end, :), ', '), strjoin(splitline, ', '));
                end
                raw = [raw; splitline];
            end
        end
    end
end