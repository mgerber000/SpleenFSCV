function fdata = filter_data1(data, sample_rate, filter_func, filter_order, low_cut_Hz, high_cut_Hz, varargin)
% filter_data(data, filter_fnc, filter_band, low_cut_Hz, high_cut_Hz, order)
%   data -- array of data to filter.
%   sample_rate -- data rate in samples per second.
%   filter_func -- allowed values are 'butter', 'cheby2', 'cheby1'
%   varargin{1} -- Optional Matlab filter pass type: 'bandpass', 'stop', 'low', 'high'.
%       Default = 'bandpass'
%
% Filter data with an IIR filter

    filter_index = find(strcmp({'butter', 'cheby2', 'cheby1'}, filter_func));
    if isempty(filter_index)
        fdata = data; % Unrecognized filter type
        return;
    end
    
    % Calculate low and high cut-offs
    highHz =  min(high_cut_Hz, 0.499 * sample_rate);
    lowHz = max(0, min(low_cut_Hz, 0.499 * highHz));
    filter_band = [lowHz highHz] * 2 / sample_rate;
    filter_pass = 'bandpass';
    if ~isempty(varargin)
        switch varargin{1}
            case 'stop'
                filter_pass = 'stop';
            case 'low' % low pass frequencies below high cut
                filter_band = filter_band(2);
                filter_pass = 'low';
            case 'high' % high pass frequencies above low cut
                filter_band = filter_band(1);
                filter_pass = 'high';
            otherwise
        end
    end
    
    % Bandpass and high pass filters remove mean offset.
    data_mean = 0;
    if any(strcmp({'bandpass', 'high'}, filter_pass))
        data_mean = mean(data);
    end
    
    % Calculate and apply filter
    switch filter_index
        case 1
            [b, a] = butter(filter_order, filter_band, filter_pass);
        case 2
            [b, a] = cheby2(filter_order, 20, filter_band, filter_pass);
        case 3
            [b, a] = cheby1(filter_order, 0.5, filter_band,filter_pass);
    end  
    fdata = filter(b, a, data - data_mean);
end

