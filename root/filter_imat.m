function imat = filter_imat(fcsvdata,filtswitch)
fscv = fcsvdata.current{1,1};

imat = [];
for i = 1:size(fscv,3)
    imat = [imat,fscv(:,:,i)];
end

if filtswitch
    d = designfilt('highpassiir', ...       % Response type
           'StopbandFrequency',0.0005, ...     % Frequency constraints
           'PassbandFrequency',0.005, ...
           'StopbandAttenuation',55, ...    % Magnitude constraints
           'PassbandRipple',4, ...
           'DesignMethod','butter', ...     % Design method
           'SampleRate',10);               % Sample rate

    for i = 1:size(imat,1)
        f = isfinite(imat(i,:));
        x = zeros(1,length(imat(i,:)))*NaN;
        y = imat(i,f);
    %     periodogram(y,rectwin(length(y)),length(y),10)
    %     hold on
        dy = filtfilt(d,y);
        x(f) = dy;
        imat(i,:) = x;
    end
end

% L = size(fscv,2);
% dimat = zeros(size(fscv,1),size(fscv,2),size(fscv,3));
% for i = 0:size(fscv,3)-1
%     dimat(:,:,i+1) = imat(:,i*L + 1:(i+1)*L);
% end
% 
% fcsvdata.current{1,1} = dimat;