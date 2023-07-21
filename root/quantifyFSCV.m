function [maxsig,tmax,NEval,dtrace] = quantifyFSCV(cmap,bands,timebounds,mxband)

if ~isnan(timebounds(1)) && ~isnan(timebounds(2))
    t = timebounds(1):timebounds(2);
    tmax = t(cmap(mxband,t) == max(cmap(mxband,t)));
    maxsig = cmap(mxband,tmax);

    fs = 1/10000;
    for i = 1:size(cmap,2)
        point1 = cmap(bands(1),i);
        mxi = cmap(mxband,i);
        dtrace(i) = fscv_peak2trough(point1,mxi);
    end
    sig = dtrace(t);
    x = 0.1*t(isfinite(sig));
    y = sig(isfinite(sig));
    NEval = trapz(x,y);
else
    tmax = NaN;
    maxsig = NaN;
    
    fs = 1/10000;
    for i = 1:size(cmap,2)
        point1 = cmap(bands(1),i);
        mxi = cmap(mxband,i);
        dtrace(i) = fscv_peak2trough(point1,mxi);
    end
    NEval = NaN;
end