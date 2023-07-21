function plot_trace(io,bounds)

plot(io,'k')
hold on
plot([1,length(io)],[0,0],'LineStyle','--')

if ~isnan(bounds(1)) && ~isnan(bounds(2))
    if ~isempty(find(isnan(io),1))
        idx = find(isnan(io));
        segs = sort([1,find(diff(idx) ~= 1),find(diff(idx) ~= 1) + 1,length(idx)]);
        segs = [segs(1:2:length(segs))',segs(2:2:length(segs))'];
        for i = 1:size(segs,1)
            x = [idx(segs(i,1)) - 1,idx(segs(i,2)) + 1];
            y = io(x);
            t = [0,length(x(1):x(2))-1];
            m = (y(2) - y(1))/(t(2) - t(1));
            b = y(1);
            inter = m*[1,2,3,4,5] + b;
            io(x(1) + 1:x(2) - 1) = inter;
        end
    end
    t = bounds(1):bounds(2);
    patch([t,flip(t,2)],...
        [io(t),zeros(1,length(t))],'k','FaceAlpha',0.5,'EdgeColor','none')
end

xlim([1,length(io)])