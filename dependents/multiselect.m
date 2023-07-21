function vals = multiselect(entry)
fig = uifigure('Position',[100 100 350 275]);


% Create Text Area
txt = uitextarea(fig,...
    'Position',[125 80 100 50]);

% Create List Box
lbox = uilistbox(fig,...
    'Position',[125 150 100 78],...
    'Multiselect','on',...
    'Items',entry,...
    'ValueChangedFcn',@selectionChanged);

btn = uibutton(fig,'push','Position',[125,40,80,40],...
    'ButtonPushedFcn',@(~,~)pushBtOk);

uiwait(fig);
% ValueChangedFcn callback
function selectionChanged(src,event)
    txt.Value = src.Value;
end

function pushBtOk(src,event)
    vals = find(contains(entry,txt.Value));
    close(fig)
end

end