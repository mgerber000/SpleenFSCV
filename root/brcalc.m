function obj = brcalc(obj)

isgood = false;
minpeakprom = 0.01;
while isgood == false
    obj = getresprate(obj, minpeakprom, 100, true, true); uiwait;
    minpeakprom = inputdlg({'Try new peak boundary? Cancel to continue'},'Choose minpeakprom',1,{num2str(minpeakprom)});
    if ~isempty(minpeakprom)
        minpeakprom = str2double(minpeakprom{1,1});
    else
        isgood = true;
    end
end
        