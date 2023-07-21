function handles = avgboundbox(handles)

load avgboundboxres.mat

handles.bands{1,1} = resboxi(1):resboxi(3);
handles.mxlocs = resboxi(2);