function varargout = FCSVCompare(varargin)
% FCSVCOMPARE MATLAB code for FCSVCompare.fig
%      FCSVCOMPARE, by itself, creates a new FCSVCOMPARE or raises the existing
%      singleton*.
%
%      H = FCSVCOMPARE returns the handle to a new FCSVCOMPARE or the handle to
%      the existing singleton*.
%
%      FCSVCOMPARE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FCSVCOMPARE.M with the given input arguments.
%
%      FCSVCOMPARE('Property','Value',...) creates a new FCSVCOMPARE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FCSVCompare_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FCSVCompare_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FCSVCompare

% Last Modified by GUIDE v2.5 19-Oct-2021 18:20:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FCSVCompare_OpeningFcn, ...
                   'gui_OutputFcn',  @FCSVCompare_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before FCSVCompare is made visible.
function FCSVCompare_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FCSVCompare (see VARARGIN)

% Choose default command line output for FCSVCompare
handles.output = hObject;
handles.data = varargin{1,1};
handles.output = varargin{1,2};
if isempty(handles.output)
    handles.output = table('Size',[1,12],'VariableTypes',{'string','string','double','double','double','double','double','double','double','double','double','double'}, ...
        'VariableNames',{'DirName','FileName','DroppedBeats','NE_vol','tbound1','tbound2','vbound1','maxsigloc','vbound2','absolutepeakA','AbsolutepeakT','AUC'});
end
handles.bands = [];
handles.mxlocs = [];
handles.bounds = [];
j = 1;
for i = 1:size(handles.data.fcsvdata.stim,1)
    if ~isempty(handles.data.fcsvdata.stim{i,1})
        seglength(j) = size(handles.data.fcsvdata.cmap{i,1},2);
        j = j+1;
    end
end

j = 1;
for i = 1:size(handles.data.fcsvdata.stim,1)
    if ~isempty(handles.data.fcsvdata.stim{i,1})
        names{j,1} = handles.data.fcsvdata.name{i,1};
        data(:,:,j) = handles.data.fcsvdata.cmap{i,1}(:,1:min(seglength));
        j = j+1;
    end
end
handles.SelectFile.String = names;
data = reshape(data,1,numel(data));
handles.Maximum.String = max(data);
handles.Minimum.String = min(data);
handles.HighPrctile.String = prctile(data,99.9);
handles.LowPrctile.String = prctile(data,0.1);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FCSVCompare wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FCSVCompare_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
close(handles.figure1)


% --- Executes on selection change in SelectFile.
function SelectFile_Callback(hObject, eventdata, handles)
% hObject    handle to SelectFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String'));
select = find(strcmp(handles.data.fcsvdata.name,handles.SelectFile.String{handles.SelectFile.Value,1}));
handles.select = select;
cmap = handles.data.fcsvdata.cmap{select,1};
V = handles.data.fcsvdata.voltage{1,1};
% if handles.data.usephysio
    t = handles.data.fcsvdata.t_fcsv{select,1};
    t = 1:length(t);
    stm = handles.data.fcsvdata.stm_fcsv{select,1};
% else
%     t = 1:size(cmap,2);
%     stm = t - 100;
%     handles.data.fcsvdata.stm_fcsv{select,1} = stm;
% end
lo_edge = str2num(handles.LowPrctile.String);
up_edge = str2num(handles.HighPrctile.String);
Vi = 1:length(V); Vt = floor(linspace(1,length(V),6));
colormap hsv
axes(handles.axes1);
imagesc(t,Vi,cmap,[lo_edge up_edge]);
colorbar
Vtl = num2cell(ndecp(V(Vt),2));
set(gca,'YDir','normal','YTick',[0 200 400 600 800 1000],'YTickLabel',Vtl)
xlabel ('time (s)')
ylabel ('voltage (V)')
if handles.data.fcsvdata.stim{select,1} == 1
    mxy = max(get(gca,'YLim')); mny = min(get(gca,'YLim'));
    t_stm(1) = t(closest(0,stm)); t_stm(2) = t(closest(8,stm));
    xp_stm = [t_stm(1) t_stm(2) t_stm(2) t_stm(1)];
    yp = [mny mny mxy mxy];
    patch(xp_stm,yp,[1 0 0],'FaceAlpha',0.3)
end

waveform = handles.data.fcsvdata.waveform{select,1};
str = strcat(num2str(waveform(1)),'uA, ',num2str(waveform(2)),'pw, ',num2str(waveform(3)),'Hz, ',num2str(waveform(4)),'stims');
title(str)

% if handles.data.usephysio == 1
%     ecg = handles.data.fcsvdata.hr{select,1};
%     t_ecg = handles.data.fcsvdata.t_ecg{select,1};
%     stm = handles.data.fcsvdata.stm_ecg{select,1};
%     axes(handles.ECGaxis)
%     plot(t_ecg,ecg)
%     xlim([t_ecg(1),t_ecg(end)])
%     if handles.data.fcsvdata.stim{select,1} == 1
%         mxy = max(get(gca,'YLim')); mny = min(get(gca,'YLim'));
%         t_stm(1) = t_ecg(closest(0,stm)); t_stm(2) = t_ecg(closest(8,stm));
%         xp_stm = [t_stm(1) t_stm(2) t_stm(2) t_stm(1)];
%         yp = [mny mny mxy mxy];
%         patch(xp_stm,yp,[1 0 0],'FaceAlpha',0.3)
%         ylim([mny,mxy])
%     end
% 
%     br = handles.data.fcsvdata.br{select,1};
%     t_br = handles.data.fcsvdata.t_br{select,1};
%     if ~isempty(br)
%         stm = handles.data.fcsvdata.stm_br{select,1};
%         axes(handles.BRaxis)
%         plot(t_br,br)
%         xlim([t_br(1),t_br(end)])
%         xlabel('Time (s)');
%         if handles.data.fcsvdata.stim{select,1} == 1
%             mxy = max(get(gca,'YLim')); mny = min(get(gca,'YLim'));
%             t_stm(1) = t_br(closest(0,stm)); t_stm(2) = t_br(closest(8,stm));
%             xp_stm = [t_stm(1) t_stm(2) t_stm(2) t_stm(1)];
%             yp = [mny mny mxy mxy];
%             patch(xp_stm,yp,[1 0 0],'FaceAlpha',0.3)
%             ylim([mny,mxy])
%         end
%     end
% else
%     axes(handles.ECGaxis)
%     plot(ones(1,size(cmap,2)))
%     
%     axes(handles.BRaxis)
%     plot(ones(1,size(cmap,2)))
% end

guidata(hObject, handles);

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% handles = avgboundbox(handles);
if isempty(handles.bands)
    for i = 1:size(handles.data.fcsvdata.cmap,1)
        if ~isempty(handles.data.fcsvdata.cmap{i,1})
            [~,~,maxsig] = bounddetect(handles.data,i);
            sigcorrmax(i) = maxsig;
        end
    end
    maxcmap = find(sigcorrmax == max(sigcorrmax));
    [lowbound,highbound,~] = bounddetect(handles.data,maxcmap);
    wv = mean(handles.data.fcsvdata.cmap{maxcmap,1}(:,lowbound:highbound),2);
    load wvfilt.mat
    fwv = filtfilt(wvD,wv);
    [handles.bands,handles.mxlocs] = dynamicband(fwv);
end

if handles.radiobutton1.Value
    if ~isempty(handles.bounds)
        lowbound = handles.bounds(1,1);
        highbound = handles.bounds(1,2);
        handles.radiobutton1.Value = false;
    else
        return
    end
else
    [lowbound,highbound,~] = bounddetect(handles.data,handles.select);
    %[lowbound,highbound] = boundtune(handles.data,handles.select,lowbound,highbound,handles.mxlocs);
end

[dbeats,sigband,maxsig,maxloc,tlowbound,thighbound,NEval,dtrace] = quantifyFSCV(handles.data,handles.select,handles.bands,lowbound:highbound,handles.mxlocs);
lowbound = tlowbound;
highbound = thighbound;
temp{1,1} = handles.data.inputdir;
temp{1,2} = handles.data.fcsvdata.name{handles.select,1};
temp{1,3} = {dbeats};
temp{1,4} = {dtrace};
temp{1,5} = {lowbound};
temp{1,6} = {highbound};
v = handles.data.fcsvdata.voltage{handles.select,1};
temp{1,7} = {v(handles.bands{1,1}(1,1),1)};
temp{1,8} = {v(handles.mxlocs,1)};
temp{1,9} = {v(handles.bands{1,1}(1,end),1)};
% temp{1,7} = {handles.bands{1,1}(1,1)};
% temp{1,8} = {handles.mxlocs};
% temp{1,9} = {handles.bands{1,1}(1,end)};
temp{1,10} = {maxsig};
temp{1,11} = {maxloc};
temp{1,12} = {NEval};
handles.output = [handles.output;temp];
handles.bounds(1,1) = lowbound;
handles.bounds(1,2) = highbound;

axes(handles.axes1);
hold on
timebounds = lowbound:highbound;
if ~isempty(lowbound)
    scatter(lowbound*ones(1000,1),1:1000,1,'k','filled');
end
if ~isempty(highbound)
    scatter(highbound*ones(1000,1),1:1000,1,'k','filled');
end
if ~isempty(lowbound) && ~isempty(highbound)
    scatter(timebounds,handles.bands{1,1}(1,1)*ones(size(timebounds,2),1),1,'k','filled');
    scatter(timebounds,handles.bands{1,1}(1,end)*ones(size(timebounds,2),1),1,'k','filled');
%     scatter(timebounds,handles.bands{1,2}(1,end)*ones(size(timebounds,2),1),1,'k','filled');
    scatter(timebounds,handles.mxlocs*ones(size(timebounds,2),1),1,'k','filled');
%     scatter(timebounds,handles.bands{1,3}(1,end)*ones(size(timebounds,2),1),1,'k','filled');
end
hold off

% axes(handles.BRaxis)
% plot(dtrace)

axes(handles.ECGaxis)
plot((1:length(dtrace))/10,dtrace); xlim([1,length(dtrace)/10])

guidata(hObject, handles);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1)

% --- Executes during object creation, after setting all properties.
function SelectFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SelectFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function MindistECG_Callback(hObject, eventdata, handles)
% hObject    handle to MindistECG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function MindistECG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MindistECG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function PeakECG_Callback(hObject, eventdata, handles)
% hObject    handle to PeakECG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function PeakECG_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PeakECG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function MindistBR_Callback(hObject, eventdata, handles)
% hObject    handle to MindistBR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function MindistBR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MindistBR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function PeakBR_Callback(hObject, eventdata, handles)
% hObject    handle to PeakBR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function PeakBR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PeakBR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function LowPrctile_Callback(hObject, eventdata, handles)
% hObject    handle to LowPrctile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LowPrctile as text
%        str2double(get(hObject,'String')) returns contents of LowPrctile as a double


% --- Executes during object creation, after setting all properties.
function LowPrctile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LowPrctile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function HighPrctile_Callback(hObject, eventdata, handles)
% hObject    handle to HighPrctile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HighPrctile as text
%        str2double(get(hObject,'String')) returns contents of HighPrctile as a double


% --- Executes during object creation, after setting all properties.
function HighPrctile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HighPrctile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Minimum_Callback(hObject, eventdata, handles)
% hObject    handle to Minimum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Minimum as text
%        str2double(get(hObject,'String')) returns contents of Minimum as a double


% --- Executes during object creation, after setting all properties.
function Minimum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Minimum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Maximum_Callback(hObject, eventdata, handles)
% hObject    handle to Maximum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Maximum as text
%        str2double(get(hObject,'String')) returns contents of Maximum as a double


% --- Executes during object creation, after setting all properties.
function Maximum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Maximum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
