function varargout = PhysioFCSV_Reg(varargin)
% PHYSIOFCSV_REG MATLAB code for PhysioFCSV_Reg.fig
%      PHYSIOFCSV_REG, by itself, creates a new PHYSIOFCSV_REG or raises the existing
%      singleton*.
%
%      H = PHYSIOFCSV_REG returns the handle to a new PHYSIOFCSV_REG or the handle to
%      the existing singleton*.
%
%      PHYSIOFCSV_REG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PHYSIOFCSV_REG.M with the given input arguments.
%
%      PHYSIOFCSV_REG('Property','Value',...) creates a new PHYSIOFCSV_REG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PhysioFCSV_Reg_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PhysioFCSV_Reg_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PhysioFCSV_Reg

% Last Modified by GUIDE v2.5 03-Nov-2020 13:44:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PhysioFCSV_Reg_OpeningFcn, ...
                   'gui_OutputFcn',  @PhysioFCSV_Reg_OutputFcn, ...
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


% --- Executes just before PhysioFCSV_Reg is made visible.
function PhysioFCSV_Reg_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PhysioFCSV_Reg (see VARARGIN)

% Choose default command line output for PhysioFCSV_Reg
handles.data = varargin{1};
handles.listbox1.String = handles.data.fcsvdata.name;
handles.BslFile.String = handles.data.fcsvdata.name;
if handles.data.usephysio == true
    axes(handles.axes1)
    plot(1:size(handles.data.physdata.trigchan,1),handles.data.physdata.trigchan);
    trigoptions = num2cell(1:size(handles.data.physdata.triggers,1));
    trigoptions{1,end+1} = 'None';
    handles.listbox2.String = trigoptions;
    handles.listbox2.Value = size(trigoptions,2);
else
    handles.listbox2.Enable = 'off';
end
handles.BslStart.String = 1;
handles.BslEnd.String = 4;
handles.TrgStart.String = 10;
handles.TrgEnd.String = 190;
handles.SegStart.String = 0;
handles.SegEnd.String = 550;
handles.Averaging.String = 1;
handles.Waveform.Value = 0;
handles.Colormap.Value = 1;
handles.Intensity.String = 50;
handles.PW.String = 500;
handles.Freq.String = 10;
handles.Stims.String = 80;
handles.output = cell(4,9);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PhysioFCSV_Reg wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PhysioFCSV_Reg_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
useBR = false;

for i = 1:size(handles.output,1)
    if ~isempty(handles.output{i,1})
%         entry = find(strcmp(handles.data.fcsvdata.name,handles.output{i,1}));
        entry = i;
        handles.data.fcsvdata.baseline{entry,1} = handles.output{i,5};
        handles.data.fcsvdata.cmap{entry,1} = handles.output{i,2};
        handles.data.fcsvdata.waveform{entry,1} = handles.output{i,8};
        if strcmp(handles.output{i,3},'None')
            handles.data.fcsvdata.stim{entry,1} = false;
        else
            handles.data.fcsvdata.stim{entry,1} = true;
        end
        
        if handles.data.usephysio == true
            trgstart = handles.output{i,9};
            stmstart = handles.output{i,4}(1,1) - handles.output{i,7}(1,1);
            fcsv_sr = handles.data.fcsvdata.sr;
            trigid = handles.data.physdata.trigid;
            train_sr = handles.data.physdata.channel_meta(trigid).dt(1,1);
            t_trg = trgstart * train_sr - train_sr;
            t_start = t_trg - (stmstart * fcsv_sr - fcsv_sr);
            seglength = (handles.output{i,7}(1,2)*fcsv_sr -fcsv_sr) - (handles.output{i,7}(1,1)*fcsv_sr -fcsv_sr);
            t_end = t_start + seglength;
            handles.data.fcsvdata.t_fcsv{entry,1} = t_start:fcsv_sr:t_end;
            handles.data.fcsvdata.stm_fcsv{entry,1} = [(t_start - t_trg):fcsv_sr:0,fcsv_sr:fcsv_sr:t_end - t_trg];

            ecgid = handles.data.physdata.ecgid;
            ecg_sr = handles.data.physdata.channel_meta(ecgid).dt(1,1);
            ecgchan = handles.data.physdata.ecgchan;
            hrchan = handles.data.physdata.hrchan;
            handles.data.fcsvdata.t_ecg{entry,1} = t_start:ecg_sr:t_end;
            handles.data.fcsvdata.stm_ecg{entry,1} = [(t_start - t_trg):ecg_sr:0,ecg_sr:ecg_sr:t_end - t_trg];
            ind_ecg = round((t_start+ecg_sr)/ecg_sr):round((t_end+ecg_sr)/ecg_sr);
            handles.data.fcsvdata.ecg{entry,1} = ecgchan(ind_ecg)';
            handles.data.fcsvdata.hr{entry,1} = hrchan(ind_ecg)';

            if useBR
                brid = find(strcmp({handles.data.physdata.channel_meta.name},'Nasal sensor'));
                brchan = handles.data.physdata.brchan;
                br_sr = handles.data.physdata.channel_meta(brid).dt(1,1);
                handles.data.fcsvdata.t_br{entry,1} = t_start:br_sr:t_end;
                handles.data.fcsvdata.stm_br{entry,1} = [(t_start - t_trg):br_sr:0,br_sr:br_sr:t_end - t_trg];
                ind_br = round((t_start+br_sr)/br_sr):round((t_end+br_sr)/br_sr);
                handles.data.fcsvdata.br{entry,1} = brchan(ind_br)';
            else
                handles.data.fcsvdata.t_br{entry,1} = [];
                handles.data.fcsvdata.br{entry,1} = [];
            end
        else
            stmstart = handles.output{i,4}(1,1) - handles.output{i,7}(1,1);
            fcsv_sr = handles.data.fcsvdata.sr;
            t_start = -(stmstart * fcsv_sr - fcsv_sr);
            seglength = (handles.output{i,7}(1,2)*fcsv_sr -fcsv_sr) - (handles.output{i,7}(1,1)*fcsv_sr -fcsv_sr);
            t_end = t_start + seglength;
            handles.data.fcsvdata.t_fcsv{entry,1} = t_start:fcsv_sr:t_end;
            handles.data.fcsvdata.stm_fcsv{entry,1} = [(t_start):fcsv_sr:0,fcsv_sr:fcsv_sr:t_end];
        end
    end
end
handles.data.fcsvdata.current = [];
if length(handles.data.fcsvdata.name) == 1 && length(handles.data.fcsvdata.cmap) > 1
    name = handles.data.fcsvdata.name;
    for i = 1:length(handles.data.fcsvdata.cmap)
        name(i,1) = {strcat(name{1,1},'_',num2str(i))};
    end
    voltage = handles.data.fcsvdata.voltage;
    voltage(1:length(handles.data.fcsvdata.cmap),1) = {voltage{1,1}};
    handles.data.fcsvdata.name = name;
    handles.data.fcsvdata.voltage = voltage;
end
varargout{1} = handles.data;
close(handles.figure1);


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String'));
handles.fcsvselect = contents{get(hObject,'Value')};
handles.BslFile.Value = get(hObject,'Value');
handles.bslselect = contents{get(hObject,'Value')};

fcsvselect = find(strcmp(handles.data.fcsvdata.name,handles.fcsvselect));
bslselect = find(strcmp(handles.data.fcsvdata.name,handles.bslselect));
I = handles.data.fcsvdata.current{fcsvselect,1};
V = handles.data.fcsvdata.voltage{fcsvselect,1};
sr = handles.data.fcsvdata.sr;

bsl(1) = str2num(handles.BslStart.String);
bsl(2) = str2num(handles.BslEnd.String);
I_bsl = handles.data.fcsvdata.current{bslselect,1};
temp = [];
for i = 1:size(I_bsl,3)
    temp = [temp, I_bsl(:,:,i)];
end
t = 0:sr:size(I_bsl,2)*size(I_bsl,3)*sr - sr;
t_bsl = [closest(bsl(1),t) closest(bsl(2),t)]; % baseline start/end indices
i_bsl = temp(:,t_bsl(1):t_bsl(2)); % all current values during baseline

t = 0:sr:size(I,2)*size(I,3)*sr - sr;
stm(1) = str2num(handles.TrgStart.String);
stm(2) = str2num(handles.TrgEnd.String);
t_stm = [closest(stm(1),t) closest(stm(2),t)];

seg(1) = str2num(handles.SegStart.String);
seg(2) = str2num(handles.SegEnd.String);
t_seg = [closest(seg(1),t) closest(seg(2),t)];

smpT = 50;

if handles.Colormap.Value == 1
    plottype = 'colormap';
elseif handles.Waveform.Value == 1
    plottype = 'waveform';
end
av = str2num(handles.Averaging.String);

axes(handles.axes1)
[~,~] = createPlot(handles,I,V,sr,i_bsl,t_stm,smpT,t_seg,plottype,av,[]);
hold off

% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox2.
function listbox2_Callback(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

contents = cellstr(get(hObject,'String'));
select = contents{get(hObject,'Value')};
if strcmp(select,'None')
    handles.trgselect = NaN;
    axes(handles.axes1)
    plot(1:size(handles.data.physdata.trigchan,1),handles.data.physdata.trigchan);
    axes(handles.axes2)
    plot(1:size(handles.data.physdata.hrchan,1),handles.data.physdata.hrchan);
    handles.Intensity.String = 'NaN';
    handles.Intensity.Enable = 'off';
    handles.PW.String = 'NaN';
    handles.PW.Enable = 'off';
    handles.Freq.String = 'NaN';
    handles.Freq.Enable = 'off';
    handles.Stims.String = 'NaN';
    handles.Stims.Enable = 'off';
%     handles.TrgStart.Enable = 'off';
%     handles.TrgEnd.Enable = 'off';
%     handles.TrgStart.String = 'NaN';
%     handles.TrgEnd.String = 'NaN';
    [x,~] = getpts(handles.axes1);
    trgstart = round(x(1));
    handles.t_trg = trgstart;
    axes(handles.axes1)
    hold on
    scatter(trgstart,handles.data.physdata.trigchan(trgstart),'r')
    hold off
    hrseg = gethrseg(handles,trgstart);
    axes(handles.axes2)
    plot(1:size(hrseg,2),hrseg); ylim([100,600])
    axes(handles.axes1)
else
    handles.trgselect = str2num(select);
    handles.Intensity.Enable = 'on';
    handles.PW.Enable = 'on';
    handles.Freq.Enable = 'on';
    handles.Stims.Enable = 'on';
%     handles.Intensity.String = 50;
%     handles.PW.String = 500;
%     handles.Freq.String = 10;
%     handles.Stims.String = 80;
    handles.TrgStart.Enable = 'on';
    handles.TrgEnd.Enable = 'on';
%     handles.TrgStart.String = '6';
%     handles.TrgEnd.String = '14';
    axes(handles.axes1)
    plot(1:size(handles.data.physdata.trigchan,1),handles.data.physdata.trigchan);
    hold on
    scatter(handles.data.physdata.triggers(handles.trgselect,1):handles.data.physdata.triggers(handles.trgselect,2),handles.data.physdata.trigchan(handles.data.physdata.triggers(handles.trgselect,1):handles.data.physdata.triggers(handles.trgselect,2)));
    hold off
    
    trgstart = handles.data.physdata.triggers(handles.trgselect,1);
    handles.t_trg = trgstart;
    hrseg = gethrseg(handles,trgstart);
    axes(handles.axes2)
    plot(1:size(hrseg,2),hrseg); ylim([100,600])
    axes(handles.axes1); 
end
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function listbox2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TrgStart_Callback(hObject, eventdata, handles)
% hObject    handle to TrgStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TrgStart as text
%        str2double(get(hObject,'String')) returns contents of TrgStart as a double


% --- Executes during object creation, after setting all properties.
function TrgStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrgStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Remove.
function Remove_Callback(hObject, eventdata, handles)
% hObject    handle to Remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.uitable1.Data(handles.uitable1.UserData,:) = [];
handles.output(handles.uitable1.UserData,:) = [];
guidata(hObject, handles);

% --- Executes on selection change in BslFile.
function BslFile_Callback(hObject, eventdata, handles)
% hObject    handle to BslFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BslFile contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BslFile


% --- Executes during object creation, after setting all properties.
function BslFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BslFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BslStart_Callback(hObject, eventdata, handles)
% hObject    handle to BslStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BslStart as text
%        str2double(get(hObject,'String')) returns contents of BslStart as a double


% --- Executes during object creation, after setting all properties.
function BslStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BslStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function BslEnd_Callback(hObject, eventdata, handles)
% hObject    handle to BslEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BslEnd as text
%        str2double(get(hObject,'String')) returns contents of BslEnd as a double


% --- Executes during object creation, after setting all properties.
function BslEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BslEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Add.
function Add_Callback(hObject, eventdata, handles)
% hObject    handle to Add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fcsvselect = handles.listbox1.String{handles.listbox1.Value,1};
trgselect = handles.listbox2.String{handles.listbox2.Value,1};
if handles.data.usephysio == true
    t_trg = handles.t_trg;
else
    t_trg = 1;
end
stm = [str2double(handles.TrgStart.String),str2double(handles.TrgEnd.String)];
bsl = [str2double(handles.BslStart.String),str2double(handles.BslEnd.String)];
seg = [str2double(handles.SegStart.String),str2double(handles.SegEnd.String)];
bslfile = handles.BslFile.String{handles.BslFile.Value,1};
intensity = str2double(handles.Intensity.String);
pw = str2double(handles.PW.String);
freq = str2double(handles.Freq.String);
stims = str2double(handles.Stims.String);

I = handles.data.fcsvdata.current{find(strcmp(handles.data.fcsvdata.name,fcsvselect)),1};
V = handles.data.fcsvdata.voltage{find(strcmp(handles.data.fcsvdata.name,fcsvselect)),1};
sr = handles.data.fcsvdata.sr;

I_bsl = handles.data.fcsvdata.current{find(strcmp(handles.data.fcsvdata.name,bslfile)),1};
temp = [];
for i = 1:size(I_bsl,3)
    temp = [temp, I_bsl(:,:,i)];
end
t = 0:sr:size(I_bsl,2)*size(I_bsl,3)*sr - sr;
t_bsl = [closest(bsl(1),t) closest(bsl(2),t)]; % baseline start/end indices
i_bsl = temp(:,t_bsl(1):t_bsl(2)); % all current values during baseline

t = 0:sr:size(I,2)*size(I,3)*sr - sr;
t_stm = [closest(stm(1),t) closest(stm(2),t)];

t_seg = [closest(seg(1),t) closest(seg(2),t)];

av = str2num(handles.Averaging.String);

smpT = 50;

[handles,cmap] = createPlot(handles,I,V,sr,i_bsl,t_stm,smpT,t_seg,'none',av,[]);

if size(handles.uitable1.Data,1) == 4
    for i = 1:size(handles.uitable1.Data(:,1),1)
        temp(i) = isempty(handles.uitable1.Data{i,1});
    end
    emptypos = find(temp,1,'first');
    if isempty(emptypos)
        emptypos = 5;
    end
    handles.uitable1.Data{emptypos,1} = fcsvselect;
    handles.uitable1.Data{emptypos,2} = trgselect;
    handles.uitable1.Data{emptypos,3} = mat2str(stm);
    handles.uitable1.Data{emptypos,4} = mat2str(bsl);
    handles.uitable1.Data{emptypos,5} = bslfile;
    handles.uitable1.Data{emptypos,6} = mat2str(seg);
    handles.uitable1.Data{emptypos,7} = mat2str([intensity,pw,freq,stims]);
    
    handles.output{emptypos,1} = fcsvselect;
    handles.output{emptypos,2} = cmap;
    handles.output{emptypos,3} = trgselect;
    handles.output{emptypos,4} = t_stm;
    handles.output{emptypos,5} = i_bsl;
    handles.output{emptypos,6} = bslfile;
    handles.output{emptypos,7} = t_seg;
    handles.output{emptypos,8} = [intensity,pw,freq,stims];
    handles.output{emptypos,9} = t_trg;
else
    emptypos = size(handles.uitable1.Data,1) + 1;
    handles.uitable1.Data{emptypos,1} = fcsvselect;
    handles.uitable1.Data{emptypos,2} = trgselect;
    handles.uitable1.Data{emptypos,3} = mat2str(stm);
    handles.uitable1.Data{emptypos,4} = mat2str(bsl);
    handles.uitable1.Data{emptypos,5} = bslfile;
    handles.uitable1.Data{emptypos,6} = mat2str(seg);
    handles.uitable1.Data{emptypos,7} = mat2str([intensity,pw,freq,stims]);
    
    handles.output{emptypos,1} = fcsvselect;
    handles.output{emptypos,2} = cmap;
    handles.output{emptypos,3} = trgselect;
    handles.output{emptypos,4} = t_stm;
    handles.output{emptypos,5} = i_bsl;
    handles.output{emptypos,6} = bslfile;
    handles.output{emptypos,7} = t_seg;
    handles.output{emptypos,8} = [intensity,pw,freq,stims];
    handles.output{emptypos,9} = t_trg;
end

guidata(hObject, handles);


% --- Executes on button press in Confirm.
function Confirm_Callback(hObject, eventdata, handles)
% hObject    handle to Confirm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);


function TrgEnd_Callback(hObject, eventdata, handles)
% hObject    handle to TrgEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TrgEnd as text
%        str2double(get(hObject,'String')) returns contents of TrgEnd as a double


% --- Executes during object creation, after setting all properties.
function TrgEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrgEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SegStart_Callback(hObject, eventdata, handles)
% hObject    handle to SegStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
segstart = str2double(get(hObject,'String'));
handles.BslStart.String = num2str(segstart + 1);
handles.BslEnd.String = num2str(segstart + 4);
handles.TrgStart.String = num2str(segstart + 10);
handles.TrgEnd.String = num2str(segstart + 30);
handles.SegStart.String = num2str(segstart);
handles.SegEnd.String = num2str(segstart + 180);

% Hints: get(hObject,'String') returns contents of SegStart as text
%        str2double(get(hObject,'String')) returns contents of SegStart as a double


% --- Executes during object creation, after setting all properties.
function SegStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SegStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SegEnd_Callback(hObject, eventdata, handles)
% hObject    handle to SegEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SegEnd as text
%        str2double(get(hObject,'String')) returns contents of SegEnd as a double


% --- Executes during object creation, after setting all properties.
function SegEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SegEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Averaging_Callback(hObject, eventdata, handles)
% hObject    handle to Averaging (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Averaging as text
%        str2double(get(hObject,'String')) returns contents of Averaging as a double


% --- Executes during object creation, after setting all properties.
function Averaging_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Averaging (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function uitable1_CellSelectionCallback(hObject, eventdata,handles)

if ~isempty(eventdata.Indices)
    handles.uitable1.UserData = eventdata.Indices(1,1);
end



function Intensity_Callback(hObject, eventdata, handles)
% hObject    handle to Intensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Intensity as text
%        str2double(get(hObject,'String')) returns contents of Intensity as a double


% --- Executes during object creation, after setting all properties.
function Intensity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Intensity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function PW_Callback(hObject, eventdata, handles)
% hObject    handle to PW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of PW as text
%        str2double(get(hObject,'String')) returns contents of PW as a double


% --- Executes during object creation, after setting all properties.
function PW_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PW (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Freq_Callback(hObject, eventdata, handles)
% hObject    handle to Freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Freq as text
%        str2double(get(hObject,'String')) returns contents of Freq as a double


% --- Executes during object creation, after setting all properties.
function Freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Stims_Callback(hObject, eventdata, handles)
% hObject    handle to Stims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Stims as text
%        str2double(get(hObject,'String')) returns contents of Stims as a double


% --- Executes during object creation, after setting all properties.
function Stims_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Stims (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
