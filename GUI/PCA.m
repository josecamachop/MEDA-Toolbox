function varargout = PCA(varargin)

%%
%GUI for Principal Components Analysis (PCA) analysis
%
%This M-file include routines from the EDA Toolbox: 
%loading_pca.m, meda_pca.m, omeda_pca.m, pca_pp.m, scores_pca.m,
%sqresiduals_pca.m and var_pca.m
%
% coded by: Elena Jim�nez Ma�as (elenajm@correo.ugr.es).
% version: 2.0
% last modification: 07/Jul/14.
%
% Copyright (C) 2014  Elena Jim�nez Ma�as
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% PCA M-file for PCA.fig
%      PCA, by itself, creates a new PCA or raises the existing
%      singleton*.
%
%      H = PCA returns the handle to a new PCA or the handle to
%      the existing singleton*.
%
%      PCA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PCA.M with the given input arguments.
%
%      PCA('Property','Value',...) creates a new PCA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PCA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PCA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PCA

% Last Modified by GUIDE v2.5 03-May-2013 19:18:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @PCA_OpeningFcn, ...
    'gui_OutputFcn',  @PCA_OutputFcn, ...
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


% --- Executes just before PCA is made visible.
function PCA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PCA (see VARARGIN)

% Choose default command line output for PCA
handles.output = hObject;

%Definici�n del estado inicial de la interfaz gr�fica PCA:

%Score plot
set(handles.text7,'Enable','off');
set(handles.text8,'Enable','off');
set(handles.popupmenu1,'Enable','off');
set(handles.popupmenu7,'Enable','off');
set(handles.text13,'Enable','off');
set(handles.text14,'Enable','off');
set(handles.popupmenu15,'Enable','off');
set(handles.popupmenu16,'Enable','off');
set(handles.pushbutton4,'Enable','off');

%MEDA
set(handles.popupmenu9,'Enable','off');
set(handles.popupmenu9,'String',' ');
set(handles.text5,'Enable','off');
set(handles.edit6,'Enable','off');
set(handles.radiobutton1,'Enable','off');
set(handles.radiobutton2,'Enable','off');
set(handles.pushbutton5,'Enable','off');
set(handles.pushbutton22,'Enable','off');

%oMEDA
set(handles.pushbutton6,'Enable','off');
set(handles.pushbutton7,'Enable','off');
set(handles.pushbutton8,'Enable','off');
set(handles.pushbutton9,'Enable','off');
set(handles.pushbutton10,'Enable','off');
set(handles.pushbutton11,'Enable','off');

%Loading plot
set(handles.text9,'Enable','off');
set(handles.text10,'Enable','off');
set(handles.popupmenu10,'Enable','off');
set(handles.popupmenu11,'Enable','off');
set(handles.text17,'Enable','off');
set(handles.text18,'Enable','off');
set(handles.popupmenu19,'Enable','off');
set(handles.popupmenu20,'Enable','off');
set(handles.pushbutton12,'Enable','off');

%Residue
set(handles.pushbutton20,'Enable','off');
set(handles.pushbutton21,'Enable','off');

%Information Panel:
handles.data.text=[];
text=sprintf('To begin the analysis:\nchoose a data matrix from the data popupmenu. If there is no data, charge data from WorkSpace by clicking on REFRESH button.');
handles.data.text=cprint(handles.text1,text,handles.data.text,0);

%Variables initialization:
handles.data.PCs=[];
handles.data.PC1=[];
handles.data.PC2=[];
handles.data.weightDummy=cell(1,1000);
handles.data.dummy={};
handles.data.sp_ID_figures=[];
handles.data.sp_matrix={};
handles.data.clean_control=zeros(1,1000);
handles.data.lp_ID_figures=[];
handles.data.lp_matrix={};
handles.data.PC1_LP=[];
handles.data.PC2_LP=[];
handles.data.control_Refresh=0;
handles.data.CORTES={};
handles.data.matrix_2PCs={};
handles.data.PCs_MEDA='';
handles.data.auxPCs=0;

%Change icon
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
javaFrame = get(hObject,'JavaFrame');
javaFrame.setFigureIcon(javax.swing.ImageIcon('icon.jpg'));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PCA wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PCA_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PCA Analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in popupmenu6.
%popupmenu6==Data
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6

incoming_data=get(hObject,'Value');%Incoming data position
string_evaluation=handles.data.WorkSpace{incoming_data};%Name of the incoming data position
data_matrix=evalin('base',string_evaluation);%Data content in that name
handles.data.data_matrix=data_matrix;

%Information Panel:
text=sprintf('Enter the number of principal components to work with, using this format: 1:x , being x the number of principal components.');
handles.data.text=cprint(handles.text1,text,handles.data.text,0);

%Initialize dummy variable:
M=size(data_matrix,1);%Number of observations
dummy=zeros(1,M);
handles.data.dummyRED=dummy;
handles.data.dummyGREEN=dummy;

set(handles.popupmenu16,'Value',1);
handles.data.label={};
set(handles.popupmenu15,'Value',1);
handles.data.classes=[];
set(handles.popupmenu20,'Value',1);
handles.data.label_LP={};
set(handles.popupmenu19,'Value',1);
handles.data.classes_LP=[];

handles.data.namePopupmenu6=string_evaluation;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
handles.data.namePopupmenu6='';
handles.data.data_matrix=[];
handles.data.WorkSpace=evalin('base','who');%name of the variables in the workspace

if ~isempty(handles.data.WorkSpace),
    set(hObject,'String',handles.data.WorkSpace);
    string_evaluation=handles.data.WorkSpace{1};%Name of the incoming data position
    data_matrix=evalin('base',string_evaluation);%Data content in that name
    handles.data.data_matrix=data_matrix;
    handles.data.namePopupmenu6=string_evaluation;
    %Initialize dummy variable:
    M=size(data_matrix,1);%Number of observations
    dummy=zeros(1,M);
    handles.data.dummyRED=dummy;
    handles.data.dummyGREEN=dummy;
else
    set(hObject,'String',' ');
end

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

guidata(hObject, handles);


% --- Executes on button press in pushbutton1.
%pushbutton1==Refresh
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data.WorkSpace=evalin('base','who');

if ~isempty(handles.data.WorkSpace),
    
    set(handles.popupmenu6, 'String', handles.data.WorkSpace);
    nombres=cellstr(get(handles.popupmenu6,'String'));
    if ~isempty(handles.data.namePopupmenu6),
        for i=1:length(nombres),
            if strcmp(nombres(i),handles.data.namePopupmenu6),
                val=i;
            end
        end
        set(handles.popupmenu6,'Value',val);
        handles.data.data_matrix=evalin('base',handles.data.WorkSpace{val});
    end
    %Para que la primera vez que se pulse Refresh con el workspace distinto
    %de vacio coja la primera matriz automaticamente
    if handles.data.control_Refresh==0 && isempty(handles.data.data_matrix),
        string_evaluation=handles.data.WorkSpace{1};
        data_matrix=evalin('base',string_evaluation);
        handles.data.data_matrix=data_matrix;
        handles.data.namePopupmenu6=string_evaluation;
    end
    
    %Refresh the Label and Classes popupmenus:
    contents=get(handles.popupmenu15,'String');
    aux=[];
    for i=1:length(handles.data.WorkSpace),
        aux=[aux handles.data.WorkSpace(i,:)];
    end
    a1=contents(1,:);
    for j=1:length(a1),
        if ~isspace(a1(j)),
            b1(j)=a1(j);
        end
    end
    aux=[b1,aux];
    set(handles.popupmenu15,'String',strvcat(aux));
    nombres=cellstr(get(handles.popupmenu15,'String'));
    if ~strcmp(handles.data.namePopupmenu15,'emptyclasses'),
        for i=1:length(nombres),
            if strcmp(nombres(i),handles.data.namePopupmenu15),
                val=i;
            end
        end
        set(handles.popupmenu15,'Value',val);
        handles.data.classes=evalin('base',handles.data.WorkSpace{val-1});    
    end
    
    contents=get(handles.popupmenu16,'String');
    aux2=[];
    for i=1:length(handles.data.WorkSpace),
        aux2=[aux2 handles.data.WorkSpace(i,:)];
    end
    a2=contents(1,:);
    for j=1:length(a2),
        if ~isspace(a2(j)),
            b2(j)=a2(j);
        end
    end
    aux2=[b2,aux2];
    set(handles.popupmenu16,'String',strvcat(aux2));
    nombres=cellstr(get(handles.popupmenu16,'String'));
    if ~strcmp(handles.data.namePopupmenu16,'emptylabel'),
        for i=1:length(nombres),
            if strcmp(nombres(i),handles.data.namePopupmenu16),
                val=i;
            end
        end
        set(handles.popupmenu16,'Value',val);
        handles.data.label=evalin('base',handles.data.WorkSpace{val-1});    
    end
    
    
    contents=get(handles.popupmenu19,'String');
    aux3=[];
    for i=1:length(handles.data.WorkSpace),
        aux3=[aux3 handles.data.WorkSpace(i,:)];
    end
    a3=contents(1,:);
    for j=1:length(a3),
        if ~isspace(a3(j)),
            b3(j)=a3(j);
        end
    end
    aux3=[b3,aux3];
    set(handles.popupmenu19,'String',strvcat(aux3));
    nombres=cellstr(get(handles.popupmenu19,'String'));
    if ~strcmp(handles.data.namePopupmenu19,'emptyclasses'),
        for i=1:length(nombres),
            if strcmp(nombres(i),handles.data.namePopupmenu19),
                val=i;
            end
        end
        set(handles.popupmenu19,'Value',val);
        handles.data.classes_LP=evalin('base',handles.data.WorkSpace{val-1});    
    end
    
    contents=get(handles.popupmenu20,'String');
    aux4=[];
    for i=1:length(handles.data.WorkSpace),
        aux4=[aux4 handles.data.WorkSpace(i,:)];
    end
    a4=contents(1,:);
    for j=1:length(a4),
        if ~isspace(a4(j)),
            b4(j)=a4(j);
        end
    end
    aux4=[b4,aux4];
    set(handles.popupmenu20,'String',strvcat(aux4));
    nombres=cellstr(get(handles.popupmenu20,'String'));
    if ~strcmp(handles.data.namePopupmenu20,'emptylabel'),
        for i=1:length(nombres),
            if strcmp(nombres(i),handles.data.namePopupmenu20),
                val=i;
            end
        end
        set(handles.popupmenu20,'Value',val);
        handles.data.label_LP=evalin('base',handles.data.WorkSpace{val-1});    
    end
    
    %Information panel:
    text=sprintf('Select the new data matrix to work with from the Data popupmenu.');
    handles.data.text=cprint(handles.text1,text,handles.data.text,0);
    
    handles.data.control_Refresh=1;
else
    set(handles.popupmenu6, 'String', ' ');
    handles.data.data_matrix=[];
    
    contents=get(handles.popupmenu15,'String');
    aux=[];
    aux=[contents(1,:),aux];
    
    contents=get(handles.popupmenu16,'String');
    aux2=[];
    aux2=[contents(1,:),aux2];
    
    contents=get(handles.popupmenu19,'String');
    aux3=[];
    aux3=[contents(1,:),aux3];
    
    contents=get(handles.popupmenu20,'String');
    aux4=[];
    aux4=[contents(1,:),aux4];
        
    %Information panel:
    text=sprintf('Warning: No data matrices in workspace.');
    handles.data.text=cprint(handles.text1,text,handles.data.text,0);
end

handles.data.new2=aux;
handles.data.new1=aux2;
handles.data.new4=aux3;
handles.data.new3=aux4;
guidata(hObject,handles);

%edit text==PCs
function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double
PCs=str2num(get(hObject,'String'));

if max(PCs)>size(handles.data.data_matrix,2),
    errordlg(sprintf('The number of PCs can not exceed the number of variables in the data matrix which is %d.',size(handles.data.data_matrix,2)));
    set(hObject,'String',strvcat(sprintf('1:%d',size(handles.data.data_matrix,2))));
    PCs=str2num(get(hObject,'String'));
end

handles.data.PCs = PCs;

%Si la variable handles.data.PCs es distinta de vac�a, imprimir en popupmenu1,
%popupmenu10, popupmenu11 y popupmenu7 los PCs posibles.
if ~isempty(handles.data.PCs),
    set(handles.popupmenu1, 'String',handles.data.PCs);
    set(handles.popupmenu7, 'String',handles.data.PCs);
    set(handles.popupmenu10, 'String',handles.data.PCs);
    set(handles.popupmenu11, 'String',handles.data.PCs);
    
    %Imprimir en popupmenu de submenu MEDA todas las combinaciones posibles
    %para hacer MEDA
    k=min(handles.data.PCs);
    options=[];
    for i=min(handles.data.PCs):max(handles.data.PCs),
        for j=k:max(handles.data.PCs),
            options=[options,i,j];
        end
        k=k+1;
    end
    
    set(handles.popupmenu9,'String','');
    for i=1:2:(length(options)-1),
        contents=get(handles.popupmenu9,'String');
        set(handles.popupmenu9,'String',strvcat(contents,sprintf('%d:%d',options(i),options(i+1))));
    end
end

if handles.data.auxPCs==0,
handles.data.PC1=min(handles.data.PCs);
handles.data.PC2=min(handles.data.PCs);
handles.data.PC1_LP=min(handles.data.PCs);
handles.data.PC2_LP=min(handles.data.PCs);
handles.data.PCs_MEDA=sprintf('%d:%d',min(handles.data.PCs),min(handles.data.PCs));
handles.data.auxPCs=1;
end

%Information Panel:
text=sprintf('Select the preprocessing of the data:\n-No preprocessing\n-Mean centering (default)\n-Auto-scaling (centers and scales data so that each variable has variance 1).');
handles.data.text=cprint(handles.text1,text,handles.data.text,0);

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton23.
%pushbutton==VAR
function pushbutton23_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x_var = var_pca(handles.data.data_matrix,max(handles.data.PCs),handles.data.prep,1);

% --- Executes on selection change in popupmenu5.
%popupmenu5==Prep
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5
nombres=cellstr(get(hObject,'String'));
val=nombres{get(hObject,'Value')};

switch val,
    case 'no preprocessing',
        prep=0;
        if handles.data.control==1,
            text=sprintf('Selected the no preprocessing option.');
            handles.data.text=cprint(handles.text1,text,handles.data.text,0);
        end
    case 'mean centering',
        prep=1;
        if handles.data.control==1,
            text=sprintf('Selected the mean centering option.');
            handles.data.text=cprint(handles.text1,text,handles.data.text,0);
        end
    case 'auto-scaling',
        prep=2;
        if handles.data.control==1,
            text=sprintf('Selected the auto-scalling option.');
            handles.data.text=cprint(handles.text1,text,handles.data.text,0);
        end
end

handles.data.prep = prep;

if handles.data.control==1,
    text=sprintf('Press on the PCA button to perform the initial analysis and activate the Score Plot, Loading Plot and MEDA menus.');
    handles.data.text=cprint(handles.text1,text,handles.data.text,1);
end
handles.data.control=1;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
set(hObject,'String',strvcat('no preprocessing','mean centering','auto-scaling'));

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

handles.data.control=0;
set(hObject, 'Value', 2);%Default value for the preprocessing method: mean-centering
popupmenu5_Callback(hObject, eventdata, handles)%Para llamar al valor por defecto


% --- Executes on button press in pushbutton2.
%pushbutton2==PCA
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Information panel:
if isempty(handles.data.data_matrix),
    errordlg('No data matrix selected, please select one.');
    return;
end
if isempty(handles.data.PCs),
    errordlg('No PCs defined, please define them.');
    return;
end

[handles.data.matrixLoadings,handles.data.matrixScores]=pca_pp(handles.data.data_matrix,max(handles.data.PCs));

%Definici�n del estado de la interfaz tras pulsar PCA:
%Score plot
set(handles.popupmenu1,'Enable','on');
set(handles.popupmenu7,'Enable','on');
set(handles.pushbutton4,'Enable','on');
set(handles.text7,'Enable','on');
set(handles.text8,'Enable','on');
set(handles.text13,'Enable','on');
set(handles.text14,'Enable','on');
set(handles.popupmenu15,'Enable','on');
set(handles.popupmenu16,'Enable','on');

%MEDA
set(handles.radiobutton1,'Enable','on');
set(handles.radiobutton2,'Enable','on');
set(handles.pushbutton5,'Enable','on');
set(handles.popupmenu9,'Enable','on');

%Loading plot
set(handles.text9,'Enable','on');
set(handles.text10,'Enable','on');
set(handles.popupmenu10,'Enable','on');
set(handles.popupmenu11,'Enable','on');
set(handles.pushbutton12,'Enable','on');
set(handles.text17,'Enable','on');
set(handles.text18,'Enable','on');
set(handles.popupmenu19,'Enable','on');
set(handles.popupmenu20,'Enable','on');

%Residue
set(handles.pushbutton20,'Enable','on');
set(handles.pushbutton21,'Enable','on');

%Information panel:
text=sprintf('Model generated successully!');
handles.data.text=cprint(handles.text1,text,handles.data.text,0);
text=sprintf('Generated from PCA model matrices p and t:\np: Loading matrix\nt: Score matrix.');
handles.data.text=cprint(handles.text1,text,handles.data.text,1);
handles.data.text=cprint(handles.text1,'',handles.data.text,-1);
text=sprintf('Plot a Score plot, a Loading plot, a MEDA plot or Residual plot, by clicking on  the proper menu.');
handles.data.text=cprint(handles.text1,text,handles.data.text,0);

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function uipanel3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%Score Plot Submenu%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in popupmenu1.
%popupmenu1==PC X-axes
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
incoming_data_PC1=get(hObject,'Value');%Incoming data position
handles.data.PC1=incoming_data_PC1;

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu7.
%popupmenu7==PC Y-axes
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu7
incoming_data_PC2=get(hObject,'Value');%Incoming data position
handles.data.PC2=incoming_data_PC2;

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu16.
%popupmenu16==Label
function popupmenu16_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu16 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu16

%Information Panel:
text=sprintf('This is an optinal field that assign a name to each of the observations.');
handles.data.text=cprint(handles.text1,text,handles.data.text,0);
text=sprintf('Label is a cell having this format: {`x`,`x`,`y`,`y`,...,`z`,`z`}, containing as many tags as the number of observations.');
handles.data.text=cprint(handles.text1,text,handles.data.text,1);
text=sprintf('To use this option, define the tags cell array and charge it from the workspace by clicking on the REFRESH button.');
handles.data.text=cprint(handles.text1,text,handles.data.text,1);

incoming_data=get(hObject,'Value');%Incoming data position
string_evaluation=handles.data.new1{incoming_data};
handles.data.namePopupmenu16=string_evaluation;
if strcmp(string_evaluation,'emptylabel'),
    label={};
    handles.data.label={};
else
    label=evalin('base',string_evaluation);
    handles.data.label=label;
end

if ~isempty(handles.data.label),
    if max(size(label))~=size(handles.data.data_matrix,1) || min(size(label))~=1,
        errordlg('Label must have as many tags as number of observations in the data matrix.');
        handles.data.namePopupmenu16='emptylabel';
        handles.data.label={};
        nombres=cellstr(get(hObject,'String'));
        for i=1:length(nombres),
            if strcmp(nombres(i),'emptylabel'),
                val=i;
            end
        end
        set(hObject,'Value',val);
    end
end

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
new1=[];
handles.data.label={};
set(hObject,'String',{'emptylabel'});

handles.data.WorkSpace=evalin('base','who');
if ~isempty(handles.data.WorkSpace),
    contents=get(hObject,'String');
    new1=[];
    for i=1:length(handles.data.WorkSpace),
        new1=[new1 handles.data.WorkSpace(i,:)];
    end
    new1=[contents,new1];
    set(hObject,'String',strvcat(new1));
end

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

nombres=cellstr(get(hObject,'String'));
for i=1:length(nombres),
    if strcmp(nombres(i),'emptylabel'),
        val=i;
    end
end

set(hObject,'Value',val);
handles.data.new1=new1;
handles.data.namePopupmenu16='emptylabel';
guidata(hObject, handles);

% --- Executes on selection change in popupmenu15.
%popupmenu15==Classes
function popupmenu15_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu15 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu15

%Information panel
text=sprintf('This is an optinal field that colour the observations according to the value assigned to each of them.');
handles.data.text=cprint(handles.text1,text,handles.data.text,0);
text=sprintf('Classes is an array having this format: [1,1,2,2,...,3,3], containing as many entries as the number of observations.');
handles.data.text=cprint(handles.text1,text,handles.data.text,1);
text=sprintf('To use this option, define the array and charge it from the work space by clocking on the REFRESH button.');
handles.data.text=cprint(handles.text1,text,handles.data.text,1);

incoming_data=get(hObject,'Value');
string_evaluation=handles.data.new2{incoming_data};
handles.data.namePopupmenu15=string_evaluation;
if strcmp(string_evaluation,'emptyclasses'),
    classes=[];
    handles.data.classes=[];
else
    classes=evalin('base',string_evaluation);
    handles.data.classes=classes;
end

if ~isempty(handles.data.classes),
    if max(size(classes))~=size(handles.data.data_matrix,1) || min(size(classes))~=1,
        errordlg('Classes must have as many tags as number of observations in the data matrix.');
        handles.data.namePopupmenu15='emptyclasses';
        handles.data.classes=[];
        nombres=cellstr(get(hObject,'String'));
        for i=1:length(nombres),
            if strcmp(nombres(i),'emptyclasses'),
                val=i;
            end
        end
        set(hObject,'Value',val);
    end
end

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popupmenu15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
new2=[];
handles.data.classes=[];
set(hObject,'String',{'emptyclasses'});

handles.data.WorkSpace=evalin('base','who');
if ~isempty(handles.data.WorkSpace),
    contents=get(hObject,'String');
    new2=[];
    for i=1:length(handles.data.WorkSpace),
        new2=[new2 handles.data.WorkSpace(i,:)];
    end
    new2=[contents,new2];
    set(hObject,'String',strvcat(new2));
end

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

nombres=cellstr(get(hObject,'String'));
for i=1:length(nombres),
    if strcmp(nombres(i),'emptyclasses'),
        val=i;
    end
end
handles.data.namePopupmenu15='emptyclasses';
set(hObject,'Value',val);
handles.data.new2=new2;
guidata(hObject, handles);

% --- Executes on button press in pushbutton4.
%Pushbutton==Plot
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.data.PC1) || isempty(handles.data.PC2),
    errordlg('Error: select the combination of PCs to plot the scores.');
end

all_opened_graphs=get(0,'Children');
new_sp_ID_figures=[];
new_sp_matrix={};
clean_ind=[];

for i=1:length(handles.data.sp_ID_figures),
    if ~isempty(find(handles.data.sp_ID_figures(i)==all_opened_graphs,1)),
        new_sp_ID_figures=[new_sp_ID_figures handles.data.sp_ID_figures(i)];
        new_sp_matrix={new_sp_matrix{:} handles.data.sp_matrix{:,i}};
    else
        clean_ind=[clean_ind i];%Identificadores de lo Score Plots cerrados
        for j=1:length(clean_ind),
            aux=clean_ind(j);
            handles.data.clean_control(aux)=0;
            handles.data.CORTES{1,aux}=[];
            handles.data.weightDummy{1,aux}=[];
            handles.data.matrix_2PCs{1,aux}=[];
            %Dummy:
            M=size(handles.data.data_matrix,1);
            dummy=zeros(1,M);
            handles.data.dummy{1,aux}=dummy;
            handles.data.dummyRED=dummy;
            handles.data.dummyGREEN=dummy;
        end
    end
end

handles.data.sp_ID_figures=new_sp_ID_figures;%Vector actualizado con los identificadores de los Score Plots abiertos 
handles.data.sp_matrix=new_sp_matrix;

if isempty(handles.data.label) && isempty(handles.data.classes),
    [T,TT]=scores_pca(handles.data.data_matrix,[handles.data.PC1 handles.data.PC2],[],handles.data.prep,1);
else if ~isempty(handles.data.label) && isempty(handles.data.classes),
        [T,TT]=scores_pca(handles.data.data_matrix,[handles.data.PC1 handles.data.PC2],[],handles.data.prep,1,handles.data.label);
    else if isempty(handles.data.label) && ~isempty(handles.data.classes),
            [T,TT]=scores_pca(handles.data.data_matrix,[handles.data.PC1 handles.data.PC2],[],handles.data.prep,1,[],handles.data.classes);
        else [T,TT]=scores_pca(handles.data.data_matrix,[handles.data.PC1 handles.data.PC2],[],handles.data.prep,1,handles.data.label,handles.data.classes);
        end
    end
end
fig=gcf;
set(fig,'Tag','ScorePlot');%En la opci�n etiqueta se indica que el gr�fico es un Score Plot

matrixPCs_oMEDA=[T(:,handles.data.PC1),T(:,handles.data.PC2)];

handles.data.sp_ID_figures=[handles.data.sp_ID_figures fig];%Vector con los identificadores de los Score Plots abiertos
handles.data.sp_matrix={handles.data.sp_matrix{:} matrixPCs_oMEDA};

%oMEDA (Select)
set(handles.pushbutton7,'Enable','on');
%Information panel
text=sprintf('To perform an oMEDA plot push on the SELECT button in the oMEDA menu (upon selection of a Score Plot).');
handles.data.text=cprint(handles.text1,text,handles.data.text,0);

guidata(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%oMEDA Submenu%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton7.
%pushbutton7==SELECT
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ID_list=get(0,'Children');
ID=ID_list(2);%Identificador de la gr�fica seleccionada (debe ser un Score Plot).

check_tag=get(ID,'Tag');
if strcmp(check_tag,'ScorePlot'),
    figure(ID);%Lanzar el Score Plot seleccionado para hacer oMEDA.
else
    errordlg('To perform oMEDA you must select a Score Plot.');
    return;
end

text=sprintf('Over the selected Score Plot draw a polinomial enclosing the required points and push on the (+) button to assign them +1 (green points) or on the (-) button to assign them -1 (red points).');
handles.data.text=cprint(handles.text1,text,handles.data.text,0);

%Es necesario recuperar los datos del Score Plot seleccionado, es decir las observaciones ploteadas en el eje x e y:
%Voy a recorrer el vector de gcfs de score plots que se llama
%handles.data.sp_ID_figures, para buscar en que posici�n est� el gcf ID.
for i=1:length(handles.data.sp_ID_figures),
    if handles.data.sp_ID_figures(i)==ID,
        % codigo de compr de que est� vacio
%         ID
%         size(handles.data.sp_matrix)
        matrix_2PCs=handles.data.sp_matrix{:,i};
    end
end

irr_pol=impoly;
vertex=getPosition(irr_pol);
N=size(vertex,1);%Tama�o de la matriz:
%filas: n�mero de v�rtices del polinomio irregular
%columnas: contiene 2 columnas: coordenada x y coordenada y de cada v�rtice.

%PASO 1:
%Calcular los par�metros A, B y C de la ecuaci�n normal de la recta, para
%todas las rectas que formen el polinomio irregular dibujado por el usuario
A=[];
B=[];
C=[];
for i=1:N,%Desde 1 hasta el n�mero de v�rtices que tenga el polinomio
    %irregular, voy a hacer lo siguiente:
    
    %Coordenadas de un v�rtice:
    x1=vertex(i,1);
    y1=vertex(i,2);
    
    %Cooredenadas del siguiente v�rtice:
    %El if controla el caso en que ya se hayan cogido todos los v�rtices,
    %el v�rtce en ese caso ser� el primero de ellos, para cerrar la figura.
    if i==N,
        x2=vertex(1,1);
        y2=vertex(1,2);
    else
        x2=vertex(i+1,1);
        y2=vertex(i+1,2);
    end
    
    %Coordenadas del vector director de la recta que une ambos v�rtices:
    u1=x2-x1;
    u2=y2-y1;
    
    A=[A,u2];%Lista de u2(segunda coordenada del vector director)
    B=[B,-u1];%Lista de u1 (primera coordenada del vector director)
    c=(u1*y1)-(u2*x1);%C�lculo del par�metro C de la ec.normal de la recta.
    C=[C,c];%Lista del par�metro C, uno por recta.
end

%PASO 2:
%Obtener los puntos de corte entre cada una de las anteriores rectas y la
%semirrecta(paralela al eje X) que forma el punto del Score matrix a estudio.
M=size(handles.data.data_matrix,1);%Number of observations in the score matrix.
X=[];
corte=0;
CORTES=[];

for j=1:M, %Se recorren todas las observaciones
    Y=matrix_2PCs(j,2);
    corte=0;
    for k=1:N,%Todas las rectas del poligono irregular
        X=(-(B(k)*Y)-C(k))/A(k);
        
        if k+1>N,
            if (Y>min(vertex(k,2),vertex(1,2)))&&(Y<max(vertex(k,2),vertex(1,2))),
                if X>matrix_2PCs(j,1),
                    corte=corte+1;
                end
            end
        else
            if (Y>min(vertex(k,2),vertex(k+1,2)))&&(Y<max(vertex(k,2),vertex(k+1,2))),
                if X>matrix_2PCs(j,1),
                    corte=corte+1;
                end
            end
        end
    end
    CORTES=[CORTES,corte];
end

set(handles.pushbutton8,'Enable','on');
set(handles.pushbutton9,'Enable','on');
set(handles.pushbutton10,'Enable','on');

handles.data.CORTES{1,ID}=CORTES;
handles.data.matrix_2PCs{1,ID}=matrix_2PCs;
handles.data.clean_control(ID)=handles.data.clean_control(ID)+1;

guidata(hObject,handles);

% --- Executes on button press in pushbutton8.
%pushbutton8==- -> RED
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ID_list=get(0,'Children');
ID=ID_list(2);%gcf del score plot que me interesa

check_tag=get(ID,'Tag');
if strcmp(check_tag,'ScorePlot'),
    figure(ID);%Ya tengo el score plot pinchado(al que le quiero hacer oMEDA) en primera plana.
    hold on;
else
    errordlg('To perform oMEDA you must select a Score Plot.');
    return;
end

M=size(handles.data.data_matrix,1);
CortesVector=handles.data.CORTES{1,ID};
matrix_2PCs=handles.data.matrix_2PCs{1,ID};

for l=1:M,
    if mod(CortesVector(l),2)==1,
        
        Xdata=matrix_2PCs(l,1);
        Ydata=matrix_2PCs(l,2);
        
        coord=plot(Xdata,Ydata);
        set(coord,'marker','o');
        set(coord,'markersize',6);
        set(coord,'markerfacecolor','r');
        set(coord,'markeredgecolor','r');
        
        %Dummy:
        handles.data.dummyRED(l)=-1;

        handles.data.clean_control(ID)=handles.data.clean_control(ID)+1;
    end
end

handles.data.dummy{1,ID}=handles.data.dummyGREEN+handles.data.dummyRED;
set(handles.pushbutton6,'Enable','on');
set(handles.pushbutton11,'Enable','on');
%Information Panel:
text=sprintf('Push on the SELECT button again to repeat this process and select all the required points.');
handles.data.text=cprint(handles.text1,text,handles.data.text,0);
text=sprintf('Optionally push on the Trend button and draw a line over the Score Plot to include weigths in the analysis.\n');
handles.data.text=cprint(handles.text1,text,handles.data.text,1);
text=sprintf('Finally push on the Plot button to obtain the oMEDA plot.');
handles.data.text=cprint(handles.text1,text,handles.data.text,1);
guidata(hObject,handles);

% --- Executes on button press in pushbutton9.
%pushbutton9==+ -> GREEN
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ID_list=get(0,'Children');
ID=ID_list(2);%gcf del score plot que me interesa

check_tag=get(ID,'Tag');
if strcmp(check_tag,'ScorePlot'),
    figure(ID);%Ya tengo el score plot pinchado(al que le quiero hacer oMEDA) en primera plana.
    hold on;
else
    errordlg('To perform oMEDA you must select a Score Plot.');
    return;
end

M=size(handles.data.data_matrix,1);
CortesVector=handles.data.CORTES{1,ID};
matrix_2PCs=handles.data.matrix_2PCs{1,ID};

for l=1:M,
    if mod(CortesVector(l),2)==1,
        Xdata=matrix_2PCs(l,1);
        Ydata=matrix_2PCs(l,2);
        
        coord=plot(Xdata,Ydata);
        set(coord,'marker','o');
        set(coord,'markersize',6);
        set(coord,'markerfacecolor','g');
        set(coord,'markeredgecolor','g');
        
        
        handles.data.dummyGREEN(l)=1;
        handles.data.clean_control(ID)=handles.data.clean_control(ID)+1;
    end
end

handles.data.dummy{1,ID}=handles.data.dummyGREEN+handles.data.dummyRED;
set(handles.pushbutton6,'Enable','on');
set(handles.pushbutton11,'Enable','on');
%Information panel:
text=sprintf('Push on the SELECT button again to repeat this process and select all the required points.');
handles.data.text=cprint(handles.text1,text,handles.data.text,0);
text=sprintf('Optionally push on the Trend button and draw a line over the Score Plot to include weigths in the analysis.\n');
handles.data.text=cprint(handles.text1,text,handles.data.text,1);
text=sprintf('Finally push on the Plot button to obtain the oMEDA plot.');
handles.data.text=cprint(handles.text1,text,handles.data.text,1);
guidata(hObject,handles);


% --- Executes on button press in pushbutton11.
%pushbutton11==TREND LINE
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ID_list=get(0,'Children');
ID=ID_list(2);%gcf del score plot que me interesa

check_tag=get(ID,'Tag');
if strcmp(check_tag,'ScorePlot'),
    figure(ID);%Ya tengo el score plot pinchado(al que le quiero hacer oMEDA) en primera plana.
    hold on;
else
    errordlg('To perform oMEDA you must select a Score Plot.');
    return;
end

matrix_2PCs=handles.data.matrix_2PCs{1,ID};
trend_line=imline;
setColor(trend_line,[0 0 0]);
vertex_line=getPosition(trend_line);

x1=vertex_line(1,1);
y1=vertex_line(1,2);
x2=vertex_line(2,1);
y2=vertex_line(2,2);

%Coordenadas del vector director de la recta que une ambos v�rtices:
u1=x2-x1;
u2=y2-y1;

%La ecuaci�n de la recta tendencia es:
A=u2;
B=-u1;
C=(u1*y1)-(u2*x1);

%Quiero el punto de corte de la tendencia con la recta que va de la observaci�n
%a la l�nea tendencia en perpendicular. Esto para cada una de las
%observaciones.
Cutoff_points=[];
M=size(handles.data.data_matrix,1);
for m=1:M,
    p1=matrix_2PCs(m,1);
    p2=matrix_2PCs(m,2);
    
    %El vector director de la recta que va de la observacion a la
    %tendencia en perpendicular es:
    v1=A;
    v2=B;
    
    %La ecuacuaci�n de la recta es:
    A2=v2;
    B2=-v1;
    C2=(v1*p2)-(v2*p1);
    
    %Ahora se obtiene el punto de corte de ambas rectas:
    %Ax+By+C=0;
    %A2x+B2y+C2=0;
    y_corte=(-C2+(A2/A)*C)/(((-A2/A)*B)+B2);
    x_corte=((-B*y_corte)-C)/A;
    Cutoff_points(m,1)=x_corte;
    Cutoff_points(m,2)=y_corte;
end
%Quedarse con la menor distancia entre un 1 y -1
lowest_dist=Inf;
ind1=1;
ind2=1;
dummy=handles.data.dummy{1,ID};
for k=1:M,
    if dummy(k)==1,
        %Coordenadas del punto que tiene asignado un 1 en la variable
        %dummy
        p1=Cutoff_points(k,1);
        p2=Cutoff_points(k,2);
        
        for l=1:M,
            if dummy(l)==-1,
                q1=Cutoff_points(l,1);
                q2=Cutoff_points(l,2);
                dist=sqrt((q1-p1)^2+(q2-p2)^2);
                if dist< lowest_dist,
                    lowest_dist=dist;
                    ind1=k;
                    ind2=l;
                end
            end
        end
        
    end
end

%Construcci�n de la nueva DUMMY con pesos:
%Calcular el punto medio entre las observaciones m�s cercanas obtenidas
%enteriormente, este ser� el nuevo cero para asignar pesos.
c1=Cutoff_points(ind1,:);
c2=Cutoff_points(ind2,:);
NewCenter=(c1+c2)/2;

%Asignaci�n de pesos
for m=1:M,
    weights(m)=sum((Cutoff_points(m,:)-NewCenter).^2);
end
weightDummy=weights.*dummy;

handles.data.weightDummy{1,ID}= weightDummy;
handles.data.clean_control(ID)=handles.data.clean_control(ID)+1;
%Information panel:
text=sprintf('Finally push on the Plot button to obtain the oMEDA plot.');
handles.data.text=cprint(handles.text1,text,handles.data.text,0);
text=sprintf('If needed, push on the SELECT button again to repeat this process and select all the required points.');
handles.data.text=cprint(handles.text1,text,handles.data.text,1);

guidata(hObject,handles);

% --- Executes on button press in pushbutton10.
%pushbutton10==CLEAN
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ID_list=get(0,'Children');
ID=ID_list(2);%gcf del score plot que me interesa

check_tag=get(ID,'Tag');
if strcmp(check_tag,'ScorePlot'),
    figure(ID);%Ya tengo el score plot pinchado(al que le quiero hacer oMEDA) en primera plana.
else
    errordlg('To clean a figure this must be a Score Plot.');
    return;
end

clean_times=handles.data.clean_control(ID);

a=get(ID,'Children');
b=a(length(a));
c=get(b,'Children');
delete(c(1:clean_times));

handles.data.clean_control(ID)=0;
handles.data.weightDummy{1,ID}=[];
%Dummy:
M=size(handles.data.data_matrix,1);
dummy=zeros(1,M);
handles.data.dummy{1,ID}=dummy;
handles.data.dummyRED=dummy;
handles.data.dummyGREEN=dummy;

guidata(hObject,handles);


% --- Executes on button press in pushbutton6.
%pushbutton6==Plot (oMEDA)
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ID_list=get(0,'Children');
ID=ID_list(2);%gcf del score plot que me interesa

check_tag=get(ID,'Tag');
if strcmp(check_tag,'ScorePlot'),

else
    errordlg('To perform oMEDA you must select a Score Plot.');
    return;
end

if isequal(handles.data.dummy{1,ID},zeros(1,size(handles.data.dummy{1,ID},2))) && isempty(handles.data.weightDummy{1,ID}),
    errordlg('To perform oMEDA you must select a Score Plot with at least one object selected.');
    return;
end

if ~isempty(handles.data.weightDummy{1,ID}),
    handles.data.weightDummy{1,ID}=handles.data.weightDummy{1,ID}./abs(max(handles.data.weightDummy{1,ID}));
    omeda_pca(handles.data.data_matrix,[min(handles.data.PC1,handles.data.PC2) max(handles.data.PC1,handles.data.PC2)],handles.data.data_matrix,handles.data.weightDummy{1,ID}',handles.data.prep,1,handles.data.label_LP);
else
    omeda_pca(handles.data.data_matrix,[min(handles.data.PC1,handles.data.PC2) max(handles.data.PC1,handles.data.PC2)],handles.data.data_matrix,handles.data.dummy{1,ID}',handles.data.prep,1,handles.data.label_LP);
end

guidata(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Residue%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton20.
%pushbutton20==Plot(residue)
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
E=sqresiduals_pca(handles.data.data_matrix,min(handles.data.PCs):max(handles.data.PCs),[],handles.data.prep,1,handles.data.label);
text=sprintf('Plotted squared residuals in the observations.');
handles.data.text=cprint(handles.text1,text,handles.data.text,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Loading Plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in popupmenu10.
%popupmenu10==PCs X-axes
function popupmenu10_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu10 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu10
incoming_data_PC1_LP=get(hObject,'Value');%Incoming data position
handles.data.PC1_LP=incoming_data_PC1_LP;

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu11.
%popupmenu11==PCs Y-axes
function popupmenu11_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu11 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu11
incoming_data_PC2_LP=get(hObject,'Value');%Incoming data position
handles.data.PC2_LP=incoming_data_PC2_LP;

guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popupmenu11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu20.
%popupmenu20==Label
function popupmenu20_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu20 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu20
text=sprintf('This is an optinal field that assign a name to each of the variables');
handles.data.text=cprint(handles.text1,text,handles.data.text,0);
text=sprintf('Label is a cell having this format: {`x`,`x`,`y`,`y`,...,`z`,`z`}, containing as many tags as the number of variables.');
handles.data.text=cprint(handles.text1,text,handles.data.text,1);
text=sprintf('To use this option, define the tags cell array and chare it from the work space by clicking on the REFRESH button.');
handles.data.text=cprint(handles.text1,text,handles.data.text,1);

incoming_data=get(hObject,'Value');
string_evaluation=handles.data.new3{incoming_data};
handles.data.namePopupmenu20=string_evaluation;
if strcmp(string_evaluation,'emptylabel'),
    label_LP={};
    handles.data.label_LP={};
else
    label_LP=evalin('base',string_evaluation);
    handles.data.label_LP=label_LP;
end

if ~isempty(handles.data.label_LP),
    if max(size(label_LP))~=size(handles.data.data_matrix,2) || min(size(label_LP))~=1,
        errordlg('Label must have as many tags as number of variables in the data matrix.');
        handles.data.namePopupmenu20='emptylabel';
        handles.data.label_LP={};
        nombres=cellstr(get(hObject,'String'));
        for i=1:length(nombres),
            if strcmp(nombres(i),'emptylabel'),
                val=i;
            end
        end
        set(hObject,'Value',val);
    end
end

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
new3=[];
handles.data.label_LP={};
set(hObject,'String',{'emptylabel'});

handles.data.WorkSpace=evalin('base','who');%nombres de las variables
if ~isempty(handles.data.WorkSpace),
    contents=get(hObject,'String');
    new3=[];
    for i=1:length(handles.data.WorkSpace),
        new3=[new3 handles.data.WorkSpace(i,:)];
    end
    new3=[contents,new3];
    set(hObject,'String',strvcat(new3));
end

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

nombres=cellstr(get(hObject,'String'));
for i=1:length(nombres),
    if strcmp(nombres(i),'emptylabel'),
        val=i;
    end
end
handles.data.namePopupmenu20='emptylabel';
set(hObject,'Value',val);
handles.data.new3=new3;
guidata(hObject, handles);

% --- Executes on selection change in popupmenu19.
%popupmenu20==Classes
function popupmenu19_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu19 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu19
text=sprintf('This is an optinal field that colour the variables according to the value assigned to each of them.');
handles.data.text=cprint(handles.text1,text,handles.data.text,0);
text=sprintf('Classes is an array having this format: [1,1,2,2,...,3,3], containing as many entries as the number of variables.');
handles.data.text=cprint(handles.text1,text,handles.data.text,1);
text=sprintf('To use this option, define the array and chare it from the work space by clicking on the REFRESH button.');
handles.data.text=cprint(handles.text1,text,handles.data.text,1);

incoming_data=get(hObject,'Value');%Incoming data position
string_evaluation=handles.data.new4{incoming_data};%Nombre correspondiente a la posici�n
handles.data.namePopupmenu19=string_evaluation;
if strcmp(string_evaluation,'emptyclasses'),
    classes_LP={};
    handles.data.classes_LP={};
else
    classes_LP=evalin('base',string_evaluation);%Contenido de ese nombre(los datos en si)
    handles.data.classes_LP=classes_LP;
end

if ~isempty(handles.data.classes_LP),
    if max(size(classes_LP))~=size(handles.data.data_matrix,2) || min(size(classes_LP))~=1,
        errordlg('Classes must have as many entries as number of variables in the data matrix.');
        handles.data.namePopupmenu19='emptyclasses';
        handles.data.classes_LP=[];
        nombres=cellstr(get(hObject,'String'));
        for i=1:length(nombres),
            if strcmp(nombres(i),'emptyclasses'),
                val=i;
            end
        end
        set(hObject,'Value',val);
    end
end

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
new4=[];
handles.data.classes_LP=[];
set(hObject,'String',{'emptyclasses'});

handles.data.WorkSpace=evalin('base','who');%nombres de las variables
if ~isempty(handles.data.WorkSpace),
    contents=get(hObject,'String');
    new4=[];
    for i=1:length(handles.data.WorkSpace),
        new4=[new4 handles.data.WorkSpace(i,:)];
    end
    new4=[contents,new4];
    set(hObject,'String',strvcat(new4));
end

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

nombres=cellstr(get(hObject,'String'));
for i=1:length(nombres),
    if strcmp(nombres(i),'emptyclasses'),
        val=i;
    end
end
handles.data.namePopupmenu19='emptyclasses';
set(hObject,'Value',val);
handles.data.new4=new4;
guidata(hObject, handles);


% --- Executes on button press in pushbutton12.
%pushbutton12==Plot (loading plot)
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isempty(handles.data.PC1_LP) || isempty(handles.data.PC2_LP),
    errordlg('Error: select the combination of PCs to plot the loadings.');
end

all_opened_graphs=get(0,'Children');
new_lp_ID_figures=[];
new_lp_matrix={};

for i=1:length(handles.data.lp_ID_figures),
    if ~isempty(find(handles.data.lp_ID_figures(i)==all_opened_graphs,1)),
        new_lp_ID_figures=[new_lp_ID_figures handles.data.lp_ID_figures(i)];
        new_lp_matrix={new_lp_matrix{:} handles.data.lp_matrix{:,i}};
    end
end

handles.data.lp_ID_figures=new_lp_ID_figures;%Identificadores de los Loadings Plots abiertos actualizado
handles.data.lp_matrix=new_lp_matrix;

if isempty(handles.data.label_LP) && isempty(handles.data.classes_LP),
    P = loadings_pca(handles.data.data_matrix,[handles.data.PC1_LP handles.data.PC2_LP],handles.data.prep,2);
else if ~isempty(handles.data.label_LP) && isempty(handles.data.classes_LP),
        P = loadings_pca(handles.data.data_matrix,[handles.data.PC1_LP handles.data.PC2_LP],handles.data.prep,2,handles.data.label_LP);
    else if isempty(handles.data.label_LP) && ~isempty(handles.data.classes_LP),
            P = loadings_pca(handles.data.data_matrix,[handles.data.PC1_LP handles.data.PC2_LP],handles.data.prep,2,' ',handles.data.classes_LP);
        else P = loadings_pca(handles.data.data_matrix,[handles.data.PC1_LP handles.data.PC2_LP],handles.data.prep,2,handles.data.label_LP,handles.data.classes_LP);
        end
    end
end
fig=gcf;
set(fig,'Tag','LoadingPlot');%A cada LoadingPlot que abro le pongo en su propiedad 'Tag' que es un LoadingPlot

matrixPCs_MEDA_LP=[P(:,handles.data.PC1_LP),P(:,handles.data.PC2_LP)];

handles.data.lp_ID_figures=[handles.data.lp_ID_figures fig];%Identificadores de los Score Plots abiertos
handles.data.lp_matrix={handles.data.lp_matrix{:} matrixPCs_MEDA_LP};

set(handles.pushbutton22,'Enable','on');

text=sprintf('To perform a MEDA plot, push on the SELECT button in the MEDA menu (upon selection of Loading Plot).');
handles.data.text=cprint(handles.text1,text,handles.data.text,0);

guidata(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%MEDA Submenu%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%edit6==thresold
function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double
thres=str2double(get(hObject,'String'));
handles.data.thres = thres;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', 0.1);
edit6_Callback(hObject, eventdata, handles);

% --- Executes on button press in radiobutton1.
%radiobutton1==thresold
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of radiobutton1
%Si radio button se�alado q ecit 6 este ON si no se�alado q este OFF
if get(handles.radiobutton1, 'Value'),
    set(handles.edit6, 'Enable', 'on');
    set(handles.text5, 'Enable', 'on');
else
    set(handles.edit6, 'Enable', 'off');
    set(handles.text5, 'Enable', 'off');
end

guidata(hObject,handles);


% --- Executes on button press in radiobutton2.
%radiobutton1==seriation
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2


% --- Executes on selection change in popupmenu9.
%popupmenu9==MEDA popupmenu
function popupmenu9_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu9 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu9

PCs_MEDA_position=get(hObject,'Value');%Incoming data position
contents=get(hObject,'String');
PCs_MEDA=contents(PCs_MEDA_position,:);%Nombre correspondiente a la posici�n

handles.data.PCs_MEDA=PCs_MEDA;

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton5.
%pushbutton5==Plot (MEDA)
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.data.PCs_MEDA),
    errordlg('To perform MEDA select the PCs from the popupmenu.');
    return;
end

if get(handles.radiobutton1,'Value')==1 && get(handles.radiobutton2,'Value')==0,
    handles.data.opt=2;
else if get(handles.radiobutton1,'Value')==0 && get(handles.radiobutton2,'Value')==1,
        handles.data.opt=3;
    else if get(handles.radiobutton2,'Value')==0 && get(handles.radiobutton2,'Value')==0,
            handles.data.opt=1;
        else handles.data.opt=4;
        end
    end
end

[meda_map,meda_dis]=meda_pca(handles.data.data_matrix,min(str2num(handles.data.PCs_MEDA)):max(str2num(handles.data.PCs_MEDA)),handles.data.prep,handles.data.thres,handles.data.opt,handles.data.label_LP);

guidata(hObject,handles);

% --- Executes on button press in pushbutton22.
%pushbutton22==Select (MEDA)
function pushbutton22_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ID_list=get(0,'Children');
ID=ID_list(2);%gcf del score plot que me interesa

check_tag=get(ID,'Tag');
if strcmp(check_tag,'LoadingPlot'),
    figure(ID);%Ya tengo el score plot pinchado(al que le quiero hacer oMEDA) en primera plana.
else
    errordlg('To perform MEDA over a Loading Plot you must select one loading plot.');
    return;
end

%Ahora vamos a recuperar su matriz:
%Voy a recorrer el vector de gcfs de score plots
%handles.data.sp_ID_figures, para buscar en que posici�n esta el gcf ID.
for i=1:length(handles.data.lp_ID_figures),
    if handles.data.lp_ID_figures(i)==ID,
        matrix_2PCs=handles.data.lp_matrix{:,i};
    end
end

irr_pol=impoly;
vertex=getPosition(irr_pol);
N=size(vertex,1);%Matrix size:
%rows: number of vertex in the irregular polynomial
%columns: this matrix contains always 2 colomns: one for the X
%coordinate and the other one for the Y coordinate of each
%vertex.

%PASO 1:
%Calcular los par�metros A, B y C de la ecuaci�n normal de la recta, para
%todas las rectas que formen el polinomio irregular dibujado por el usuario
A=[];
B=[];
C=[];
for i=1:N,%Desde 1 hasta el n�mero de v�rtices que tenga el polinomio
    %irregular, voy a hacer lo siguiente:
    
    %Coordenadas de un v�rtice
    x1=vertex(i,1);
    y1=vertex(i,2);
    
    %Cooredenadas del siguiente v�rtice:
    %El if controla el caso en que ya se hayan cogido todos los v�rtices,
    %el v�rtce en ese caso ser� el primero de ellos, para cerrar la figura.
    if i==N,
        x2=vertex(1,1);
        y2=vertex(1,2);
    else
        x2=vertex(i+1,1);
        y2=vertex(i+1,2);
    end
    
    %Coordenadas del vector director de la recta que une ambos v�rtices:
    u1=x2-x1;
    u2=y2-y1;
    
    A=[A,u2];%Lista de u2(segunda coordenada del vector director)
    B=[B,-u1];%Lista de u1 (primera coordenada del vector director)
    c=(u1*y1)-(u2*x1);%C�lculo del par�metro C de la ec.normal de la recta.
    C=[C,c];%Lista del par�metro C, uno por recta.
end

%PASO 2:
%Obtener los puntos de corte entre cada una de las anteriores rectas y la
%semirrecta(paralela al eje X) que forma el punto del Score matrix a estudio.
M=size(handles.data.data_matrix,2);%Number of observations in the score matrix.
X=[];
corte=0;
CORTES=[];

for j=1:M, %All the observations from the Score Matrix: t
    Y=matrix_2PCs(j,2);
    corte=0;
    for k=1:N,%Todas las rectas del poligono irregular
        X=(-(B(k)*Y)-C(k))/A(k);
        
        if k+1>N,
            if (Y>min(vertex(k,2),vertex(1,2)))&&(Y<max(vertex(k,2),vertex(1,2))),
                if X>matrix_2PCs(j,1),
                    corte=corte+1;
                end
            end
        else
            if (Y>min(vertex(k,2),vertex(k+1,2)))&&(Y<max(vertex(k,2),vertex(k+1,2))),
                if X>matrix_2PCs(j,1),
                    corte=corte+1;
                end
            end
        end
    end
    CORTES=[CORTES,corte];
end

CortesVector=CORTES;
vector_vars=[];
for l=1:M,
    if mod(CortesVector(l),2)==1,
        Xdata=matrix_2PCs(l,1);
        Ydata=matrix_2PCs(l,2);
        
        coord=plot(Xdata,Ydata);
        set(coord,'marker','o');
        set(coord,'markersize',6);
        set(coord,'markerfacecolor','b');
        set(coord,'markeredgecolor','b');
        
        vector_vars=[vector_vars l];
    end
end

if get(handles.radiobutton1,'Value')==1 && get(handles.radiobutton2,'Value')==0,
    handles.data.opt=2;
else if get(handles.radiobutton1,'Value')==0 && get(handles.radiobutton2,'Value')==1,
        handles.data.opt=3;
    else if get(handles.radiobutton2,'Value')==0 && get(handles.radiobutton2,'Value')==0,
            handles.data.opt=1;
        else handles.data.opt=4;
        end
    end
end

[meda_map,meda_dis]=meda_pca(handles.data.data_matrix,[min(handles.data.PC1_LP,handles.data.PC2_LP) max(handles.data.PC1_LP,handles.data.PC2_LP)],handles.data.prep,handles.data.thres,handles.data.opt,handles.data.label_LP,vector_vars);

guidata(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Residue%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton21.
%pushbutton21==Plot (residue)
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
E=sqresiduals_pca(handles.data.data_matrix,min(handles.data.PCs):max(handles.data.PCs),[],handles.data.prep,2,handles.data.label_LP);
text=sprintf('Plotted squared residuals in the variables.');
handles.data.text=cprint(handles.text1,text,handles.data.text,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%Information panel%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%text1==Static text-Information panel
% --- Executes during object creation, after setting all properties.
function text1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
