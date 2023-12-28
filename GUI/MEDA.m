function varargout = MEDA(varargin)
%%
%Initialize MEDA Graphic User Interface
%
% coded by: Elena JimÃ©nez MaÃ±as (elenajm@correo.ugr.es)
%           Rafael Rodriguez Gomez (rodgom@ugr.es)
% last modification: 31/Jan/15.
%
% Copyright (C) 2016  University of Granada, Granada
% Copyright (C) 2016 Elena Jiménez Mañas, Rafael Rodriguez Gomez
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

%%
% MEDA M-file for MEDA.fig
%      MEDA, by itself, creates a new MEDA or raises the existing
%      singleton*.
%
%      H = MEDA returns the handle to a new MEDA or the handle to
%      the existing singleton*.
%
%      MEDA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MEDA.M with the given input arguments.
%
%      MEDA('Property','Value',...) creates a new MEDA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MEDA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MEDA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MEDA

% Last Modified by GUIDE v2.5 29-Jun-2018 22:28:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @MEDA_OpeningFcn, ...
    'gui_OutputFcn',  @MEDA_OutputFcn, ...
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

% --- Executes just before MEDA is made visible.
function MEDA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MEDA (see VARARGIN)

% Choose default command line output for MEDA
handles.output = hObject;
handles.data.version='1.3';
%Change icon
%warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
%javaFrame = get(hObject,'JavaFrame');
%javaFrame.setFigureIcon(javax.swing.ImageIcon('icon2.png'));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MEDA wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%Load the logo image
axes(handles.axesMEDA);
% kk1 = get(gca, 'position');
% disp(kk1)
% kk2 = get(gcf, 'position');
% disp(kk2)
% set(gcf, 'position', [kk2(1) kk2(2) kk2(3)*1.3 kk2(4)*1.3]);
% set(gca, 'position', [kk1(1) kk1(2) kk1(3)*1.2 kk1(4)*1.2]);
logo=image(imread('Logo.png'));
axis image
set(get(logo,'Parent'),'YTick',[]);
set(get(logo,'Parent'),'XTick',[]);

% --- Outputs from this function are returned to the command line.
function varargout = MEDA_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


% --- Executes on button press in pcaButton.
function pcaButton_Callback(hObject, eventdata, handles)
% hObject    handle to pcaButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PCA;

% --- Executes on button press in plsButton.
function plsButton_Callback(hObject, eventdata, handles)
% hObject    handle to plsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
PLS;


% --------------------------------------------------------------------
function help_Callback(hObject, eventdata, handles)
% hObject    handle to help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function about_menu_Callback(hObject, eventdata, handles)
% hObject    handle to about_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
myicon = imread('icon.png');
h = msgbox({strcat('Version number: ',handles.data.version),'','Main Author: José Camacho (jcamacho@ugr.es)'},'About','custom',myicon);
