function varargout = EDA(varargin)
%%
%Initialize EDA Graphic User Interface
%
% coded by: Elena Jiménez Mañas (elenajm@correo.ugr.es)
%           Rafael Rodriguez Gomez (rodgom@ugr.es)
% version: 2.0
% last modification: 31/Jan/15.
%
% Copyright (C) 2014  Elena Jiménez Mañas
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
% EDA M-file for EDA.fig
%      EDA, by itself, creates a new EDA or raises the existing
%      singleton*.
%
%      H = EDA returns the handle to a new EDA or the handle to
%      the existing singleton*.
%
%      EDA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EDA.M with the given input arguments.
%
%      EDA('Property','Value',...) creates a new EDA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EDA_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EDA_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EDA

% Last Modified by GUIDE v2.5 25-Jan-2015 13:26:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @EDA_OpeningFcn, ...
    'gui_OutputFcn',  @EDA_OutputFcn, ...
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

% --- Executes just before EDA is made visible.
function EDA_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EDA (see VARARGIN)

% Choose default command line output for EDA
handles.output = hObject;

%Change icon
warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
javaFrame = get(hObject,'JavaFrame');
javaFrame.setFigureIcon(javax.swing.ImageIcon('icon.jpg'));

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EDA wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%Load the logo image
axes(handles.axesEDA);
% kk1 = get(gca, 'position');
% disp(kk1)
% kk2 = get(gcf, 'position');
% disp(kk2)
% set(gcf, 'position', [kk2(1) kk2(2) kk2(3)*1.3 kk2(4)*1.3]);
% set(gca, 'position', [kk1(1) kk1(2) kk1(3)*1.2 kk1(4)*1.2]);
logo=image(imread('logo2.png'));
axis image
set(get(logo,'Parent'),'YTick',[]);
set(get(logo,'Parent'),'XTick',[]);

% --- Outputs from this function are returned to the command line.
function varargout = EDA_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

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
