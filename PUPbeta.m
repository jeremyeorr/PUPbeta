function varargout = PUPbeta(varargin)
% PUPBETA MATLAB code for PUPbeta.fig
%      PUPBETA, by itself, creates a new PUPBETA or raises the existing
%      singleton*.
%
%      H = PUPBETA returns the handle to a new PUPBETA or the handle to
%      the existing singleton*.
%
%      PUPBETA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PUPBETA.M with the given input arguments.
%
%      PUPBETA('Property','Value',...) creates a new PUPBETA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PUPbeta_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PUPbeta_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PUPbeta

% Last Modified by GUIDE v2.5 15-Apr-2015 16:48:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @PUPbeta_OpeningFcn, ...
    'gui_OutputFcn',  @PUPbeta_OutputFcn, ...
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


% --- Executes just before PUPbeta is made visible.
function PUPbeta_OpeningFcn(hObject, eventdata, handles, varargin)

global ExportDataSpreadsheet AnalyzeDataSpreadsheet handletext;
AnalyzeDataSpreadsheet = get(handles.edit2,'String');
ExportDataSpreadsheet = get(handles.edit1,'String');
handletext = handles.text5;
set(handletext,'String','Ready');
[currentdir,~,~] = fileparts(mfilename('fullpath')); cd(currentdir); %make the filepath of this file the current directory
drawnow;

% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PUPbeta (see VARARGIN)

% Choose default command line output for PUPbeta
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PUPbeta wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PUPbeta_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ExportDataSpreadsheet
ExportDataSpreadsheet = get(hObject,'String');
drawnow;
disp(['Export data using: ' ExportDataSpreadsheet]);
% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global AnalyzeDataSpreadsheet
AnalyzeDataSpreadsheet = get(hObject,'String');
disp(['Analyze data using: ' AnalyzeDataSpreadsheet]);
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%channel check
global ExportDataSpreadsheet handletext
%handletext = handles.text5;
displaytext = 'Getting channel data';
set(handletext,'String',displaytext);
try
    getEDFchannelnames();
    displaytext = ['Channel data copied to ' ExportDataSpreadsheet];
catch me
    displaytext = me.message;
end
disp(displaytext); set(handletext,'String',displaytext);


% --- Executes on button press in pushbutton2 "Run Export"
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ExportDataSpreadsheet handletext
disp(['Converting data using: ' ExportDataSpreadsheet]);
try
    ConvertToMat();
catch me
    displaytext=me.message;
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global AnalyzeDataSpreadsheet handletext
AnalyzeDataSpreadsheet=get(handles.edit2,'String');
try
    Analysis();
catch me
    displaytext=me.message;
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global AnalyzeDataSpreadsheet handletext
AnalyzeDataSpreadsheet=get(handles.edit2,'String');
try
    RatingsPUP();
catch me
    displaytext=me.message;
    disp(displaytext); set(handletext,'String',displaytext); drawnow;
end
