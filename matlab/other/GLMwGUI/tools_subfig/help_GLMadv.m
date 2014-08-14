function varargout = help_GLMadv(varargin)
% HELP_GLMADV MATLAB code for help_GLMadv.fig
%      HELP_GLMADV, by itself, creates a new HELP_GLMADV or raises the existing
%      singleton*.
%
%      H = HELP_GLMADV returns the handle to a new HELP_GLMADV or the handle to
%      the existing singleton*.
%
%      HELP_GLMADV('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HELP_GLMADV.M with the given input arguments.
%
%      HELP_GLMADV('Property','Value',...) creates a new HELP_GLMADV or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before help_GLMadv_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to help_GLMadv_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help help_GLMadv

% Last Modified by GUIDE v2.5 19-Sep-2013 11:47:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @help_GLMadv_OpeningFcn, ...
                   'gui_OutputFcn',  @help_GLMadv_OutputFcn, ...
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


% --- Executes just before help_GLMadv is made visible.
function help_GLMadv_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to help_GLMadv (see VARARGIN)

% Choose default command line output for help_GLMadv
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes help_GLMadv wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = help_GLMadv_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
global TAS PPMS Basispars

val = get(hObject,'Value');
switch val
    case 1
        ncols = 0;
    case 2
        ncols = 2;
    case 3
        ncols = 4;
    case 4
        ncols = 6;
    case 5
        ncols = 8;
    case 6
        ncols = 10;
    case 7
        ncols = 12;
end

Basispars(1) = ncols;
ihbasprs.ncols = ncols;  % Number of basis vectors for post-spike kernel
% SangWook's change 8/13
DTsim = 0.1;
ppms = PPMS;
timeafterspk = TAS;
        
ihbasprs.hpeaks = [DTsim*ppms*4 DTsim*ppms*timeafterspk*0.9];  % Peak location for first and last vectors

ihbasprs.b = .4;  % How nonlinear to make spacings
ihbasprs.absref = DTsim*ppms; % absolute refractory period 
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,DTsim,timeafterspk,ppms);
plot(linspace(0,TAS,length(ihbasis)),ihbasis)

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


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2
global TAS PPMS Basispars

val = get(hObject,'Value');
switch val
    case 1
        maxloc = 0;
    case 2
        maxloc = .7;
    case 3
        maxloc = .75;
    case 4
        maxloc = .8;
    case 5
        maxloc = .85;
    case 6
        maxloc = .9;
    case 7
        maxloc = .95;
end

Basispars(2) = maxloc;
ihbasprs.ncols = Basispars(1);  % Number of basis vectors for post-spike kernel
% SangWook's change 8/13
DTsim = 0.1;
ppms = PPMS;
timeafterspk = TAS;
        
ihbasprs.hpeaks = [DTsim*ppms*Basispars(3) DTsim*ppms*timeafterspk*maxloc];  % Peak location for first and last vectors

ihbasprs.b = .4;  % How nonlinear to make spacings
ihbasprs.absref = DTsim*ppms; % absolute refractory period 
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,DTsim,timeafterspk,ppms);
plot(linspace(0,TAS,length(ihbasis)),ihbasis)

% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3
global TAS PPMS Basispars

val = get(hObject,'Value');
switch val
    case 1
        minloc = 0;
    case 2
        minloc = 1;
    case 3
        minloc = 2;
    case 4
        minloc = 4;
    case 5
        minloc = 8;
    case 6
        minloc = 16;
end

Basispars(3) = minloc;
maxloc = Basispars(2);
ihbasprs.ncols = Basispars(1);  % Number of basis vectors for post-spike kernel
% SangWook's change 8/13
DTsim = 0.1;
ppms = PPMS;
timeafterspk = TAS;
        
ihbasprs.hpeaks = [DTsim*ppms*minloc DTsim*ppms*timeafterspk*maxloc];  % Peak location for first and last vectors

ihbasprs.b = .4;  % How nonlinear to make spacings
ihbasprs.absref = DTsim*ppms; % absolute refractory period 
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,DTsim,timeafterspk,ppms);
plot(linspace(0,TAS,length(ihbasis)),ihbasis)

% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Basispars

A = num2str(Basispars(1));
B = num2str(100*Basispars(2));
C = num2str(Basispars(3));

disp(['Number of basis functions : ' A])
disp(['Maximum peak location of basis function : ' B '%'])
disp(['Minimim location of basis functions : ' C])

close(handles.figure1)

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Basispars

Basispars = [8 .9 4];
close(handles.figure1)

% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1
global TAS PPMS Basispars

ihbasprs.ncols = Basispars(1);  % Number of basis vectors for post-spike kernel
% SangWook's change 8/13
DTsim = 0.1;
ppms = PPMS;
timeafterspk = TAS;
        
ihbasprs.hpeaks = [DTsim*ppms*Basispars(3) DTsim*ppms*timeafterspk*Basispars(2)];  % Peak location for first and last vectors

ihbasprs.b = .4;  % How nonlinear to make spacings
ihbasprs.absref = DTsim*ppms; % absolute refractory period 
[iht,ihbas,ihbasis] = makeBasis_PostSpike(ihbasprs,DTsim,timeafterspk,ppms);
plot(linspace(0,TAS,length(ihbasis)),ihbasis)
