function varargout = help_ISIs(varargin)
% HELP_ISIS MATLAB code for help_ISIs.fig
%      HELP_ISIS, by itself, creates a new HELP_ISIS or raises the existing
%      singleton*.
%
%      H = HELP_ISIS returns the handle to a new HELP_ISIS or the handle to
%      the existing singleton*.
%
%      HELP_ISIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HELP_ISIS.M with the given input arguments.
%
%      HELP_ISIS('Property','Value',...) creates a new HELP_ISIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before help_ISIs_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to help_ISIs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help help_ISIs

% Last Modified by GUIDE v2.5 16-Jan-2014 21:46:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @help_ISIs_OpeningFcn, ...
                   'gui_OutputFcn',  @help_ISIs_OutputFcn, ...
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


% --- Executes just before help_ISIs is made visible.
function help_ISIs_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to help_ISIs (see VARARGIN)

% Choose default command line output for help_ISIs
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes help_ISIs wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = help_ISIs_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function axes5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes5
load test_spike_time_data
ppms = 10;
ISIcutoffL = 0;
ISIcutoffR = max(fliplr(sort(diff(spt),'descend'))./ppms);
 
ISIs = fliplr(sort(diff(spt),'descend'))./ppms;
ISIsdom = linspace(0,max(ISIs),50);
PISIs = histc(ISIs,ISIsdom)./norm(histc(ISIs,ISIsdom));
    
ISIAreadom = linspace(0,ISIsdom(end),length(ISIsdom).*100);
PISIArea = interp1(ISIsdom,PISIs,ISIAreadom);
[tempval,tempind] = min(abs(ISIAreadom - ISIcutoffL));
PISIArea(1:tempind) = 0;
[tempval,tempind] = min(abs(ISIAreadom - ISIcutoffR));
PISIArea(tempind:end) = 0;
 
plot(ISIsdom,PISIs,'LineWidth',4,'color',[0 0 0])
hold on
area(ISIAreadom,PISIArea,'LineWidth',1)
hold off
xlabel('ISI (ms)')
ylabel('Normalized P(ISI)')
axis tight


% --- Executes during object creation, after setting all properties.
function axes6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes6
load test_spike_time_data
ppms = 10;
ISIcutoffL = 35;
ISIcutoffR = 65;
 
ISIs = fliplr(sort(diff(spt),'descend'))./ppms;
ISIsdom = linspace(0,max(ISIs),50);
PISIs = histc(ISIs,ISIsdom)./norm(histc(ISIs,ISIsdom));
    
ISIAreadom = linspace(0,ISIsdom(end),length(ISIsdom).*100);
PISIArea = interp1(ISIsdom,PISIs,ISIAreadom);
[tempval,tempind] = min(abs(ISIAreadom - ISIcutoffL));
PISIArea(1:tempind) = 0;
[tempval,tempind] = min(abs(ISIAreadom - ISIcutoffR));
PISIArea(tempind:end) = 0;
 
plot(ISIsdom,PISIs,'LineWidth',4,'color',[0 0 0])
hold on
area(ISIAreadom,PISIArea,'LineWidth',1)
hold off
xlabel('ISI (ms)')
ylabel('Normalized P(ISI)')
axis tight
