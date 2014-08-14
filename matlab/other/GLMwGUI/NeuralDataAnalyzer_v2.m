function varargout = NeuralDataAnalyzer_v2(varargin)
% NEURALDATAANALYZER_V2 MATLAB code for NeuralDataAnalyzer_v2.fig
%      NEURALDATAANALYZER_V2, by itself, creates a new NEURALDATAANALYZER_V2 or raises the existing
%      singleton*.
%
%      H = NEURALDATAANALYZER_V2 returns the handle to a new NEURALDATAANALYZER_V2 or the handle to
%      the existing singleton*.
%
%      NEURALDATAANALYZER_V2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEURALDATAANALYZER_V2.M with the given input arguments.
%
%      NEURALDATAANALYZER_V2('Property','Value',...) creates a new NEURALDATAANALYZER_V2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NeuralDataAnalyzer_v2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NeuralDataAnalyzer_v2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NeuralDataAnalyzer_v2

% Last Modified by GUIDE v2.5 16-Jan-2014 12:10:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NeuralDataAnalyzer_v2_OpeningFcn, ...
                   'gui_OutputFcn',  @NeuralDataAnalyzer_v2_OutputFcn, ...
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


% --- Executes just before NeuralDataAnalyzer_v2 is made visible.
function NeuralDataAnalyzer_v2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NeuralDataAnalyzer_v2 (see VARARGIN)

% Choose default command line output for NeuralDataAnalyzer_v2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NeuralDataAnalyzer_v2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NeuralDataAnalyzer_v2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%SAVE BUTTON%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
disp('Saving...')
basedir = pwd;

FStim = get(handles.text1,'String');
[FStimM,FStimN] = size(FStim);

FSpt = get(handles.text2,'String');
[FSptM,FSptN] = size(FSpt);

if FStim(1) == 'n' % the user did not specify the file
    errordlg('The stimulus file has not been specified','Error')
    return
elseif FSpt(1) == 'n' % the user did not specify the file
    errordlg('The spike time file has not been specified','Error')
    return
elseif FStimM == 2 % the user tried to input the data once but did not
    errordlg('The stimulus file has not been specified','Error')
    return
elseif FSptM == 2 % the user tried to input the data once but did not
    errordlg('The spike time file has not been specified','Error')
    return
else
    
    newFile = generateMATfile(FStim,FSpt); % new file name assignment for the title
    
    % repeated stimulus data file loading
    RStim = get(handles.text3,'String');
    [RStimM,RStimN] = size(RStim);
 
    % repeated spike time data file loading
    RSpt = get(handles.text4,'String');
    [RSptM,RSptN] = size(RSpt);
    
    if RStim(1) == 'n' % the user did not specify the file
    elseif RSpt(1) == 'n' % the user did not specify the file
    elseif RStimM == 2 % the user tried to input the data once but did not
    elseif RSptM == 2 % the user tried to input the data once but did not
    else
        % load the stimulus file.    
        unknownvar = load(RStim);
        fieldunknown = fieldnames(unknownvar);

        % if the input stimulus mat file contains more than 1 file, stop the
        % program and cast the error.
        if length(fieldunknown) ~= 1
            errordlg('Make sure that input data contains only one data file : Repeated stimulus data','Error')
            return
        end

        firstel = fieldunknown{1};
        rstim = getfield(unknownvar,firstel); % Stimulus file assignment
    
        % load the repeated spike time file.    
        unknownvar = load(RSpt);
        fieldunknown = fieldnames(unknownvar);

        % if the input stimulus mat file contains more than 1 file, stop the
        % program and cast the error.
        if length(fieldunknown) ~= 1
            errordlg('Make sure that input data contains only one data file : Spike time data','Error')
            return
        end  

        firstel = fieldunknown{1};
        rspt = getfield(unknownvar,firstel); % spike time file assignment
        
        RepeatedData.stim = rstim;
        RepeatedData.spike_time = rspt;
        
        save([newFile '.mat'],'RepeatedData','-append')
    end
    
    if isfield(handles,'autocorr') % see if auto correlation analysis has done
        autocorr.autocorr = handles.autocorr.autocorr;
        save([newFile '.mat'],'autocorr','-append')
    end

    if isfield(handles,'ISI') % see if ISI analysis has been done
        ISI.ISIs = handles.ISI.ISIs;
        ISI.ISIsdom = handles.ISI.ISIsdom;
        ISI.PISIs = handles.ISI.PISIs;
        save([newFile '.mat'],'ISI','-append')
    end

    if isfield(handles,'sta_analysis')  % see if STA analysis has been done
        sta_analysis.STA = handles.sta_analysis.STA;
        sta_analysis.ProjectionOntoSTA = handles.sta_analysis.Proj;
        sta_analysis.PriorDist = handles.sta_analysis.PriorDist;
        sta_analysis.FiltteredDist = handles.sta_analysis.FiltteredDist;
        sta_analysis.STAThresh = handles.sta_analysis.Thresh;
        save([newFile '.mat'],'sta_analysis','-append')
    end

    if isfield(handles,'stc_analysis') % see if STC analysis has been done
        stc_analysis.Mode1 = handles.stc_analysis.Mode1;
        stc_analysis.Mode2 = handles.stc_analysis.Mode2;
        stc_analysis.Covdiff = handles.stc_analysis.Covdiff;
        stc_analysis.ProjectionOntoMode1 = handles.stc_analysis.Proj1;
        stc_analysis.ProjectionOntoMode2 = handles.stc_analysis.Proj2;
        stc_analysis.Mode1Thresh = handles.stc_analysis.Mode1Thresh;
        stc_analysis.Mode2Thresh = handles.stc_analysis.Mode2Thresh;
        save([newFile '.mat'],'stc_analysis','-append')
    end

    if isfield(handles,'glm_analysis') % see if GLM analysis has been done
        glm_analysis.Stimfilt = handles.glm_analysis.Stimfilt;
        glm_analysis.Histfilt = handles.glm_analysis.Histfilt;
        glm_analysis.negloglival = handles.glm_analysis.negloglival;
        save([newFile '.mat'],'glm_analysis','-append')
    end

    if isfield(handles,'PSTHpredict') % see if PSTH prediction has been done
        PSTH_prediction.STApsth = handles.PSTHpredict.STApsth;
        PSTH_prediction.STCpsth = handles.PSTHpredict.STCpsth;
        PSTH_prediction.GLMpsth = handles.PSTHpredict.GLMpsth;
        save([newFile '.mat'],'PSTH_prediction','-append')
    end

    if isfield(handles,'ISIpredict') % see if ISI prediction has been done
        ISI_prediction.SimPISIs = handles.ISIpredict.SimPISI;
        ISI_prediction.STAPISIs = handles.ISIpredict.STAPISI;
        ISI_prediction.STCPISIs = handles.ISIpredict.STCPISI;
        ISI_prediction.GLMPISIs = handles.ISIpredict.GLMPISI;
        save([newFile '.mat'],'ISI_prediction','-append')
    end

    disp(['The data has been saved as ' newFile '.mat'])
end

function newfile = generateMATfile(input1,input2)
% This function generates duplicated file that contains analyzed data.
input1 = input1(1:length(input1)-4);
load(input1)
input2 = input2(1:length(input2)-4);
load(input2)
str = '_analyzed';
input1 = input1(length(pwd)+2:end);
input2 = input2(length(pwd)+2:end);
save([input1 '_' input2 str])
newfile = [input1 '_' input2 str];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%QUIT BUTTON%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close(handles.figure1)
display('Quit plotting gui!')

% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in checkbox7.
function checkbox7_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox7


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%STA ANALYSIS BUTTON%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
basedir = pwd;
addpath([basedir '/tools_revcorrelations/']);

% Stimulus data file loading
FStim = get(handles.text1,'String');
[FStimM,FStimN] = size(FStim);

% Spike time data file loading
FSpt = get(handles.text2,'String');
[FSptM,FSptN] = size(FSpt);

if FStim(1) == 'n' % the user did not specify the file
    errordlg('The stimulus file has not been specified','Error')
    return
elseif FSpt(1) == 'n' % the user did not specify the file
    errordlg('The spike time file has not been specified','Error')
    return
elseif FStimM == 2 % the user tried to input the data once but did not
    errordlg('The stimulus file has not been specified','Error')
    return
elseif FSptM == 2 % the user tried to input the data once but did not
    errordlg('The spike time file has not been specified','Error')
    return
else
    % load the stimulus file.    
    Stim = handles.stim;
    len = length(Stim);
    
    % load the spike time file.    
    spt = handles.spike_time;

    ppms = str2num(get(handles.edit1,'String')); % points per ms
    if isempty(ppms) % if ppms is empty or not an integer
        errordlg('Please insert positive integer value for a display time','Error')
        return
    elseif ppms <= 0 % ppms is not a positive integer
        errordlg('Points per ms must be positive value','Error')
        return
    end
    handles.ppms = ppms;
    
    dt = str2num(get(handles.edit2,'String'));
    if isempty(dt) % if dt is empty or not an integer
        errordlg('Please insert positive number between 0 and 1 for the sampling rate','Error')
        return
    elseif dt <= 0 % ppms is not a positive integer
        errordlg('Sampling rate must be positive value','Error')
        return
    end
    handles.dt = dt;
    
    timeb4spk = str2num(get(handles.edit3,'String'));
    if isempty(timeb4spk) % if ppms is empty or not an integer
        errordlg('Please insert positive integer value for a display time(time before the spike)','Error')
        return
    elseif timeb4spk <= 0 % ppms is not a positive integer
        errordlg('Display time(time before the spike) must be positive value','Error')
        return
    end
    handles.timeb4spk = timeb4spk;
    
    % If the user did not check spike time sorting check box or user did
    % check that but did not specify any value, just go ahead the analysis
    % with all spikes.
    All = 0;
    if get(handles.checkbox2,'Value') == 0
        All = 1;
    else
        % if the user checked the isolated spike,

        % min ISI assignment. If the user did not specify, the default value
        % will be the zero.
        ISIcutoffL = str2num(get(handles.edit6,'String'));
        if isempty(ISIcutoffL)
            errordlg('Please insert positive integer value for a minimum ISI cutoff','Error')
            return
        elseif ISIcutoffL < 0
            errordlg('Minimum ISI cutoff time must be positive value','Error')
            return
        end
        Minval = handles.ISI.ISIcutoffL;

        % max ISI assignment. If the user did not specify, the default value
        % will be the maximum ISI possible.
        cutoffR = get(handles.edit7,'String');
        if cutoffR(1) == 'm' && cutoffR(2) == 'a' && cutoffR(3) == 'x'
            ISIcutoffR = max(fliplr(sort(diff(spt),'descend'))./ppms);
        elseif isempty(cutoffR) || str2num(cutoffR) < 0
            errordlg('Please insert "max" or positive integer value for a maximum ISI cutoff','Error')
            return
        elseif ISIcutoffL >= str2num(cutoffR)
            errordlg('Right bound of ISI cutoff should be bigger than the left bound of ISI cutoff','Error')
            return
        else
            ISIcutoffR = str2num(cutoffR);
        end
        Maxval = handles.ISI.ISIcutoffR ;

        Isolated_spike_time = [];

        Isolated_spike_time(1) = spt(1);

        for ii = 2:length(spt)
            if (spt(ii) - spt(ii-1))/handles.ppms >= Minval &&...
                    (spt(ii) - spt(ii-1))/handles.ppms <= Maxval
                Isolated_spike_time = [Isolated_spike_time,spt(ii)];
            end
        end
    end
    
    FR = length(spt)/(len/(1000*handles.ppms));

    % starting STA analysis
    StaStatBar = waitbar(0,'Calculating spike triggered average');
    if All
        STA = getSTA(handles.stim,handles.spike_time,handles.timeb4spk,...
            handles.ppms,handles.dt); % feature calculation
        waitbar(0.5,StaStatBar,'Calculating spike triggered average');
        [FiltteredResponse,PriorDist,FiltteredDist] = FeatSelect(...
            handles.stim,spt,STA); % filttered response and prior calculation
    else
        STA = getSTA(handles.stim,Isolated_spike_time,...
            handles.timeb4spk,handles.ppms,dt); % feature calculation
        waitbar(0.5,StaStatBar,'Calculating spike triggered average');
        [FiltteredResponse,PriorDist,FiltteredDist] = FeatSelect(...
            handles.stim,Isolated_spike_time,STA); % filttered response and prior calculation
        handles.sta_analysis.Isolated_spike_time = Isolated_spike_time;
    end

    waitbar(1,StaStatBar,'Calculating spike triggered average');
    Thresh = NLest(FR,PriorDist,FiltteredDist); % nonlinearity
    Thresh = Thresh./norm(Thresh);

    % saving files
    handles.sta_analysis.STA = STA;
    handles.sta_analysis.Proj = FiltteredResponse;
    handles.sta_analysis.PriorDist = PriorDist;
    handles.sta_analysis.FiltteredDist = FiltteredDist;
    handles.sta_analysis.Thresh = Thresh;
    handles.FR = FR;

    % creating string objects for plot title
    FStim = FStim(1:length(FStim)-4);
    FSpt = FSpt(1:length(FSpt)-4);
    FStim = FStim(length(pwd)+2:end);
    FSpt = FSpt(length(pwd)+2:end);
    newFile = [FStim ' and ' FSpt];

    guidata(hObject,handles)

    close(StaStatBar)

    % Plotting preparation
    minFiltteresResp = min(FiltteredResponse);
    maxFiltteresResp = max(FiltteredResponse);
    LNdom = linspace(minFiltteresResp,maxFiltteresResp,length(PriorDist));

    figure(3)
    set(figure(3),'Position',[.5 .5 1000 1000])
    subplot(2,2,1)
    plot(linspace(-handles.timeb4spk,0,length(handles.sta_analysis.STA)),...
        handles.sta_analysis.STA,'LineWidth',4)
    set(gca,'FontSize',13)
    title('Spike triggered average')
    xlabel('Time (ms)')
    axis tight
    subplot(2,2,2)
    plot(LNdom,PriorDist,LNdom,FiltteredDist)
    set(gca,'FontSize',13)
    xlabel('Projections onto STA')
    legend('P(s), prior','P(s|spike)','Location','Northwest')
    subplot(2,2,3)
    set(gca,'FontSize',13)
    plot(LNdom,Thresh)
    xlabel('Projection onto STA')
    ylabel('P(spike | s)')
    title('STA threshold')
    annotation(figure(3),'textbox',...
        [0.161 0.969278033794163 0.715 0.0261136712749616],...
        'String',{[newFile ' STA analysis' ' ' datestr(now)]},...
        'HorizontalAlignment','center',...
        'FontSize',14,...
        'FitBoxToText','off',...
        'LineStyle','none');
    if All
        annotation(figure(3),'textbox',...
            [0.563000000000001 0.346738059944806 0.334 0.0967741935483871],...
            'String',{['Number of spike used : ' num2str(length(spt)) ...
            ' / ' num2str(length(spt))]},...
            'HorizontalAlignment','center',...
            'FontSize',12,...
            'FitBoxToText','off',...
            'LineStyle','none');
    else
        annotation(figure(3),'textbox',...
            [0.563000000000001 0.346738059944806 0.334 0.0967741935483871],...
            'String',{['Number of spike used : ' num2str(length(Isolated_spike_time)) ...
            ' / ' num2str(length(spt))]},...
            'HorizontalAlignment','center',...
            'FontSize',12,...
            'FitBoxToText','off',...
            'LineStyle','none');
    end

    AutoSave = get(handles.checkbox3,'Value');
    if AutoSave
        saveas(figure(3),[FStim '_' FSpt '_analyzed_STA_analysis_' datestr(now)],'fig');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%STC ANALYSIS BUTTON%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
basedir = pwd;
addpath([basedir '/tools_revcorrelations/']);

% Stimulus data file loading
FStim = get(handles.text1,'String');
[FStimM,FStimN] = size(FStim);
 
% Spike time data file loading
FSpt = get(handles.text2,'String');
[FSptM,FSptN] = size(FSpt);

if FStim(1) == 'n' % the user did not specify the file
    errordlg('The stimulus file has not been specified','Error')
    return
elseif FSpt(1) == 'n' % the user did not specify the file
    errordlg('The spike time file has not been specified','Error')
    return
elseif FStimM == 2 % the user tried to input the data once but did not
    errordlg('The stimulus file has not been specified','Error')
    return
elseif FSptM == 2 % the user tried to input the data once but did not
    errordlg('The spike time file has not been specified','Error')
    return
else
    % load the stimulus file.    
    Stim = handles.stim;
    len = length(Stim);
    
    % load the spike time file.    
    spt = handles.spike_time;
 
    ppms = str2num(get(handles.edit1,'String')); % points per ms
    if isempty(ppms) % if ppms is empty or not an integer
        errordlg('Please insert positivee integer value for the number of points per ms','Error')
        return
    elseif ppms <= 0 % ppms is not a positive integer
        errordlg('The number of points per ms must be positive value','Error')
        return
    end
    handles.ppms = ppms;
    
    dt = str2num(get(handles.edit2,'String'));
    if isempty(dt) % if dt is empty or not an integer
        errordlg('Please insert positive number between 0 and 1 for the sampling rate','Error')
        return
    elseif dt <= 0 % ppms is not a positive integer
        errordlg('Sampling rate must be positive value','Error')
        return
    end
    handles.dt = dt;
    
    timeb4spk = str2num(get(handles.edit3,'String'));
    if isempty(timeb4spk) % if ppms is empty or not an integer
        errordlg('Please insert positive integer value for a display time(time before the spike)','Error')
        return
    elseif timeb4spk <= 0 % ppms is not a positive integer
        errordlg('Display time(time before the spike) must be positive value','Error')
        return
    end
    handles.timeb4spk = timeb4spk;
    
    % If the user did not check spike time sorting check box or user did
    % check that but did not specify any value, just go ahead the analysis
    % with all spikes.
    All = 0;
    if get(handles.checkbox2,'Value') == 0
        All = 1;
    else
        % if the user checked the isolated spike,
 
        % min ISI assignment. If the user did not specify, the default value
        % will be the zero.
        ISIcutoffL = str2num(get(handles.edit6,'String'));
        if isempty(ISIcutoffL)
            errordlg('Please insert positive integer value for a minimum ISI cutoff','Error')
            return
        elseif ISIcutoffL < 0
            errordlg('Minimum ISI cutoff time must be positive value','Error')
            return
        end
        Minval = handles.ISI.ISIcutoffL;
 
        % max ISI assignment. If the user did not specify, the default value
        % will be the maximum ISI possible.
        cutoffR = get(handles.edit7,'String');
        if cutoffR(1) == 'm' && cutoffR(2) == 'a' && cutoffR(3) == 'x'
            ISIcutoffR = max(fliplr(sort(diff(spt),'descend'))./ppms);
        elseif isempty(cutoffR) || str2num(cutoffR) < 0
            errordlg('Please insert "max" or positive integer value for a maximum ISI cutoff','Error')
            return
        elseif ISIcutoffL >= str2num(cutoffR)
            errordlg('Right bound of ISI cutoff should be bigger than the left bound of ISI cutoff','Error')
            return
        else
            ISIcutoffR = str2num(cutoffR);
        end
        Maxval = handles.ISI.ISIcutoffR ;
 
        Isolated_spike_time = [];
 
        Isolated_spike_time(1) = spt(1);
 
        for ii = 2:length(spt)
            if (spt(ii) - spt(ii-1))/handles.ppms >= Minval &&...
                    (spt(ii) - spt(ii-1))/handles.ppms <= Maxval
                Isolated_spike_time = [Isolated_spike_time,spt(ii)];
            end
        end
    end
    
    FR = length(spt)/(len/(1000*handles.ppms));
    
    % starting STC analysis
    StcStatBar = waitbar(0,'Calculating covariance modes');
    if All % if we use all spikes
        [Mode1,Mode2,covdiff] = getSTC(handles.stim,handles.spike_time,...
            handles.timeb4spk,handles.ppms,handles.dt); % features calculation
        [eigenvec,eigenval] = eig(-covdiff);
        SortedEigenval = diag(eigenval);

        waitbar(0.25,StcStatBar,'Calculating covariance modes');
        [Mode1FiltteredResponse,Mode1PriorDist,Mode1FiltteredDist] = ...
            FeatSelect(handles.stim,handles.spike_time,Mode1); % filttered response and prior calculation for the first mode
        waitbar(0.5,StcStatBar,'Calculating covariance modes');

        [Mode2FiltteredResponse,Mode2PriorDist,Mode2FiltteredDist] = ...
            FeatSelect(handles.stim,handles.spike_time,Mode2); % filttered response and prior calculation for the second mode
    else
        [Mode1,Mode2,covdiff] = getSTC(handles.stim,Isolated_spike_time,...
            handles.timeb4spk,handles.ppms,handles.dt); % features calculation
        [eigenvec,eigenval] = eig(-covdiff);
        SortedEigenval = diag(eigenval);

        waitbar(0.25,StcStatBar,'Calculating covariance modes');
        [Mode1FiltteredResponse,Mode1PriorDist,Mode1FiltteredDist] = ...
            FeatSelect(handles.stim,Isolated_spike_time,Mode1); % filttered response and prior calculation for the first mode
        waitbar(0.5,StcStatBar,'Calculating covariance modes');

        [Mode2FiltteredResponse,Mode2PriorDist,Mode2FiltteredDist] = ...
            FeatSelect(handles.stim,Isolated_spike_time,Mode2); % filttered response and prior calculation for the second mode
    end
    waitbar(0.75,StcStatBar,'Calculating covariance modes');
    Mode1Thresh = NLest(FR,Mode1PriorDist,Mode1FiltteredDist); % nonlinearity
    Mode1Thresh = Mode1Thresh./norm(Mode1Thresh);
    
    Mode2Thresh = NLest(FR,Mode2PriorDist,Mode2FiltteredDist); % nonlinearity
    Mode2Thresh = Mode2Thresh./norm(Mode2Thresh);
    waitbar(1,StcStatBar,'Calculating covariance modes');
    
    % saving files
    handles.stc_analysis.Mode1 = Mode1;
    handles.stc_analysis.Mode2 = Mode2;
    handles.stc_analysis.Covdiff = covdiff;
    handles.stc_analysis.Proj1 = Mode1FiltteredResponse;
    handles.stc_analysis.Proj2 = Mode2FiltteredResponse;
    handles.stc_analysis.Mode1Thresh = Mode1Thresh;
    handles.stc_analysis.Mode2Thresh = Mode2Thresh;
    handles.FR = FR;

    % creating string objects for plot title
    FStim = FStim(1:length(FStim)-4);
    FSpt = FSpt(1:length(FSpt)-4);
    FStim = FStim(length(pwd)+2:end);
    FSpt = FSpt(length(pwd)+2:end);
    newFile = [FStim ' and ' FSpt];

    guidata(hObject,handles)

    close(StcStatBar)

    % Plotting preparation
    Mode1FSlb = min(Mode1FiltteredResponse);
    Mode1FSub = max(Mode1FiltteredResponse);
    Mode2FSlb = min(Mode2FiltteredResponse);
    Mode2FSub = max(Mode2FiltteredResponse);
    Mode1CT3dom = linspace(Mode1FSlb,Mode1FSub,100);
    Mode2CT3dom = linspace(Mode2FSlb,Mode2FSub,100);

    Mode1FiltteredDist = Mode1FiltteredDist(:);
    Mode2FiltteredDist = Mode2FiltteredDist(:);
    Mode1PriorDist = Mode1PriorDist(:);
    Mode2PriorDist = Mode2PriorDist(:);
    
    TDFiltteredDist = sqrt(Mode1FiltteredDist*Mode2FiltteredDist');
    TDPriorDist = sqrt(Mode1PriorDist*Mode2PriorDist');

    for i = 1:length(Mode1FiltteredDist)
        if Mode1FiltteredDist(i) >= max(Mode1FiltteredDist)*.1
            Mode1sflb = i;
            break
        end
    end

    for i = length(Mode1FiltteredDist):-1:1
        if Mode1FiltteredDist(i) >= max(Mode1FiltteredDist)*.1
            Mode1sfub = i;
            break
        end
    end

    for i = 1:length(Mode2FiltteredDist)
        if Mode2FiltteredDist(i) >= max(Mode2FiltteredDist)*.1
            Mode2sflb = i;
            break
        end
    end

    for i = length(Mode2FiltteredDist):-1:1
        if Mode2FiltteredDist(i) >= max(Mode2FiltteredDist)*.1
            Mode2sfub = i;
            break
        end
    end

    if (Mode1sfub-Mode1sflb)>(Mode2sfub-Mode2sflb)
        sfub = Mode1sfub;
        sflb = Mode1sflb;
    else
        sfub = Mode2sfub;
        sflb = Mode2sflb;
    end

    clear Mode1sfub Mode1sflb Mode2sfub Mode2sflb

    Mode1Thresh = Mode1Thresh(:);
    Mode2Thresh = Mode2Thresh(:);
    
    TDThresh = Mode1Thresh*Mode2Thresh';
    TDThresh = TDThresh./norm(TDThresh);

    figure(4)
    set(figure(4),'Position',[.5 .5 1200 1000])
    subplot(2,3,1)
    plot(linspace(-handles.timeb4spk,0,handles.ppms*handles.timeb4spk*dt),...
        Mode1,'b','Linewidth',4)
    hold on
    plot(linspace(-handles.timeb4spk,0,handles.ppms*handles.timeb4spk*dt),...
        Mode2,'g','Linewidth',4)
    hold off
    set(gca,'Fontsize',13)
    title('Covariance mode')
    xlabel('Time (ms)')
    axis tight
    legend('First mode','Second mode','location','Northwest')
    subplot(2,3,2)
    imagesc(covdiff)
    set(gca,'Fontsize',13,'XTick',[],'YTick',[])
    title('Covariance difference matrix')
    xlabel([num2str(handles.timeb4spk) '(ms)'])
    ylabel([num2str(handles.timeb4spk) '(ms)'])
    subplot(2,3,3)
    plot(1:length(SortedEigenval),SortedEigenval,'ko')
    set(gca,'Fontsize',13)
    title('Eigenvalue spectrum')
    xlabel('Eigenvalue number')
    ylabel('Eigenvalue')
    axis tight
    subplot(2,3,4)
    contour3(Mode1CT3dom,Mode2CT3dom,TDFiltteredDist,10);
    hold on
    contour(Mode1CT3dom,Mode2CT3dom,TDPriorDist,5);
    hold on
    mesh(Mode1CT3dom(sflb:4:sfub),Mode2CT3dom(sflb:4:sfub),...
        TDFiltteredDist(sflb:4:sfub,sflb:4:sfub),'EdgeColor',[.7 .7 .7],'FaceColor','none');
    hold off
    set(gca,'Fontsize',13)
    view(-40,30)
    xlabel('P(S_{2} | spike)')
    ylabel('P(S_{1} | spike)')
    zlabel('P(S_{1},S_{2} | spikes)')
    title('Feature selectivity')
    subplot(2,3,6)
    mesh(Mode1CT3dom(1:4:end),Mode2CT3dom(1:4:end),TDThresh(1:4:end,1:4:end))
    set(gca,'Fontsize',13)
    xlabel('Projection onto Mode 2')
    ylabel('Projection onto Mode 1')
    zlabel('P(spikes | S_{1},S_{2})')
    title('2 dimensional representation of nonlinearity')
    axis([Mode1FSlb Mode1FSub Mode2FSlb Mode2FSub 0 max(max(TDThresh))])
    annotation(figure(4),'textbox',...
        [0.161 0.969278033794163 0.715 0.0261136712749616],...
        'String',{[newFile ' Covariance analysis' ' ' datestr(now)]},...
        'HorizontalAlignment','center',...
        'FontSize',14,...
        'FitBoxToText','off',...
        'LineStyle','none');
    if All
        annotation(figure(4),'textbox',...
            [0.434333333333335 0.449943757030371 0.191499999999999 0.0391402687567301],...
            'String',{['Number of spike used : ' num2str(length(spt)) ...
            ' / ' num2str(length(spt))]},...
            'HorizontalAlignment','center',...
            'FontSize',12,...
            'FitBoxToText','off',...
            'LineStyle','none');
    else
        annotation(figure(4),'textbox',...
            [0.434333333333335 0.449943757030371 0.191499999999999 0.0391402687567301],...
            'String',{['Number of spike used : ' num2str(length(Isolated_spike_time)) ...
            ' / ' num2str(length(spt))]},...
            'HorizontalAlignment','center',...
            'FontSize',12,...
            'FitBoxToText','off',...
            'LineStyle','none');
    end

    AutoSave = get(handles.checkbox3,'Value');
    if AutoSave
        saveas(figure(4),[FStim '_' FSpt '_STC_analysis_' datestr(now)],'fig');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%GLM ANALYSIS BUTTON%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
basedir = pwd;
addpath([basedir '/tools_revcorrelations/']);
addpath([basedir '/tools_newglm/']);

% Stimulus data file loading
FStim = get(handles.text1,'String');
[FStimM,FStimN] = size(FStim);
 
% Spike time data file loading
FSpt = get(handles.text2,'String');
[FSptM,FSptN] = size(FSpt);
 
if FStim(1) == 'n' % the user did not specify the file
    errordlg('The stimulus file has not been specified','Error')
    return
elseif FSpt(1) == 'n' % the user did not specify the file
    errordlg('The spike time file has not been specified','Error')
    return
elseif FStimM == 2 % the user tried to input the data once but did not
    errordlg('The stimulus file has not been specified','Error')
    return
elseif FSptM == 2 % the user tried to input the data once but did not
    errordlg('The spike time file has not been specified','Error')
    return
else

    % load the stimulus file.    
    Stim = handles.stim;
    len = length(Stim);
    
    % load the spike time file.    
    spt = handles.spike_time;
 
    ppms = str2num(get(handles.edit1,'String')); % points per ms
    if isempty(ppms) % if ppms is empty or not an integer
        errordlg('Please insert positivee integer value for the number of points per ms','Error')
        return
    elseif ppms <= 0 % ppms is not a positive integer
        errordlg('The number of points per ms must be positive value','Error')
        return
    end
    handles.ppms = ppms;
    
    timeb4spk = str2num(get(handles.edit3,'String'));
    if isempty(timeb4spk) % if ppms is empty or not an integer
        errordlg('Please insert positive integer value for a display time(time before the spike)','Error')
        return
    elseif timeb4spk <= 0 % ppms is not a positive integer
        errordlg('Display time(time before the spike) must be positive value','Error')
        return
    end
    handles.timeb4spk = timeb4spk;

    timeafterspk = str2num(get(handles.edit4,'String'));
    if isempty(timeafterspk) % if timeafterspk is empty or not an integer
        errordlg('Please insert positive integer value for a display time(timeafterspk)','Error')
        return
    elseif timeafterspk <= 0 % timeafterspk is not a positive integer
        errordlg('Display time(timeafterspk) must be positive value','Error')
        return
    end
    handles.timeafterspk = timeafterspk;

    % If the user did not check spike time sorting check box or user did
    % check that but did not specify any value, just go ahead the analysis
    % with all spikes.
    All = 0;
    if get(handles.checkbox3,'Value') == 0
        All = 1;
    else
        % if the user checked the isolated spike,
 
        % min ISI assignment. If the user did not specify, the default value
        % will be the zero.
        ISIcutoffL = str2num(get(handles.edit6,'String'));
        if isempty(ISIcutoffL)
            errordlg('Please insert positive integer value for a minimum ISI cutoff','Error')
            return
        elseif ISIcutoffL < 0
            errordlg('Minimum ISI cutoff time must be positive value','Error')
            return
        end
        Minval = handles.ISI.ISIcutoffL;
 
        % max ISI assignment. If the user did not specify, the default value
        % will be the maximum ISI possible.
        cutoffR = get(handles.edit7,'String');
        if cutoffR(1) == 'm' && cutoffR(2) == 'a' && cutoffR(3) == 'x'
            ISIcutoffR = max(fliplr(sort(diff(spt),'descend'))./ppms);
        elseif isempty(cutoffR) || str2num(cutoffR) < 0
            errordlg('Please insert "max" or positive integer value for a maximum ISI cutoff','Error')
            return
        elseif ISIcutoffL >= str2num(cutoffR)
            errordlg('Right bound of ISI cutoff should be bigger than the left bound of ISI cutoff','Error')
            return
        else
            ISIcutoffR = str2num(cutoffR);
        end
        Maxval = handles.ISI.ISIcutoffR ;
 
        Isolated_spike_time = [];
 
        Isolated_spike_time(1) = spt(1);
 
        for ii = 2:length(spt)
            if (spt(ii) - spt(ii-1))/handles.ppms >= Minval &&...
                    (spt(ii) - spt(ii-1))/handles.ppms <= Maxval
                Isolated_spike_time = [Isolated_spike_time,spt(ii)];
            end
        end
    end
    
    FR = length(spt)/(len/(1000*handles.ppms));
    if str2num(get(handles.edit2,'String')) >= 0.1
        dt = 0.1;
    else
        dt = 0.01;
    end

    % starting GLM analysis
    DataStruct.timeb4spk = handles.timeb4spk;
    DataStruct.timeafterspk = handles.timeafterspk;
    DataStruct.dt = dt;
    DataStruct.spike_time = handles.spike_time;
    DataStruct.stim = handles.stim;
    DataStruct.ppms = handles.ppms;

    [Stimfilt,Histfilt,negloglival] = GLM_gen(DataStruct);

    % saving files
    handles.glm_analysis.Stimfilt = Stimfilt;
    handles.glm_analysis.Histfilt = Histfilt;
    handles.glm_analysis.negloglival = negloglival;
    handles.FR = FR;

    % creating string objects for plot title
    FStim = FStim(1:length(FStim)-4);
    FSpt = FSpt(1:length(FSpt)-4);
    FStim = FStim(length(pwd)+2:end);
    FSpt = FSpt(length(pwd)+2:end);
    newFile = [FStim ' and ' FSpt];
 
    guidata(hObject,handles)

    figure(5)
    set(figure(5),'Position',[.5 .5 1000 400])
    subplot(1,2,1)
    plot(linspace(-handles.timeb4spk,0,length(Stimfilt)),Stimfilt,'LineWidth',...
        4,'Color',[0 0 0])
    set(gca,'FontSize',13)
    title('Stimulus filter')
    xlabel('Time (ms)')
    axis tight
    subplot(1,2,2)
    plot(linspace(0,handles.timeafterspk,length(Histfilt)),Histfilt,'LineWidth',...
        4,'Color',[0 0 0])
    set(gca,'FontSize',13)
    title('Spike history filter')
    xlabel('Time (ms)')
    axis tight
    annotation(figure(5),'textbox',...
        [0.161 0.969278033794163 0.715 0.0261136712749616],...
        'String',{[newFile ' GLM analysis' ' ' datestr(now)]},...
        'HorizontalAlignment','center',...
        'FontSize',14,...
        'FitBoxToText','off',...
        'LineStyle','none');
    if All
        annotation(figure(5),'textbox',[0.403156035328749 0.00476190476190476 0.235705593719336 0.064404761904762],...
            'String',{['Number of spike used : ' num2str(length(handles.spike_time))...
            ' / ' num2str(length(handles.spike_time))]},...
            'HorizontalAlignment','center',...
            'FontSize',13,...
            'FitBoxToText','off',...
            'LineStyle','none');
    else
        annotation(figure(5),'textbox',[0.403156035328748 0.014933090024331 0.231780176643767 0.04],...
            'String',{['Number of spike used : ' num2str(length(Isolated_spike_time))...
            ' / ' num2str(length(handles.spike_time))]},...
            'HorizontalAlignment','center',...
            'FontSize',13,...
            'FitBoxToText','off',...
            'LineStyle','none');
    end
    AutoSave = get(handles.checkbox3,'Value');
    if AutoSave
        saveas(figure(5),[FStim '_' FSpt '_GLM_analysis_' datestr(now)],'fig');
    end
end


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%GLM ADVANCE SETTINGS BUTTON%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
basedir = pwd;
addpath([basedir '/tools_subfig/']);
addpath([basedir '/tools_newglm/']);

strTAS = get(handles.edit4,'String');
strPPMS = get(handles.edit1,'String');

if isempty(str2num(strTAS))
    errordlg('Please specify the time after spike value to display before GLM advance setting','Error')
    return
else
    global TAS PPMS Basispars flagval
    flagval = 1;
    TAS = str2num(strTAS);
    PPMS = str2num(strPPMS);
    Basispars = [8,0.9,4];
    help_GLMadv
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%FIRST STIMULUS ENTERING BUTTON%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname] = uigetfile('*.mat','Select the stimulus data file');
FirstStimPath = [pathname,filename];

[FSPm,FSPn] = size(FirstStimPath);

if FSPn == 2 % user did not specify the input stimulus file
    % do nothing
else
    % load the stimulus file.    
    unknownvar = load(FirstStimPath);
    fieldunknown = fieldnames(unknownvar);

    % if the input stimulus mat file contains more than 1 file, stop the
    % program and cast the error.
    if length(fieldunknown) ~= 1
        errordlg('Make sure that input data contains only one data file : Stimulus data','Error')
        return
    end
    
    firstel = fieldunknown{1};
    Stim = getfield(unknownvar,firstel); % Stimulus file assignment

    handles.stim = Stim;

    handles.FirstStimPath = FirstStimPath;
    guidata(hObject,handles);
    set(handles.text1,'String',FirstStimPath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%FIRST SPIKE TIME ENTERING BUTTON%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname] = uigetfile('*.mat','Select the spike time data file');
FirstSptPath = [pathname,filename];

[FSPm,FSPn] = size(FirstSptPath);

if FSPn == 2 % user did not specify the input stimulus file
    % do nothing
else
    % load the stimulus file.    
    unknownvar = load(FirstSptPath);
    fieldunknown = fieldnames(unknownvar);

    % if the input stimulus mat file contains more than 1 file, stop the
    % program and cast the error.
    if length(fieldunknown) ~= 1
        errordlg('Make sure that input data contains only one data file : Spike time data','Error')
        return
    end  

    firstel = fieldunknown{1};
    spt = getfield(unknownvar,firstel); % spike time file assignment
    
    handles.spike_time = spt;
    
    handles.FirstSptPath = FirstSptPath;
    guidata(hObject,handles);
    set(handles.text2,'String',FirstSptPath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%REPEATED INPUT/OUTPUT ENTERING BUTTON%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
if get(handles.checkbox1,'Value') == 1
    set(handles.pushbutton3,'enable','on');
    set(handles.text3,'enable','on');
    set(handles.pushbutton4,'enable','on');
    set(handles.text4,'enable','on');
else
    set(handles.pushbutton3,'enable','off');
    set(handles.text3,'enable','off');
    set(handles.text3,'String','no file selected');
    set(handles.pushbutton4,'enable','off');
    set(handles.text4,'enable','off');
    set(handles.text4,'String','no file selected');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%REPEATED STIMULUS ENTERING BUTTON%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname] = uigetfile('*.mat','Select the stimulus data file');
ReapeatedStimPath = [pathname,filename];

[RSPm,RSPn] = size(ReapeatedStimPath);

if RSPn == 2 % user did not specify the input stimulus file
    % do nothing
else
    % load the stimulus file.    
    unknownvar = load(ReapeatedStimPath);
    fieldunknown = fieldnames(unknownvar);

    % if the input stimulus mat file contains more than 1 file, stop the
    % program and cast the error.
    if length(fieldunknown) ~= 1
        errordlg('Make sure that input data contains only one data file : Repeated stimulus data','Error')
        return
    end
    
    firstel = fieldunknown{1};
    rstim = getfield(unknownvar,firstel); % Stimulus file assignment

    handles.ReapeatedStim = rstim;

    handles.ReapeatedStimPath = ReapeatedStimPath;
    guidata(hObject,handles);
    set(handles.text3,'String',ReapeatedStimPath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%REPEATED SPIKE TIME ENTERING BUTTON%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename,pathname] = uigetfile('*.mat','Select the spike time data file');
RepeatedSptPath = [pathname,filename];

[RSPm,RSPn] = size(RepeatedSptPath);

if RSPn == 2 % user did not specify the input stimulus file
    % do nothing
else
    % load the spike time file.    
    unknownvar = load(RepeatedSptPath);
    fieldunknown = fieldnames(unknownvar);

    % if the input stimulus mat file contains more than 1 file, stop the
    % program and cast the error.
    if length(fieldunknown) ~= 1
        errordlg('Make sure that input data contains only one data file : Spike time data','Error')
        return
    end  

    firstel = fieldunknown{1};
    rspt = getfield(unknownvar,firstel); % spike time file assignment
    
    handles.RepeatedSpt = rspt;
    
    handles.RepeatedSptPath = RepeatedSptPath;
    guidata(hObject,handles);
    set(handles.text4,'String',RepeatedSptPath);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%USE ISOLATED SPIKES CHECKBOX%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2
if get(handles.checkbox2,'Value') == 1
    set(handles.text15,'enable','on');
    set(handles.edit6,'enable','on');
    set(handles.edit7,'enable','on');
    set(handles.text16,'enable','on');
    set(handles.text17,'enable','on');
    set(handles.pushbutton7,'enable','on');
else
    set(handles.text15,'enable','off');
    set(handles.edit6,'enable','off');
    set(handles.edit7,'enable','off');
    set(handles.text16,'enable','off');
    set(handles.text17,'enable','off');
    set(handles.pushbutton7,'enable','off');
    set(handles.edit6,'String','0');
    set(handles.edit7,'String','max');
end

function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


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



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%ISI HELP BUTTON%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
basedir = pwd;
addpath([basedir '/tools_subfig/']);
help_ISIs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%INPUT ISI DISTRIBUTION BUTTON%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
basedir = pwd;
addpath([basedir '/tools_revcorrelations/']);

% Spike time data file loading
FSpt = get(handles.text2,'String');
[FSptM,FSptN] = size(FSpt);

if FSpt(1) == 'n' % the user did not specify the file
    errordlg('The spike time file has not been specified','Error')
    return
elseif FSptM == 2 % the user tried to input the data once but did not
    errordlg('The spike time file has not been specified','Error')
    return
else
    % load the spike time file.    
    spt = handles.spike_time;
    
    % file sampling rate
    ppms = str2num(get(handles.edit1,'String'));
    if isempty(ppms) % if ppms is empty or not an integer
        errordlg('Please insert positive integer value for number of points per ms','Error')
        return
    elseif ppms <= 0 % ppms is not a positive integer
        errordlg('Number of points per ms must be positive value','Error')
        return
    end
    
    handles.ppms = ppms;
    
    % minimum spike time used
    ISIcutoffL = str2num(get(handles.edit6,'String'));
    if isempty(ISIcutoffL)
        errordlg('Please insert positive integer value for a minimum ISI cutoff','Error')
        return
    elseif ISIcutoffL < 0
        errordlg('Minimum ISI cutoff time must be positive value','Error')
        return
    end
    
    % maximum spike time used
    cutoffR = get(handles.edit7,'String');
    if cutoffR(1) == 'm' && cutoffR(2) == 'a' && cutoffR(3) == 'x'
        ISIcutoffR = max(fliplr(sort(diff(spt),'descend'))./ppms);
    elseif isempty(cutoffR) || str2num(cutoffR) < 0
        errordlg('Please insert "max" or positive integer value for a maximum ISI cutoff','Error')
        return
    elseif ISIcutoffL >= str2num(cutoffR)
        errordlg('Right bound of ISI cutoff should be bigger than the left bound of ISI cutoff','Error')
        return
    else
        ISIcutoffR = str2num(cutoffR);
    end
    
    handles.ISI.ISIcutoffL = ISIcutoffL;
    handles.ISI.ISIcutoffR = ISIcutoffR;
    
    % interspike interval distribution calculation
    ISIs = fliplr(sort(diff(spt),'descend'))./ppms;
    ISIsdom = linspace(0,max(ISIs),50);
    PISIs = histc(ISIs,ISIsdom)./norm(histc(ISIs,ISIsdom));
    
    % area plotting for the identical interspike interval distribution
    % using specified spike time cutoff
    ISIAreadom = linspace(0,ISIsdom(end),length(ISIsdom).*100);
    PISIArea = interp1(ISIsdom,PISIs,ISIAreadom);
    [tempval,tempind] = min(abs(ISIAreadom - ISIcutoffL));
    PISIArea(1:tempind) = 0;
    [tempval,tempind] = min(abs(ISIAreadom - ISIcutoffR));
    PISIArea(tempind:end) = 0;
    
    handles.ISI.ISIs = ISIs;
    handles.ISI.ISIsdom = ISIsdom;
    handles.ISI.PISIs = PISIs;
    handles.ISI.ISIAreadom = ISIAreadom;
    handles.ISI.PISIArea = PISIArea;
    
    % new file name assignment for the title
    newFSpt = FSpt(1:end-4);
    newFSpt = newFSpt(length(pwd)+2:end);

    guidata(hObject,handles)

    figure(2)
    set(figure(2),'Position',[.5 .5 1000 400])
    subplot(1,2,1)
    plot(ISIsdom,PISIs,'LineWidth',4,'color',[0 0 0])
    hold on
    area(ISIAreadom,PISIArea,'LineWidth',1)
    hold off
    set(gca,'FontSize',13)
    xlabel('ISI (ms)')
    ylabel('Normalized P(ISI)')
    subplot(1,2,2)
    semilogy(ISIsdom,PISIs,'LineWidth',4)
    set(gca,'FontSize',13)
    xlabel('ISI (ms)')
    ylabel('log(Normalized P(ISI))')
    annotation(figure(2),'textbox',...
        [0.399429833169773 0.955775035863101 0.242375858684985 0.0399061032863854],...
        'String',{[newFSpt ' ISI distribution' ' ' datestr(now)]},...
        'HorizontalAlignment','center',...
        'FontSize',13,...
        'FitBoxToText','off',...
        'LineStyle','none');
end


function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%AUTO CORRELATION FUNCTION BUTTON%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
basedir = pwd;
addpath([basedir '/tools_revcorrelations/']);

% waitbar loading
AutoStatBar = waitbar(0,'Calculating autocorrelation');

% Stimulus data file loading
FStim = get(handles.text1,'String');
[FStimM,FStimN] = size(FStim);

if FStim(1) == 'n' % the user did not specify the file
    close(AutoStatBar)
    errordlg('The stimulus file has not been specified','Error')
    return
elseif FStimM == 2 % the user tried to input the data once but did not
    close(AutoStatBar)
    errordlg('The stimulus file has not been specified','Error')
    return
else
    % load the stimulus file.    
    Stim = handles.stim;
    
    disptime = str2num(get(handles.edit5,'String')); % display time
    ppms = str2num(get(handles.edit1,'String')); % points per ms
    
    handles.ppms = ppms;
    
    if isempty(disptime) % if display time is empty or not an integer
        close(AutoStatBar)
        errordlg('Please insert positivee integer value for a display time','Error')
        return
    elseif disptime <= 0 % display time is not a positive integer
        close(AutoStatBar)
        errordlg('Display time must be positive value','Error')
        return
    end
    
    if isempty(ppms) % if ppms is empty or not an integer
        close(AutoStatBar)
        errordlg('Please insert positivee integer value for the points per ms','Error')
        return
    elseif ppms <= 0 % ppms is not a positive integer
        close(AutoStatBar)
        errordlg('Number of points per ms must be positive value','Error')
        return
    end
    
    disptime = disptime.*ppms;

    waitbar(0.5,AutoStatBar,'Calculating autocorrelation');
    
    % autocorrelation of stimulus and normalization
    AutoCorr = xcorr(Stim-mean(Stim),Stim-mean(Stim),disptime)./...
        norm(xcorr(Stim,Stim,disptime));
    AutoCorr = AutoCorr(1:end-1);

    waitbar(1,AutoStatBar,'Calculating autocorrelation');

    % updating handle structure for the new file 
    handles.autocorr.autocorr = AutoCorr;

    % closing the statbar 
    close(AutoStatBar)

    % updating handle structure
    guidata(hObject,handles)

    % plotting the autocorrelation function
    figure(1)
    plot(linspace(-disptime/ppms,disptime/ppms,length(AutoCorr)),AutoCorr,...
        'LineWidth',4)
    set(gca,'FontSize',13)
    xlabel('Time (ms)')
    title(['Normalized autocorrelation ' datestr(now)])
end


function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%ISI PREDICTION BUTTON%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
basedir = pwd;
addpath([basedir '/tools_revcorrelations/']);
 
if get(handles.checkbox1,'Value') % if user specified the new I/O set for the prediction
    
    % new stimulus data file loading
    RStim = get(handles.text3,'String');
    [RStimM,RStimN] = size(RStim);
 
    % new spike time data file loading
    RSpt = get(handles.text4,'String');
    [RSptM,RSptN] = size(RSpt);
    
    if RStim(1) == 'n' % the user did not specify the file
        errordlg('The repeated stimulus file has not been specified','Error')
        return
    elseif RSpt(1) == 'n' % the user did not specify the file
        errordlg('The repeated spike time file has not been specified','Error')
        return
    elseif RStimM == 2 % the user tried to input the data once but did not
        errordlg('The repeated stimulus file has not been specified','Error')
        return
    elseif RSptM == 2 % the user tried to input the data once but did not
        errordlg('The repeated spike time file has not been specified','Error')
        return
    else
        clear RStimM RStimN RSptM RSptN
        
        rstim = handles.ReapeatedStim;
    
        handles.RepeatedStimPath = RStim;
        
        rspt = handles.RepeatedSpt;
    
        handles.RepeatedSptPath = RSpt;
        
        ppms = str2num(get(handles.edit1,'String')); % points per ms
        if isempty(ppms) % if ppms is empty or not an integer
            errordlg('Please insert positivee integer value for the number of points per ms','Error')
            return
        elseif ppms <= 0 % ppms is not a positive integer
            errordlg('The number of points per ms must be positive value','Error')
            return
        end
        handles.ppms = ppms;

        dt = str2num(get(handles.edit2,'String'));
        if isempty(dt) % if dt is empty or not an integer
            errordlg('Please insert positive number between 0 and 1 for the sampling rate','Error')
            return
        elseif dt <= 0 % ppms is not a positive integer
            errordlg('Sampling rate must be positive value','Error')
            return
        end
        handles.dt = dt;

        timeb4spk = str2num(get(handles.edit3,'String'));
        if isempty(timeb4spk) % if ppms is empty or not an integer
            errordlg('Please insert positive integer value for a display time(time before the spike)','Error')
            return
        elseif timeb4spk <= 0 % ppms is not a positive integer
            errordlg('Display time(time before the spike) must be positive value','Error')
            return
        end
        handles.timeb4spk = timeb4spk;
 
        timeafterspk = str2num(get(handles.edit4,'String'));
        if isempty(timeafterspk) % if timeafterspk is empty or not an integer
            errordlg('Please insert positive integer value for a display time(timeafterspk)','Error')
            return
        elseif timeafterspk <= 0 % timeafterspk is not a positive integer
            errordlg('Display time(timeafterspk) must be positive value','Error')
            return
        end
        handles.timeafterspk = timeafterspk;
    
        handles.PredBinSize = str2num(get(handles.edit11,'String'));
        if isempty(handles.PredBinSize)
            errordlg('Please insert integer value for a sampling bin size','Error')
            return
        elseif handles.PredBinSize <= 0
            errordlg('Sampling bin size must be positive value','Error')
            return
        end
        
        % creating string objects for plot title
        RStim = RStim(1:length(RStim)-4);
        RSpt = RSpt(1:length(RSpt)-4);
        RStim = RStim(length(pwd)+2:end);
        RSpt = RSpt(length(pwd)+2:end);
        newRFile = [RStim ' and ' RSpt];
        
        BinSize = handles.PredBinSize*handles.ppms;
        SimISIs = fliplr(sort(diff(handles.RepeatedSpt),'descend'))./handles.ppms;
        NumbBin = round(length(SimISIs)./BinSize);
        SimISIsdom = linspace(0,max(SimISIs),NumbBin);
        SimPISIs = histc(SimISIs,SimISIsdom);
        SimPISIs = SimPISIs./norm(SimPISIs); % normalization
        
        Bar = waitbar(0,'Calculating ISI distribution');
        
        % what models user selected?
        ModelFlag = [get(handles.checkbox5,'Value'),get(handles.checkbox6,'Value'),...
            get(handles.checkbox7,'Value')];
 
        if ModelFlag(1) == 1 && ~isfield(handles,'sta_analysis')
            close(Bar)
            errordlg('No STA analysis has been done','Error')
            return
        elseif ModelFlag(2) == 1 && ~isfield(handles,'stc_analysis')
            close(Bar)
            errordlg('No STC analysis has been done','Error')
            return
        elseif ModelFlag(3) == 1 && ~isfield(handles,'glm_analysis')
            close(Bar)
            errordlg('No GLM analysis has been done','Error')
            return
        end
        
        handles.ISIpredict.SimPISI = SimPISIs;
        handles = ISIgen(rstim,rspt,newRFile,NumbBin,handles,Bar,ModelFlag);
        guidata(hObject,handles)
    end
else % if user did not specify the repeated input/output, go ahead with existing data
    Q = questdlg('No repeated input/output set has been specified. Initially given input/output set will be used for ISI distribution plotting. Would you like to proceed?','PSTH');
    switch Q
        case 'Yes'
            % Stimulus data file loading
            FStim = get(handles.text1,'String');
            [FStimM,FStimN] = size(FStim);

            % Spike time data file loading
            FSpt = get(handles.text2,'String');
            [FSptM,FSptN] = size(FSpt);

            if FStim(1) == 'n' % the user did not specify the file
                errordlg('The stimulus file has not been specified','Error')
                return
            elseif FSpt(1) == 'n' % the user did not specify the file
                errordlg('The spike time file has not been specified','Error')
                return
            elseif FStimM == 2 % the user tried to input the data once but did not
                errordlg('The stimulus file has not been specified','Error')
                return
            elseif FSptM == 2 % the user tried to input the data once but did not
                errordlg('The spike time file has not been specified','Error')
                return
            else
                clear FStimM FStimN FSptM FSptN
                
                stim = handles.stim;
                spt = handles.spike_time;
                
                ppms = str2num(get(handles.edit1,'String')); % points per ms
                if isempty(ppms) % if ppms is empty or not an integer
                    errordlg('Please insert positivee integer value for the number of points per ms','Error')
                    return
                elseif ppms <= 0 % ppms is not a positive integer
                    errordlg('The number of points per ms must be positive value','Error')
                    return
                end
                handles.ppms = ppms;

                dt = str2num(get(handles.edit2,'String'));
                if isempty(dt) % if dt is empty or not an integer
                    errordlg('Please insert positive number between 0 and 1 for the sampling rate','Error')
                    return
                elseif dt <= 0 % ppms is not a positive integer
                    errordlg('Sampling rate must be positive value','Error')
                    return
                end
                handles.dt = dt;

                timeb4spk = str2num(get(handles.edit3,'String'));
                if isempty(timeb4spk) % if ppms is empty or not an integer
                    errordlg('Please insert positive integer value for a display time(time before the spike)','Error')
                    return
                elseif timeb4spk <= 0 % ppms is not a positive integer
                    errordlg('Display time(time before the spike) must be positive value','Error')
                    return
                end
                handles.timeb4spk = timeb4spk;

                timeafterspk = str2num(get(handles.edit4,'String'));
                if isempty(timeafterspk) % if timeafterspk is empty or not an integer
                    errordlg('Please insert positive integer value for a display time(timeafterspk)','Error')
                    return
                elseif timeafterspk <= 0 % timeafterspk is not a positive integer
                    errordlg('Display time(timeafterspk) must be positive value','Error')
                    return
                end
                handles.timeafterspk = timeafterspk;

                handles.PredBinSize = str2num(get(handles.edit11,'String'));
                if isempty(handles.PredBinSize)
                    errordlg('Please insert integer value for a sampling bin size','Error')
                    return
                elseif handles.PredBinSize <= 0
                    errordlg('Sampling bin size must be positive value','Error')
                    return
                end

                % creating string objects for plot title
                FStim = FStim(1:length(FStim)-4);
                FSpt = FSpt(1:length(FSpt)-4);
                FStim = FStim(length(pwd)+2:end);
                FSpt = FSpt(length(pwd)+2:end);
                newFFile = [FStim ' and ' FSpt];

                BinSize = handles.PredBinSize*handles.ppms;
                SimISIs = fliplr(sort(diff(spt),'descend'))./handles.ppms;
                NumbBin = round(length(SimISIs)./BinSize);
                SimISIsdom = linspace(0,max(SimISIs),NumbBin);
                SimPISIs = histc(SimISIs,SimISIsdom);
                SimPISIs = SimPISIs./norm(SimPISIs); % normalization

                Bar = waitbar(0,'Calculating ISI distribution');

                % what models user selected?
                ModelFlag = [get(handles.checkbox5,'Value'),get(handles.checkbox6,'Value'),...
                    get(handles.checkbox7,'Value')];

                if ModelFlag(1) == 1 && ~isfield(handles,'sta_analysis')
                    close(Bar)
                    errordlg('No STA analysis has been done','Error')
                    return
                elseif ModelFlag(2) == 1 && ~isfield(handles,'stc_analysis')
                    close(Bar)
                    errordlg('No STC analysis has been done','Error')
                    return
                elseif ModelFlag(3) == 1 && ~isfield(handles,'glm_analysis')
                    close(Bar)
                    errordlg('No GLM analysis has been done','Error')
                    return
                end

                handles.ISIpredict.SimPISI = SimPISIs;
                handles = ISIgen(stim,spt,newFFile,NumbBin,handles,Bar,ModelFlag);
                guidata(hObject,handles)
            end
        case 'No'
            return
        case 'Cancel'
            return
    end
end


function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%RASTER AND PSTH ANALYSIS BUTTON%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
basedir = pwd;
addpath([basedir '/tools_revcorrelations/']);

if get(handles.checkbox1,'Value') % if user specified the new I/O set for the prediction
    
    % new stimulus data file loading
    RStim = get(handles.text3,'String');
    [RStimM,RStimN] = size(RStim);
 
    % new spike time data file loading
    RSpt = get(handles.text4,'String');
    [RSptM,RSptN] = size(RSpt);
    
    if RStim(1) == 'n' % the user did not specify the file
        errordlg('The repeated stimulus file has not been specified','Error')
        return
    elseif RSpt(1) == 'n' % the user did not specify the file
        errordlg('The repeated spike time file has not been specified','Error')
        return
    elseif RStimM == 2 % the user tried to input the data once but did not
        errordlg('The repeated stimulus file has not been specified','Error')
        return
    elseif RSptM == 2 % the user tried to input the data once but did not
        errordlg('The repeated spike time file has not been specified','Error')
        return
    else
        clear RStimM RStimN RSptM RSptN
        
        % load the repeated stimulus file.    
        unknownvar = load(RStim);
        fieldunknown = fieldnames(unknownvar);

        % if the repeated input stimulus mat file contains more than 1 
        % file, stop the program and cast the error.
        if length(fieldunknown) ~= 1
            errordlg('Make sure that repeated input data contains only one data file : Stimulus data','Error')
            return
        end  

        firstel = fieldunknown{1};
        rstim = getfield(unknownvar,firstel); % repeated stimulus file assignment
    
        handles.RepeatedStim = rstim;
    
        handles.RepeatedStimPath = RStim;
        
        % load the repeated spike time file.    
        unknownvar = load(RSpt);
        fieldunknown = fieldnames(unknownvar);

        % if the repeated spike time mat file contains more than 1 
        % file, stop the program and cast the error.
        if length(fieldunknown) ~= 1
            errordlg('Make sure that repeated spike time data contains only one data file : Spike time data','Error')
            return
        end  

        firstel = fieldunknown{1};
        rspt = getfield(unknownvar,firstel); % repeated stimulus file assignment
    
        handles.RepeatedSpt = rspt;
    
        handles.RepeatedSptPath = RSpt;
        
        ppms = str2num(get(handles.edit1,'String')); % points per ms
        if isempty(ppms) % if ppms is empty or not an integer
            errordlg('Please insert positivee integer value for the number of points per ms','Error')
            return
        elseif ppms <= 0 % ppms is not a positive integer
            errordlg('The number of points per ms must be positive value','Error')
            return
        end
        handles.ppms = ppms;

        dt = str2num(get(handles.edit2,'String'));
        if isempty(dt) % if dt is empty or not an integer
            errordlg('Please insert positive number between 0 and 1 for the sampling rate','Error')
            return
        elseif dt <= 0 % ppms is not a positive integer
            errordlg('Sampling rate must be positive value','Error')
            return
        end
        handles.dt = dt;

        timeb4spk = str2num(get(handles.edit3,'String'));
        if isempty(timeb4spk) % if ppms is empty or not an integer
            errordlg('Please insert positive integer value for a display time(time before the spike)','Error')
            return
        elseif timeb4spk <= 0 % ppms is not a positive integer
            errordlg('Display time(time before the spike) must be positive value','Error')
            return
        end
        handles.timeb4spk = timeb4spk;
 
        timeafterspk = str2num(get(handles.edit4,'String'));
        if isempty(timeafterspk) % if timeafterspk is empty or not an integer
            errordlg('Please insert positive integer value for a display time(timeafterspk)','Error')
            return
        elseif timeafterspk <= 0 % timeafterspk is not a positive integer
            errordlg('Display time(timeafterspk) must be positive value','Error')
            return
        end
    
        handles.timeafterspk = timeafterspk;
    
        % Take the right repetition value or format from the input
        handles.NumbRep = str2num(get(handles.edit8,'String'));
        if isempty(handles.NumbRep)
            errordlg('Please insert integer value for a number of repetition','Error')
            return
        elseif handles.NumbRep <= 0
            errordlg('Number of repetition must be positive value','Error')
            return
        end
                
        % creating string objects for plot title
        RStim = RStim(1:length(RStim)-4);
        RSpt = RSpt(1:length(RSpt)-4);
        RStim = RStim(length(pwd)+2:end);
        RSpt = RSpt(length(pwd)+2:end);
        newRFile = [RStim ' and ' RSpt];
    
        Rep = handles.NumbRep;
        
        % take the left bound value or format from the input
        handles.LTime = str2num(get(handles.edit9,'String'));
        if isempty(handles.LTime)
            errordlg('Please insert integer value for starting time to display','Error')
            return
        elseif handles.NumbRep < 0
            errordlg('Starting time to display must be positive value','Error')
            return
        end

        % take the right bound value or format from the input
        handles.RTime = str2num(get(handles.edit10,'String'));
        if isempty(handles.RTime)
            errordlg('Please insert integer value for ending time to display','Error')
        elseif handles.RTime <= 0
            errordlg('Ending time to display must be positive value','Error')
        elseif handles.RTime.*1000*handles.ppms > length(rstim)
            errordlg('Right bound of time to display is longer than the existing data','Error')
        end

        % left display time handling
        SortedSptLInd = 0;
        if handles.LTime == 0 % if the user wants to see from the first.
            Lpos = 1;
            SortedSptLInd = 1;
        else
            Lpos = (handles.ppms.*handles.LTime.*1000)+1;
        end
        Rpos = handles.ppms.*handles.RTime.*1000;

        % sorting the stimulus in terms of the time display user specified
        SortedStim = rstim(Lpos:Rpos);

        % sorting the spike time data in terms of the time display user specified
        if SortedSptLInd == 0
            for i = 1:length(rspt)
                if rspt(i) > Lpos
                    break
                end
            end
            SortedSptLInd = i;
        end
        if Rpos ~= length(rstim)
            for i = SortedSptLInd:length(rspt)
                if rspt(i) > Rpos
                    break
                end
            end
            SortedSptRInd = i-1;
        else
            SortedSptRInd = length(rspt);
        end
        Sortedspt = rspt(SortedSptLInd:SortedSptRInd);

        % status bar
        Bar = waitbar(0,'Calculating psth');
        
        % what models user selected?
        ModelFlag = [get(handles.checkbox5,'Value'),get(handles.checkbox6,'Value'),...
            get(handles.checkbox7,'Value')];

        if ModelFlag(1) == 1 && ~isfield(handles,'sta_analysis')
            close(Bar)
            errordlg('No STA analysis has been done','Error')
            return
        elseif ModelFlag(2) == 1 && ~isfield(handles,'stc_analysis')
            close(Bar)
            errordlg('No STC analysis has been done','Error')
            return
        elseif ModelFlag(3) == 1 && ~isfield(handles,'glm_analysis')
            close(Bar)
            errordlg('No GLM analysis has been done','Error')
            return
        end
        
        % psth function
        handles = PSTHgen(SortedStim,Sortedspt,handles,newRFile,Rpos,Lpos,Bar,ModelFlag);

        % updating handles object
        guidata(hObject,handles);
    end
else
    Q = questdlg('No repeated input/output set has been specified. Initially given input/output set will be used for PSTH prediction plotting. Would you like to proceed?','PSTH');
    switch Q
        case 'Yes'
            
            % Stimulus data file loading
            FStim = get(handles.text1,'String');
            [FStimM,FStimN] = size(FStim);

            % Spike time data file loading
            FSpt = get(handles.text2,'String');
            [FSptM,FSptN] = size(FSpt);

            if FStim(1) == 'n' % the user did not specify the file
                errordlg('The stimulus file has not been specified','Error')
                return
            elseif FSpt(1) == 'n' % the user did not specify the file
                errordlg('The spike time file has not been specified','Error')
                return
            elseif FStimM == 2 % the user tried to input the data once but did not
                errordlg('The stimulus file has not been specified','Error')
                return
            elseif FSptM == 2 % the user tried to input the data once but did not
                errordlg('The spike time file has not been specified','Error')
                return
            else
                clear FStimM FStimN FSptM FSptN

                % load the stimulus file.    
                stim = handles.stim;

                % load the spike time file.    
                spt = handles.spike_time;

                ppms = str2num(get(handles.edit1,'String')); % points per ms
                if isempty(ppms) % if ppms is empty or not an integer
                    errordlg('Please insert positivee integer value for the number of points per ms','Error')
                    return
                elseif ppms <= 0 % ppms is not a positive integer
                    errordlg('The number of points per ms must be positive value','Error')
                    return
                end
                handles.ppms = ppms;

                dt = str2num(get(handles.edit2,'String'));
                if isempty(dt) % if dt is empty or not an integer
                    errordlg('Please insert positive number between 0 and 1 for the sampling rate','Error')
                    return
                elseif dt <= 0 % ppms is not a positive integer
                    errordlg('Sampling rate must be positive value','Error')
                    return
                end
                handles.dt = dt;

                timeb4spk = str2num(get(handles.edit3,'String'));
                if isempty(timeb4spk) % if ppms is empty or not an integer
                    errordlg('Please insert positive integer value for a display time(time before the spike)','Error')
                    return
                elseif timeb4spk <= 0 % ppms is not a positive integer
                    errordlg('Display time(time before the spike) must be positive value','Error')
                    return
                end
                handles.timeb4spk = timeb4spk;

                timeafterspk = str2num(get(handles.edit4,'String'));
                if isempty(timeafterspk) % if timeafterspk is empty or not an integer
                    errordlg('Please insert positive integer value for a display time(timeafterspk)','Error')
                    return
                elseif timeafterspk <= 0 % timeafterspk is not a positive integer
                    errordlg('Display time(timeafterspk) must be positive value','Error')
                    return
                end
                handles.timeafterspk = timeafterspk;

                % Take the right repetition value or format from the input
                handles.NumbRep = str2num(get(handles.edit8,'String'));
                if isempty(handles.NumbRep)
                    errordlg('Please insert integer value for a number of repetition','Error')
                    return
                elseif handles.NumbRep <= 0
                    errordlg('Number of repetition must be positive value','Error')
                    return
                end

                % creating string objects for plot title
                FStim = FStim(1:length(FStim)-4);
                FSpt = FSpt(1:length(FSpt)-4);
                FStim = FStim(length(pwd)+2:end);
                FSpt = FSpt(length(pwd)+2:end);
                newFFile = [FStim ' and ' FSpt];

                Rep = handles.NumbRep;

                % take the left bound value or format from the input
                handles.LTime = str2num(get(handles.edit9,'String'));
                if isempty(handles.LTime)
                    errordlg('Please insert integer value for starting time to display','Error')
                    return
                elseif handles.NumbRep < 0
                    errordlg('Starting time to display must be positive value','Error')
                    return
                end

                % take the right bound value or format from the input
                handles.RTime = str2num(get(handles.edit10,'String'));
                if isempty(handles.RTime)
                    errordlg('Please insert integer value for ending time to display','Error')
                elseif handles.RTime <= 0
                    errordlg('Ending time to display must be positive value','Error')
                elseif handles.RTime.*1000*handles.ppms > length(stim)
                    errordlg('Right bound of time to display is longer than the existing data','Error')
                end

                % left display time handling
                SortedSptLInd = 0;
                if handles.LTime == 0 % if the user wants to see from the first.
                    Lpos = 1;
                    SortedSptLInd = 1;
                else
                    Lpos = (handles.ppms.*handles.LTime.*1000)+1;
                end
                Rpos = handles.ppms.*handles.RTime.*1000;

                % sorting the stimulus in terms of the time display user specified
                SortedStim = stim(Lpos:Rpos);

                % sorting the spike time data in terms of the time display user specified
                if SortedSptLInd == 0
                    for i = 1:length(spt)
                        if spt(i) > Lpos
                            break
                        end
                    end
                    SortedSptLInd = i;
                end
                if Rpos ~= length(stim)
                    for i = SortedSptLInd:length(spt)
                        if spt(i) > Rpos
                            break
                        end
                    end
                    SortedSptRInd = i-1;
                else
                    SortedSptRInd = length(spt);
                end
                Sortedspt = spt(SortedSptLInd:SortedSptRInd);

                % status bar
                Bar = waitbar(0,'Calculating psth');

                % what models user selected?
                ModelFlag = [get(handles.checkbox5,'Value'),get(handles.checkbox6,'Value'),...
                    get(handles.checkbox7,'Value')];

                if ModelFlag(1) == 1 && ~isfield(handles,'sta_analysis')
                    close(Bar)
                    errordlg('No STA analysis has been done','Error')
                    return
                elseif ModelFlag(2) == 1 && ~isfield(handles,'stc_analysis')
                    close(Bar)
                    errordlg('No STC analysis has been done','Error')
                    return
                elseif ModelFlag(3) == 1 && ~isfield(handles,'glm_analysis')
                    close(Bar)
                    errordlg('No GLM analysis has been done','Error')
                    return
                end

                % psth function
                handles = PSTHgen(SortedStim,Sortedspt,handles,newFFile,Rpos,Lpos,Bar,ModelFlag);

                % updating handles object
                guidata(hObject,handles);
            end
        case 'No'
            return
        case 'Cancle'
            return
    end
end
