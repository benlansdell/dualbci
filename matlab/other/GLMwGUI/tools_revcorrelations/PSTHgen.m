function handles = PSTHgen(Stim,spt,handles,newfi,Rpos,Lpos,ModelStatBar,flag)

global sta_scale_factor stc_scale_factor glm_scale_factor ppms

sta_scale_factor = [];
stc_scale_factor = [];
glm_scale_factor = [];
ppms = handles.ppms;

GaussScale = 0.007;
if sum(flag) == 1
    handles = PSTHgen_one(Stim,spt,handles,newfi,Rpos,Lpos,ModelStatBar,flag,GaussScale);
elseif sum(flag) == 2
    handles = PSTHgen_two(Stim,spt,handles,newfi,Rpos,Lpos,ModelStatBar,flag,GaussScale);
elseif sum(flag) == 3
    handles = PSTHgen_three(Stim,spt,handles,newfi,Rpos,Lpos,ModelStatBar,flag,GaussScale);
else
    close(ModelStatBar)
    errordlg('Check at least one model to use.','Error');
    return
end
end

function handles = PSTHgen_one(Stim,spt,handles,newfi,Rpos,Lpos,ModelStatBar,flag,GaussScale)

% sliding gaussian window setting.
GaWid = GaussScale.*(Rpos - Lpos + 1); % the std of the Gaussian filter window
GaDom = 3*GaWid; % the width of the Gaussian filter window
Ga = exp(-0.5*((-GaDom:GaDom)./GaWid).^2)./(sqrt(2*pi).*GaWid); % Gaussian window

Rep = handles.NumbRep;

if flag(1) == 1
    model = 'STA';
elseif flag(2) == 1
    model = 'STC';
elseif flag(3) == 1
    model = 'GLM';
end

for i = 1:Rep
    waitbar(i/Rep,ModelStatBar,['Applying the stimulus to ' model ...
        ' and generating spike time']);
    if flag(1) == 1
        Modelspt{i} = STAModel(Stim,handles.sta_analysis,handles.FR);
    elseif flag(2) == 1
        Modelspt{i} = STCModel(Stim,handles.stc_analysis,handles.FR);
    elseif flag(3) == 1
        Modelspt{i} = GLMModel(Stim,handles.glm_analysis,handles.FR);
    end
end
        
close(ModelStatBar)
        
figure(6)
set(figure(6),'Position',[.5 .5 1000 500])
for ii = 1:length(spt)
            
    % simulation line..displayed by black line
    subplot(2,1,1)
    line([spt(ii)/(1000*handles.ppms) spt(ii)/(1000*handles.ppms)],[0 length(Modelspt)+1],'Color',[.5 .5 .5]);
    set(gca,'FontSize',13)
    xlabel('Time')
    ylabel('Trial')
    if flag(1) == 1
        title([newfi ' PSTH analysis (Grey : simulation, Blue : STA) ' ' ' datestr(now)])
    elseif flag(2) == 1
        title([newfi ' PSTH analysis (Grey : simulation, Red : STC) ' ' ' datestr(now)])
    elseif flag(3) == 1
        title([newfi ' PSTH analysis (Grey : simulation, Green : GLM) ' ' ' datestr(now)])
    end
    axis([handles.LTime handles.RTime 0 Rep+1])
end
        
% Model raster displayed by blue line
for i = 1:length(Modelspt)
    tempModelspt = Modelspt{i};
    for ii = 1:length(tempModelspt)
        if Lpos == 0
            Xval = tempModelspt(ii);
        else
            Xval = tempModelspt(ii) + Lpos - 1;
        end
        subplot(2,1,1)
        if flag(1) == 1
            line([Xval/(1000*handles.ppms) Xval/(1000*handles.ppms)],[i-1 i],'Color','b');
        elseif flag(2) == 1
            line([Xval/(1000*handles.ppms) Xval/(1000*handles.ppms)],[i-1 i],'Color','r');
        elseif flag(3) == 1
            line([Xval/(1000*handles.ppms) Xval/(1000*handles.ppms)],[i-1 i],'Color','g');
        end
    end
    hold on
    if i == length(Modelspt)
        psth = SpkTimeCell2Psth(Modelspt,length(Stim));
        psth = psth(1:(Rpos-Lpos+1));
        Modelpsth = conv2(psth,Ga);
        Modelpsth = Modelpsth(GaDom+1:length(Modelpsth)-GaDom);
        Temppsth = Modelpsth;
        Modelpsth = Modelpsth./max(Modelpsth);
                
        % Model PSTH display
        subplot(2,1,2)
        if flag(1) == 1
            area(linspace(handles.LTime,handles.RTime,length(Modelpsth)),...
            Modelpsth,'LineStyle',':','FaceColor','b')
        elseif flag(2) == 1
            area(linspace(handles.LTime,handles.RTime,length(Modelpsth)),...
            Modelpsth,'LineStyle',':','FaceColor','r')
        elseif flag(3) == 1
            area(linspace(handles.LTime,handles.RTime,length(Modelpsth)),...
            Modelpsth,'LineStyle',':','FaceColor','g')
        end
        hold on
        for ii = 1:length(spt)
            line([spt(ii)/(1000*handles.ppms) spt(ii)/(1000*handles.ppms)],[0 1],'Color',[.5 .5 .5]);
        end
        set(gca,'FontSize',13)
        xlabel('Time (sec)')
        axis tight
        hold on
    end
end
hold off
AutoSave = get(handles.checkbox4,'Value');
if AutoSave
    midind = round(length(newfi)/2);
    firstfi = newfi(1:midind-3);
    secondfi = newfi(midind+3:end);
    saveas(figure(6),[firstfi '_' secondfi '_PSTH_prediction_' model '_' datestr(now)],'fig');
end

if flag(1) == 1
     handles.PSTHpredict.STApsth = Temppsth;
     handles.PSTHpredict.STCpsth = [];
     handles.PSTHpredict.GLMpsth = [];
elseif flag(2) == 1
     handles.PSTHpredict.STApsth = [];
     handles.PSTHpredict.STCpsth = Temppsth;
     handles.PSTHpredict.GLMpsth = [];
elseif flag(3) == 1
     handles.PSTHpredict.STApsth = [];
     handles.PSTHpredict.STCpsth = [];
     handles.PSTHpredict.GLMpsth = Temppsth;
end

end

function handles = PSTHgen_two(Stim,spt,handles,newfi,Rpos,Lpos,ModelStatBar,flag,GaussScale)
% sliding gaussian window setting.
GaWid = GaussScale.*(Rpos - Lpos + 1); % the std of the Gaussian filter window
GaDom = 3*GaWid; % the width of the Gaussian filter window
Ga = exp(-0.5*((-GaDom:GaDom)./GaWid).^2)./(sqrt(2*pi).*GaWid); % Gaussian window

Rep = handles.NumbRep;

if flag(3) == 0
    modelone = 'STA';
    modeltwo = 'STC';
elseif flag(2) == 0
    modelone = 'STA';
    modeltwo = 'GLM';
elseif flag(1) == 0
    modelone = 'STC';
    modeltwo = 'GLM';
end

for i = 1:Rep
    waitbar(i/Rep,ModelStatBar,['Applying the stimulus to ' modelone ' model and ' modeltwo ' and generating spike time']);
    if flag(3) == 0
        Modelonespt{i} = STAModel(Stim,handles.sta_analysis,handles.FR);
        Modeltwospt{i} = STCModel(Stim,handles.stc_analysis,handles.FR);
    elseif flag(2) == 0
        Modelonespt{i} = STAModel(Stim,handles.sta_analysis,handles.FR);
        Modeltwospt{i} = GLMModel(Stim,handles.glm_analysis,handles.FR);
    elseif flag(1) == 0
        Modelonespt{i} = STCModel(Stim,handles.stc_analysis,handles.FR);
        Modeltwospt{i} = GLMModel(Stim,handles.glm_analysis,handles.FR);
    end
end
        
close(ModelStatBar)
        
figure(6)
set(figure(6),'Position',[.5 .5 1300 500])
if flag(3) == 0
    annotation(figure(6),'textbox',...
        [0.24045044160942 0.959899749373433 0.557390578999013 0.0325814536340841],...
        'String',{[newfi ' PSTH analysis (Grey : simulation, Blue : STA, Red : STC) ' ' ' datestr(now)]},...
        'HorizontalAlignment','center',...
        'FontSize',13,...
        'FitBoxToText','off',...
        'LineStyle','none');
elseif flag(2) == 0
    annotation(figure(6),'textbox',...
        [0.24045044160942 0.959899749373433 0.557390578999013 0.0325814536340841],...
        'String',{[newfi ' PSTH analysis (Grey : simulation, Blue : STA, Green : GLM) ' ' ' datestr(now)]},...
        'HorizontalAlignment','center',...
        'FontSize',13,...
        'FitBoxToText','off',...
        'LineStyle','none');
elseif flag(1) == 0
    annotation(figure(6),'textbox',...
        [0.24045044160942 0.959899749373433 0.557390578999013 0.0325814536340841],...
        'String',{[newfi ' PSTH analysis (Grey : simulation, Red : STC, Green : GLM) ' ' ' datestr(now)]},...
        'HorizontalAlignment','center',...
        'FontSize',13,...
        'FitBoxToText','off',...
        'LineStyle','none');
end

for ii = 1:length(spt)
            
    % simulation line..displayed by black line
    subplot(2,2,1)
    line([spt(ii)/(1000*handles.ppms) spt(ii)/(1000*handles.ppms)],[0 length(Modelonespt)+1],'Color',[.5 .5 .5]);
    set(gca,'FontSize',13)
    xlabel('Time')
    ylabel('Trial')
    axis([handles.LTime handles.RTime 0 Rep+1])
    hold on
            
    subplot(2,2,3)
    line([spt(ii)/(1000*handles.ppms) spt(ii)/(1000*handles.ppms)],[0 length(Modeltwospt)+1],'Color',[.5 .5 .5]);
    set(gca,'FontSize',13)
    xlabel('Time')
    ylabel('Trial')
    axis([handles.LTime handles.RTime 0 Rep+1])
    hold on
end
        
% First model raster displayed by blue line
for i = 1:length(Modelonespt)
	tempFirstmodelspt = Modelonespt{i};
    for ii = 1:length(tempFirstmodelspt)
        if Lpos == 0
            Xval = tempFirstmodelspt(ii);
        else
            Xval = tempFirstmodelspt(ii) + Lpos - 1;
        end
        subplot(2,2,1)
        if flag(3) == 0
            line([Xval/(1000*handles.ppms) Xval/(1000*handles.ppms)],[i-1 i],'Color','b');
        elseif flag(2) == 0
            line([Xval/(1000*handles.ppms) Xval/(1000*handles.ppms)],[i-1 i],'Color','b');
        elseif flag(1) == 0
            line([Xval/(1000*handles.ppms) Xval/(1000*handles.ppms)],[i-1 i],'Color','r');
        end
    end
    hold on
    if i == length(Modelonespt)
        psth = SpkTimeCell2Psth(Modelonespt,length(Stim));
        psth = psth(1:(Rpos-Lpos+1));
        Firstmodelpsth = conv2(psth,Ga);
        Firstmodelpsth = Firstmodelpsth(GaDom+1:length(Firstmodelpsth)-GaDom);
                
        % First model PSTH display
        subplot(2,2,4)
        for ii = 1:length(spt)
            line([spt(ii)/(1000*handles.ppms) spt(ii)/(1000*handles.ppms)],[-1 1],'Color',[.5 .5 .5]);
        end
        set(gca,'FontSize',13)
        xlabel('Time (sec)')
        axis tight
        hold on
    end
end
        
% Second model raster displayed by blue line
for i = 1:length(Modeltwospt)
	tempSecondmodelspt = Modeltwospt{i};
    for ii = 1:length(tempSecondmodelspt)
        if Lpos == 0
            Xval = tempSecondmodelspt(ii);
        else
            Xval = tempSecondmodelspt(ii) + Lpos - 1;
        end
        subplot(2,2,3)
        if flag(3) == 0
            line([Xval/(1000*handles.ppms) Xval/(1000*handles.ppms)],[i-1 i],'Color','r');
        elseif flag(2) == 0
            line([Xval/(1000*handles.ppms) Xval/(1000*handles.ppms)],[i-1 i],'Color','g');
        elseif flag(1) == 0
            line([Xval/(1000*handles.ppms) Xval/(1000*handles.ppms)],[i-1 i],'Color','g');
        end
    end
    hold on
    if i == length(Modeltwospt)
        psth = SpkTimeCell2Psth(Modeltwospt,length(Stim));
        psth = psth(1:(Rpos-Lpos+1));
        Secondmodelpsth = conv2(psth,Ga);
        Secondmodelpsth = Secondmodelpsth(GaDom+1:length(Secondmodelpsth)-GaDom);
    end
end     

TempFirstPsth = Firstmodelpsth;
TempSecondPsth = Secondmodelpsth;

Firstpsthmax = max(Firstmodelpsth);
Secondpsthmax = max(Secondmodelpsth);

totalmax = max([Firstpsthmax,Secondpsthmax]);

Firstmodelpsth = Firstmodelpsth./totalmax;
Secondmodelpsth = Secondmodelpsth./totalmax;

Secondmodelpsth = -Secondmodelpsth;

if flag(3) == 0
    subplot(2,2,4)
    area(linspace(handles.LTime,handles.RTime,length(Firstmodelpsth)),...
        Firstmodelpsth,'LineStyle',':','FaceColor','b')
    hold on
    area(linspace(handles.LTime,handles.RTime,length(Secondmodelpsth)),...
        Secondmodelpsth,'LineStyle',':','FaceColor','r')
elseif flag(2) == 0
    subplot(2,2,4)
    area(linspace(handles.LTime,handles.RTime,length(Firstmodelpsth)),...
        Firstmodelpsth,'LineStyle',':','FaceColor','b')
    hold on
    area(linspace(handles.LTime,handles.RTime,length(Secondmodelpsth)),...
        Secondmodelpsth,'LineStyle',':','FaceColor','g')
elseif flag(1) == 0
    subplot(2,2,4)
    area(linspace(handles.LTime,handles.RTime,length(Firstmodelpsth)),...
        Firstmodelpsth,'LineStyle',':','FaceColor','r')
    hold on
    area(linspace(handles.LTime,handles.RTime,length(Secondmodelpsth)),...
        Secondmodelpsth,'LineStyle',':','FaceColor','g')
end
hold off
AutoSave = get(handles.checkbox4,'Value');
if AutoSave
    midind = round(length(newfi)/2);
    firstfi = newfi(1:midind-3);
    secondfi = newfi(midind+3:end);
    saveas(figure(6),[firstfi '_' secondfi '_PSTH_prediction_' modelone '_' modeltwo '_' datestr(now)],'fig');
end

if flag(3) == 0
    handles.PSTHpredict.STApsth = TempFirstPsth;
    handles.PSTHpredict.STCpsth = TempSecondPsth;
    handles.PSTHpredict.GLMpsth = [];
elseif flag(2) == 0
    handles.PSTHpredict.STApsth = TempFirstPsth;
    handles.PSTHpredict.STCpsth = [];
    handles.PSTHpredict.GLMpsth = TempSecondPsth;
elseif flag(1) == 0
    handles.PSTHpredict.STApsth = [];
    handles.PSTHpredict.STCpsth = TempFirstPsth;
    handles.PSTHpredict.GLMpsth = TempSecondPsth;
end

end

function handles = PSTHgen_three(Stim,spt,handles,newfi,Rpos,Lpos,ModelStatBar,flag,GaussScale)

% sliding gaussian window setting.
GaWid = GaussScale.*(Rpos - Lpos + 1); % the std of the Gaussian filter window
GaDom = 3*GaWid; % the width of the Gaussian filter window
Ga = exp(-0.5*((-GaDom:GaDom)./GaWid).^2)./(sqrt(2*pi).*GaWid); % Gaussian window

Rep = handles.NumbRep;

for i = 1:Rep
    waitbar(i/Rep,ModelStatBar,'Applying the stimulus to STA,STC, and GLM model and generating spike time');
    Modelonespt{i} = STAModel(Stim,handles.sta_analysis,handles.FR);
    Modeltwospt{i} = STCModel(Stim,handles.stc_analysis,handles.FR);
    Modelthreespt{i} = GLMModel(Stim,handles.glm_analysis,handles.FR);
end
        
close(ModelStatBar)
        
figure(6)
set(figure(6),'Position',[.5 .5 1300 500])
annotation(figure(6),'textbox',...
    [0.204111647920283 0.950705882352941 0.63396489771269 0.0420000000000001],...
    'String',{[newfi ' PSTH analysis (Grey : simulation, Blue : STA, Red : STC, Green : GLM) ' ' ' datestr(now)]},...
    'HorizontalAlignment','center',...
    'FontSize',13,...
    'FitBoxToText','off',...
    'LineStyle','none');
for ii = 1:length(spt)
    % simulation line..displayed by black line
    subplot(3,2,1)
    line([spt(ii)/(1000*handles.ppms) spt(ii)/(1000*handles.ppms)],[0 length(Modelonespt)+1],'Color',[.5 .5 .5]);
    set(gca,'FontSize',13)
    xlabel('Time')
    ylabel('Trial')
    axis([handles.LTime handles.RTime 0 Rep+1])
    hold on
            
    subplot(3,2,3)
    line([spt(ii)/(1000*handles.ppms) spt(ii)/(1000*handles.ppms)],[0 length(Modeltwospt)+1],'Color',[.5 .5 .5]);
    set(gca,'FontSize',13)
    xlabel('Time')
    ylabel('Trial')
    axis([handles.LTime handles.RTime 0 Rep+1])
    hold on
    
    subplot(3,2,5)
    line([spt(ii)/(1000*handles.ppms) spt(ii)/(1000*handles.ppms)],[0 length(Modelthreespt)+1],'Color',[.5 .5 .5]);
    set(gca,'FontSize',13)
    xlabel('Time')
    ylabel('Trial')
    axis([handles.LTime handles.RTime 0 Rep+1])
    hold on
end
        
% First model raster displayed by blue line
for i = 1:length(Modelonespt)
	tempFirstmodelspt = Modelonespt{i};
    for ii = 1:length(tempFirstmodelspt)
        if Lpos == 0
            Xval = tempFirstmodelspt(ii);
        else
            Xval = tempFirstmodelspt(ii) + Lpos - 1;
        end
        subplot(3,2,1)
        line([Xval/(1000*handles.ppms) Xval/(1000*handles.ppms)],[i-1 i],'Color','b');
    end
    hold on
    if i == length(Modelonespt)
        psth = SpkTimeCell2Psth(Modelonespt,length(Stim));
        psth = psth(1:(Rpos-Lpos+1));
        Firstmodelpsth = conv2(psth,Ga);
        Firstmodelpsth = Firstmodelpsth(GaDom+1:length(Firstmodelpsth)-GaDom);
                
        % First model PSTH display
        subplot(3,2,[4,6])
        for ii = 1:length(spt)
            line([spt(ii)/(1000*handles.ppms) spt(ii)/(1000*handles.ppms)],[-2 1],'Color',[.5 .5 .5]);
        end
        set(gca,'FontSize',13)
        xlabel('Time (sec)')
        axis tight
        hold on
    end
end
        
% Second model raster displayed by blue line
for i = 1:length(Modeltwospt)
	tempSecondmodelspt = Modeltwospt{i};
    for ii = 1:length(tempSecondmodelspt)
        if Lpos == 0
            Xval = tempSecondmodelspt(ii);
        else
            Xval = tempSecondmodelspt(ii) + Lpos - 1;
        end
        subplot(3,2,3)
        line([Xval/(1000*handles.ppms) Xval/(1000*handles.ppms)],[i-1 i],'Color','r');
    end
    hold on
    if i == length(Modeltwospt)
        psth = SpkTimeCell2Psth(Modeltwospt,length(Stim));
        psth = psth(1:(Rpos-Lpos+1));
        Secondmodelpsth = conv2(psth,Ga);
        Secondmodelpsth = Secondmodelpsth(GaDom+1:length(Secondmodelpsth)-GaDom);
    end
end

% Third model raster displayed by blue line
for i = 1:length(Modelthreespt)
	tempThirdmodelspt = Modelthreespt{i};
    for ii = 1:length(tempThirdmodelspt)
        if Lpos == 0
            Xval = tempThirdmodelspt(ii);
        else
            Xval = tempThirdmodelspt(ii) + Lpos - 1;
        end
        subplot(3,2,5)
        line([Xval/(1000*handles.ppms) Xval/(1000*handles.ppms)],[i-1 i],'Color','g');
    end
    hold on
    if i == length(Modelthreespt)
        psth = SpkTimeCell2Psth(Modelthreespt,length(Stim));
        psth = psth(1:(Rpos-Lpos+1));
        Thirdmodelpsth = conv2(psth,Ga);
        Thirdmodelpsth = Thirdmodelpsth(GaDom+1:length(Thirdmodelpsth)-GaDom);
    end
end

TempFirstPsth = Firstmodelpsth;
TempSecondPsth = Secondmodelpsth;
TempThirdPsth = Thirdmodelpsth;

Firstmodelpsth = Firstmodelpsth./max(Firstmodelpsth);
Secondmodelpsth = Secondmodelpsth./max(Secondmodelpsth);
Thirdmodelpsth = Thirdmodelpsth./max(Thirdmodelpsth);

Firstpsthmax = max(Firstmodelpsth);
Secondpsthmax = max(Secondmodelpsth);
Thirdpsthmax = max(Thirdmodelpsth);

totalmax = max([Firstpsthmax,Secondpsthmax,Thirdpsthmax]);

Firstmodelpsth = Firstmodelpsth./totalmax;
Secondmodelpsth = Secondmodelpsth./totalmax;
Thirdmodelpsth = Thirdmodelpsth./totalmax;

Secondmodelpsth = -Secondmodelpsth;
Thirdmodelpsth = Thirdmodelpsth - 2;
basevalue = -2;

subplot(3,2,[4,6])
area(linspace(handles.LTime,handles.RTime,length(Firstmodelpsth)),...
    Firstmodelpsth,'LineStyle',':','FaceColor','b')
hold on
subplot(3,2,[4,6])
area(linspace(handles.LTime,handles.RTime,length(Secondmodelpsth)),...
    Secondmodelpsth,'LineStyle',':','FaceColor','r')  
hold on
subplot(3,2,[4,6])
area(linspace(handles.LTime,handles.RTime,length(Thirdmodelpsth)),...
    Thirdmodelpsth,basevalue,'LineStyle',':','FaceColor','g')
hold off

AutoSave = get(handles.checkbox4,'Value');
if AutoSave
    midind = round(length(newfi)/2);
    firstfi = newfi(1:midind-3);
    secondfi = newfi(midind+3:end);
    saveas(figure(6),[firstfi '_' secondfi '_PSTH_prediction_STA_STC_GLM_' datestr(now)],'fig');
end

Secondmodelpsth = Secondmodelpsth./max(Secondmodelpsth);

handles.PSTHpredict.STApsth = TempFirstPsth;
handles.PSTHpredict.STCpsth = TempSecondPsth;
handles.PSTHpredict.GLMpsth = TempThirdPsth;

end

function psth = SpkTimeCell2Psth(spk,len)
% This function takes the spike time data and the length of the spike and
% produce the psth.

psth = 0;

for i = 1:length(spk)
    spike_time = spk{i};
    spt(spike_time) = 1;
    spt(end:len) = 0;
    psth = psth + spt;
end
end