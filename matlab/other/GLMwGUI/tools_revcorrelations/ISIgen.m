function handles = ISIgen(Stim,spt,newfi,NumbBin,handles,ModelStatBar,flag)

global sta_scale_factor stc_scale_factor glm_scale_factor ppms

sta_scale_factor = [];
stc_scale_factor = [];
glm_scale_factor = [];
ppms = handles.ppms;

if sum(flag) == 1
    handles = ISIgen_one(Stim,spt,newfi,NumbBin,handles,ModelStatBar,flag);
elseif sum(flag) == 2
    handles = ISIgen_two(Stim,spt,newfi,NumbBin,handles,ModelStatBar,flag);
elseif sum(flag) == 3
    handles = ISIgen_three(Stim,spt,newfi,NumbBin,handles,ModelStatBar,flag);
else
    close(ModelStatBar)
    error('Check at least one model to use.');
end
end

function handles = ISIgen_one(Stim,spt,newfi,NumbBin,handles,ModelStatBar,flag)

if flag(1) == 1
    model = 'STA';
elseif flag(2) == 1
    model = 'STC';
elseif flag(3) == 1
    model = 'GLM';
end

SimISIs = fliplr(sort(diff(spt),'descend'))./handles.ppms;
SimISIsdom = linspace(0,max(SimISIs),NumbBin);
SimPISIs = histc(SimISIs,SimISIsdom);
SimPISIs = SimPISIs./norm(SimPISIs); % normalization

waitbar(0.5,ModelStatBar,['Applying the stimulus to ' model ' model']);

if flag(1) == 1
    Modelspt = STAModel(Stim,handles.sta_analysis,handles.FR);
elseif flag(2) == 1
    Modelspt = STCModel(Stim,handles.stc_analysis,handles.FR);
elseif flag(3) == 1
    Modelspt = GLMModel(Stim,handles.glm_analysis,handles.FR);
end
        
waitbar(0.7,ModelStatBar,'Calculating ISI distribution');
        
ModelISIs = fliplr(sort(diff(Modelspt),'descend'))./handles.ppms;

% default absolute refractory period for GLM
if flag(3) == 1
    GLMISIs = [];
    for i = 1:length(ModelISIs)
        if ModelISIs(i) > 20/handles.ppms
            GLMISIs = [GLMISIs,ModelISIs(i)];
        end
    end
    ModelISIs = GLMISIs;
end

ModelISIsdom = linspace(0,max(ModelISIs),NumbBin);
ModelPISIs = histc(ModelISIs,ModelISIsdom);
ModelPISIs = ModelPISIs./norm(ModelPISIs); % normalization

waitbar(1,ModelStatBar,'Calculating ISI distribution');
close(ModelStatBar)
        
if flag(1) == 1
    handles.ISIpredict.STAPISI = ModelPISIs;
    handles.ISIpredict.STCPISI = [];
    handles.ISIpredict.GLMPISI = [];
elseif flag(2) == 1
    handles.ISIpredict.STAPISI = [];
    handles.ISIpredict.STCPISI = ModelPISIs;
    handles.ISIpredict.GLMPISI = [];
elseif flag(3) == 1
    handles.ISIpredict.STAPISI = [];
    handles.ISIpredict.STCPISI = [];
    handles.ISIpredict.GLMPISI = ModelPISIs;
end
        
figure(7)
semilogy(SimISIsdom,SimPISIs,'LineWidth',4,'Color',[0.5 0.5 0.5])
hold on
if flag(1) == 1
    semilogy(ModelISIsdom,ModelPISIs,'LineWidth',4,'Color','b')
elseif flag(2) == 1
    semilogy(ModelISIsdom,ModelPISIs,'LineWidth',4,'Color','r')
elseif flag(3) == 1
    semilogy(SimISIsdom,ModelPISIs,'LineWidth',4,'Color','g')
end
hold off
set(gca,'FontSize',13)
title([newfi ' ISI prediction ' model ' ' datestr(now)])
xlabel('\Delta t | spike (ms)')
ylabel('log P( \Delta t | spike)')
legend('Target',[model ' Model'])
        
AutoSave = get(handles.checkbox4,'Value');
if AutoSave
    midind = round(length(newfi)/2);
    firstfi = newfi(1:midind-3);
    secondfi = newfi(midind+3:end);
    saveas(figure(7),[firstfi '_' secondfi '_ISI_prediction_' model '_' datestr(now)],'fig');
end
end

function handles = ISIgen_two(Stim,spt,newfi,NumbBin,handles,ModelStatBar,flag)

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

SimISIs = fliplr(sort(diff(spt),'descend'))./handles.ppms;
SimISIsdom = linspace(0,max(SimISIs),NumbBin);
SimPISIs = histc(SimISIs,SimISIsdom);
SimPISIs = SimPISIs./norm(SimPISIs); % normalization

waitbar(0.3,ModelStatBar,['Applying the stimulus to ' modelone ' and ' modeltwo ' model']);

if flag(3) == 0
    Modelonespt = STAModel(Stim,handles.sta_analysis,handles.FR);
    waitbar(0.5,ModelStatBar,['Applying the stimulus to ' modelone ' and ' modeltwo ' model']);
    Modeltwospt = STCModel(Stim,handles.stc_analysis,handles.FR);
elseif flag(2) == 0
    Modelonespt = STAModel(Stim,handles.sta_analysis,handles.FR);
    waitbar(0.5,ModelStatBar,['Applying the stimulus to ' modelone ' and ' modeltwo ' model']);
    Modeltwospt = GLMModel(Stim,handles.glm_analysis,handles.FR);
elseif flag(1) == 0
    Modelonespt = STCModel(Stim,handles.stc_analysis,handles.FR);
    waitbar(0.5,ModelStatBar,['Applying the stimulus to ' modelone ' and ' modeltwo ' model']);
    Modeltwospt = GLMModel(Stim,handles.glm_analysis,handles.FR);
end
        
waitbar(0.7,ModelStatBar,'Calculating ISI distribution');
        
ModeloneISIs = fliplr(sort(diff(Modelonespt),'descend'))./handles.ppms;

ModeloneISIsdom = linspace(0,max(ModeloneISIs),NumbBin);
ModelonePISIs = histc(ModeloneISIs,ModeloneISIsdom);
ModelonePISIs = ModelonePISIs./norm(ModelonePISIs); % normalization

ModeltwoISIs = fliplr(sort(diff(Modeltwospt),'descend'))./handles.ppms;

% default absolute refractory period for GLM
if flag(2) == 0 || flag(1) == 0
    GLMISIs = [];
    for i = 1:length(ModeltwoISIs)
        if ModeltwoISIs(i) > 20/handles.ppms
            GLMISIs = [GLMISIs,ModeltwoISIs(i)];
        end
    end
    ModeltwoISIs = GLMISIs;
end

ModeltwoISIsdom = linspace(0,max(ModeltwoISIs),NumbBin);
ModeltwoPISIs = histc(ModeltwoISIs,ModeltwoISIsdom);
ModeltwoPISIs = ModeltwoPISIs./norm(ModeltwoPISIs); % normalization

waitbar(1,ModelStatBar,'Calculating ISI distribution');
close(ModelStatBar)
        
if flag(3) == 0
    handles.ISIpredict.STAPISI = ModelonePISIs;
    handles.ISIpredict.STCPISI = ModeltwoPISIs;
    handles.ISIpredict.GLMPISI = [];
elseif flag(2) == 0
    handles.ISIpredict.STAPISI = ModelonePISIs;
    handles.ISIpredict.STCPISI = [];
    handles.ISIpredict.GLMPISI = ModeltwoPISIs;
elseif flag(1) == 0
    handles.ISIpredict.STAPISI = [];
    handles.ISIpredict.STCPISI = ModelonePISIs;
    handles.ISIpredict.GLMPISI = ModeltwoPISIs;
end
        
figure(7)
semilogy(SimISIsdom,SimPISIs,'LineWidth',4,'Color',[0.5 0.5 0.5])
hold on
if flag(3) == 0
    semilogy(ModeloneISIsdom,ModelonePISIs,'LineWidth',4,'Color','b')
    hold on
    semilogy(ModeltwoISIsdom,ModeltwoPISIs,'LineWidth',4,'Color','r')
elseif flag(2) == 0
    semilogy(ModeloneISIsdom,ModelonePISIs,'LineWidth',4,'Color','b')
    hold on
    semilogy(ModeltwoISIsdom,ModeltwoPISIs,'LineWidth',4,'Color','g')
elseif flag(1) == 0
    semilogy(ModeloneISIsdom,ModelonePISIs,'LineWidth',4,'Color','r')
    hold on
    semilogy(ModeltwoISIsdom,ModeltwoPISIs,'LineWidth',4,'Color','g')
end
hold off
set(gca,'FontSize',13)
title([newfi ' ISI prediction ' modelone ' and ' modeltwo ' ' datestr(now)])
xlabel('\Delta t | spike (ms)')
ylabel('log P( \Delta t | spike)')
legend('Target',[modelone ' Model'],[modeltwo ' Model'])
        
AutoSave = get(handles.checkbox4,'Value');
if AutoSave
    midind = round(length(newfi)/2);
    firstfi = newfi(1:midind-3);
    secondfi = newfi(midind+3:end);
    saveas(figure(7),[firstfi '_' secondfi '_ISI_prediction_' modelone '_' modeltwo '_' datestr(now)],'fig');
end
end

function handles = ISIgen_three(Stim,spt,newfi,NumbBin,handles,ModelStatBar,flag)

SimISIs = fliplr(sort(diff(spt),'descend'))./handles.ppms;
SimISIsdom = linspace(0,max(SimISIs),NumbBin);
SimPISIs = histc(SimISIs,SimISIsdom);
SimPISIs = SimPISIs./norm(SimPISIs); % normalization

waitbar(0.3,ModelStatBar,'Applying the stimulus to STA, STC, and GLM model');

STAspt = STAModel(Stim,handles.sta_analysis,handles.FR);

waitbar(0.4,ModelStatBar,'Applying the stimulus to STA, STC, and GLM model');

STCspt = STCModel(Stim,handles.stc_analysis,handles.FR);

waitbar(0.5,ModelStatBar,'Applying the stimulus to STA, STC, and GLM model');

GLMspt = GLMModel(Stim,handles.glm_analysis,handles.FR);

waitbar(0.6,ModelStatBar,'Calculating ISI distribution');
        
STAISIs = fliplr(sort(diff(STAspt),'descend'))./handles.ppms;
STCISIs = fliplr(sort(diff(STCspt),'descend'))./handles.ppms;
GLMISIs = fliplr(sort(diff(GLMspt),'descend'))./handles.ppms;

% default absolute refractory period for GLM
tempGLMISIs = [];
for i = 1:length(GLMISIs)
    if GLMISIs(i) > 20/handles.ppms
        tempGLMISIs = [tempGLMISIs,GLMISIs(i)];
    end
end
GLMISIs = tempGLMISIs;

STAISIsdom = linspace(0,max(STAISIs),NumbBin);
STAPISIs = histc(STAISIs,STAISIsdom);
STAPISIs = STAPISIs./norm(STAPISIs); % normalization

STCISIsdom = linspace(0,max(STCISIs),NumbBin);
STCPISIs = histc(STCISIs,STCISIsdom);
STCPISIs = STCPISIs./norm(STCPISIs); % normalization

GLMISIsdom = linspace(0,max(GLMISIs),NumbBin);
GLMPISIs = histc(GLMISIs,GLMISIsdom);
GLMPISIs = GLMPISIs./norm(GLMPISIs); % normalization

waitbar(1,ModelStatBar,'Calculating ISI distribution');
close(ModelStatBar)
        
handles.ISIpredict.STAPISI = STAPISIs;
handles.ISIpredict.STCPISI = STCPISIs;
handles.ISIpredict.GLMPISI = GLMPISIs;
        
figure(7)
semilogy(SimISIsdom,SimPISIs,'LineWidth',4,'Color',[0.5 0.5 0.5])
hold on
semilogy(STAISIsdom,STAPISIs,'LineWidth',4,'Color','b')
hold on
semilogy(STCISIsdom,STCPISIs,'LineWidth',4,'Color','r')
hold on
semilogy(GLMISIsdom,GLMPISIs,'LineWidth',4,'Color','g')
hold off
set(gca,'FontSize',13)
title([newfi ' ISI prediction STA, STC, and GLM ' datestr(now)])
xlabel('\Delta t | spike (ms)')
ylabel('log P( \Delta t | spike)')
legend('Target','STA','STC','GLM')
        
AutoSave = get(handles.checkbox4,'Value');
if AutoSave
    midind = round(length(newfi)/2);
    firstfi = newfi(1:midind-3);
    secondfi = newfi(midind+3:end);
    saveas(figure(7),[firstfi '_' secondfi '_ISI_prediction_STA_STC_GLM_' datestr(now)],'fig');
end
end