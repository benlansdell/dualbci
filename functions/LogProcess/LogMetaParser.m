%
%
% State machine for finding experiment intervals of interest
%
% 1) Find run/pause section
%  
% 2) Find last parameter update during this segment
%
%   
%
% 3) Find Run/Pause -> Run time immediately following time from 1)
% 4) Find Run/Pause -> Pause following 3)

% Revised seqence:
% 1) same, find push/push or push pull setting
% 2) find following run/pause start & run/pause stop
% 3) Verify no changes to pvmap in between

% Revised again sequence:
% 1) find *first* push/push or push/pull setting
% 2) find following change *away* from this setting, if any
% 3) ensure non-zero trials in between

function trials = LogMetaParser(data, tstart, tstop, dbFlag)
global matpath
try
if (nargin < 5)
    dbFlag = false;
end
if (nargin == 1)
    tstart = data.messages.time(1);
    tstop = data.messages.time(end);
end

% times of run events
tRun = data.messages.time(...
            data.messages.type == 1 & ...
            data.messages.data(:,1) == 1);
tRun = tRun(tRun >= tstart & tRun <= tstop)';
        
tPause = data.messages.time(...
            data.messages.type == 1 & ...
            data.messages.data(:,1) == 0);
tPause = tPause(tPause >= tstart & tPause <= tstop)';

% iterate over segments

segment = 1;
count = 0;
trials = data.trials;
N = length(data.trials.index);  

if (N <= 50)
    fprintf('Less than 50 trials, ignoring.\n\n')
    trials = [];
    return
end

trials.segment = zeros(N,1);        % index continuous segments w/ same map
trials.Success = trials.success;
trials.Map = cell(N,1);
trials.Map(1:N) = {''};
trials.MoveDir(N) = ' ';
trials.Neurons = zeros(N,4); % will be populated with variable-length arrays
trials.Date = cell(N,1);
trials.Order = zeros(N,1);
trials.FEGain = zeros(N,1);
%trials.velControl = zeros(N,1);
trials.W = zeros(N,1);
trials.D = zeros(N,1);
trials.IDW = zeros(N,1);
trials.IDM = zeros(N,1);
trials.ID = zeros(N,1);
trials.MT = zeros(N,1);
trials.HT = zeros(N,1);
trials.endPos = zeros(N,2);
trials.overshoots = zeros(N,1);
trials.TidxSets = cell(N,1);
trials.reward = zeros(N,1);

try 
    if (length(data.nev) > 0)
        nevFlag = true;
    end
catch me
    nevFlag = false;
end
if nevFlag
NEVtimes = cell2mat({data.nev(:).Toffset}');    % NEV start times +/- err in sample periods
NEVtimes = NEVtimes(:,1)/data.sampleRate;       % just avg start time, in sec
end

segment = 0;
lastMap = '';
neurons = [];    % Will contain neurons active during this trial

for tr = tRun
    count = count + 1;
    % sanity check: skip any segments lasting less than 2 min
    tp = tPause(find(tPause > tr,1,'first'));
    if (isempty(tp))
        tp = data.messages.time(end);
    end
    if ((tp - tr) < 2*60)
        sprintf('Warning: skipping %1.2f min section\n',(tp-tr)/60);
        continue
    end
    
    
    
    %tD = SessionToTrialData(data, tr, tp);
    
    % find parameter settings for segment
    pvmapIdx = find(data.pvmap.time < tp,1,'last');
    paramsIdx = find(data.params.time < tp,1, 'last');
    warnFlag = false;
    if (data.pvmap.time(pvmapIdx) > tr || data.params.time(paramsIdx) > tr)
        sprintf('Warning: settings changed during segment %d',count);
        warnFlag = true;
    end
    
    % figure out meta data: 
    % mapping type (push/push, push/pull, manual 1D, manual 2D)
    % active neurons
    enabled = data.pvmap.enabled(pvmapIdx,:);
    angles = data.pvmap.angles(pvmapIdx,logical(enabled));
    %angles = [angles, zeros(1,size(enabled,2)-size(angles,2))];  (pads with zeros, dumb idea)
    intDisable = data.pvmap.intDisable(pvmapIdx,:);

% mapDirs = [find(enabled); angles];  % enabled channel indices & angles
    
% offXEnabled = (sum((enabled .* angles ~= 180) & (enabled.*angles ~= 0) > 0));
    
    
    % detect dual control, if not dual then much easier from there
    if ( sum(enabled(1:4)) && sum(enabled(5:end)))
        map = 'Dual Control';
    else
        
        % 1D vs. 2D: use trig functions, assume 1D is horiz or vert
        vmap = sum(abs(sin(angles*pi/180))) > 1-eps;
        hmap = sum(abs(cos(angles*pi/180))) > 1-eps;
        
        if (xor(vmap, hmap))  % vertical or horizontal but not both
            if (vmap)
                map = '1D Vert  ';
            else
                map = '1D Horiz ';
            end
        else
            map = '2D ';
        end
        
        % Determine brain/manual
        if (sum(enabled(1:4)))
            map = [map,'Brain  '];
        else
            map = [map,'Manual '];
        end
        
        % Determine velocity/position in controlled direction
        if ( ((vmap) && intDisable(2)) || ((hmap) && intDisable(1)) )
            map = [map,'Position'];
        else 
            map = [map,'Velocity'];
        end
    end

        
        idxSet = find(data.trials.time >= tr & data.trials.time <= tp)';

        %activeNeurons = neurons.*data.pvmap.enabled(pvmapIdx,1:4);
    
    if (dbFlag)
        fprintf('Processing trials %d - %d\n',idxSet(1),idxSet(end))
    end
    
%%% This logic intended to lump sequential segments w/ same map together,
%%% but does not work as designed. Instead let every segment be separate.
    if (~strcmp(map,lastMap))       % first possibly sub-segment of new seg
            segment = segment + 1;
            lastMap = map;
                trials.segInfo(segment).idxSet = idxSet;
    trials.segInfo(segment).trun = tr;
    trials.segInfo(segment).tpause = tp;
    trials.segInfo(segment).map = map;
    trials.segInfo(segment).enabled = enabled;
    trials.segInfo(segment).pvmapIdx = pvmapIdx;
    trials.segInfo(segment).warn = warnFlag;
            
    else % for subsequent sub-segments, update total included trials + end time
        trials.segInfo(segment).idxSet = [trials.segInfo(segment).idxSet, idxSet];
        trials.segInfo(segment).tpause = tp;
        trials.segInfo(segment).warn = trials.segInfo(segment).warn || warnFlag;
    end

    % meta data that can be batch-assigned
    trials.W(idxSet) = 2*data.params.successRadius(paramsIdx);
    trials.HT(idxSet) = data.params.dwellTime(paramsIdx);
    %trials.Neurons(idxSet,:) = activeNeurons;
    trials.Order(idxSet) = count;
    trials.MT(idxSet) = trials.duration(idxSet) - trials.HT(idxSet);
    trials.Map(idxSet') = {map};
    trials.segment(idxSet) = segment;
    if (size(data.pvmap.gains,2) > 4)
        trials.FEGain(idxSet) = data.pvmap.gains(pvmapIdx,5);
    end
    
    % add meta data to each trial 
    for i = idxSet
        % Keep neurons in sync with nev info
        try
        if nevFlag && ~isempty(find(NEVtimes < trials.time(i),1))
            
            neurons = data.nev(find(NEVtimes < trials.time(i),1,'last')).chans;
            
        else
            neurons = -[1 1 1 1];
        end
        
        trials.D(i) = norm(data.trials.targetPos(i,:) - ...
                          data.trials.startPos(i,:),2);
        trials.Neurons(i,:) = neurons.*enabled(1:4);    
        trials.Date{i} = data.fname;
        catch me
            me.getReport
            keyboard
        end
        %trials.ID(i) = log2(trials.A(i)/trials.W(i) + 1);
        %trials.IDW(i) = log2((trials.A(i) - trials.W(i)/2)/trials.W(i) + 1);

        if (trials.targetPos(i,1) - trials.startPos(i,1) > 0)
            trials.MoveDir(i) = 'L';
        else
            trials.MoveDir(i) = 'R';
        end
        
        % take last sampled cursor position before trial end as ending
        % cursor position
        timeIdx = find(data.stateHist.time <= trials.time(i),1,'last');
        trials.endPos(i,:) = data.stateHist.cursor(timeIdx,:);
        
        % Count overshoots: find trajectory segment of trial, then error
        idxRange = find(data.stateHist.time < data.trials.time(i) - data.trials.duration(i),1,'last'):...
                find(data.stateHist.time < data.trials.time(i),1,'last');
        trials.TidxSets{i} = idxRange;
        errorNorm = sqrt(data.stateHist.error(idxRange,1).^2 + data.stateHist.error(idxRange,2).^2);
        trials.overshoots(i) = floor(sum(abs(diff(sign(errorNorm-trials.W(i)/2)/2)))/2);
        if (dbFlag > 1)
            figure(20)
            clf
            plot(errorNorm, 'b-')
            hold on
            plot([1 length(errorNorm)],[1 1]*trials.W(i)/2,'k:');
            fprintf('Overshoots: %d\n',trials.overshoots(i))
            pause
        end
            
        % determine whether reward was given, if reward event happened
        % within 0.1s of the trial end
        [tOff,ir] = min(abs(data.rewards.time - trials.time(i)));
        trials.reward(i) = tOff < 0.1 && data.rewards.given(ir);
    end
    
    % more batch assignment, dependent on A
    
    % Standard Shannon ID
    trials.ID(idxSet) = log2(trials.D(idxSet)./trials.W(idxSet) + ones(length(idxSet),1));    
    
    % Shenoy ID, compensated for width
    trials.IDW(idxSet) = log2((trials.D(idxSet) - trials.W(idxSet)/2)./trials.W(idxSet) + ones(length(idxSet),1));
    
    % Meyer et al 1988 definition, compensated for target width
    trials.IDM(idxSet) = ((trials.D(idxSet) - trials.W(idxSet)/2)./trials.W(idxSet)).^(0.5);
    
end

% FUCK YOU MATLAB (undo magical transpose that comes from character assign)
trials.Map = trials.Map';
trials.MoveDir = trials.MoveDir';

% Write the segment information to the MAT file
data.trialBlocks = trials.segInfo;
save([matpath filesep data.fname '.mat'],'data');

if (dbFlag > 1)
    H = MATOverview(data.fname);
    fprintf('Time window is from %2.2f to %2.2f, %2.2f min long\n',...
        tstart/60, tstop/60, (tstop-tstart)/60);
    fprintf('Trial type is %s\n',map);      % what the hell is this talking about? -CM 10/27/2013
    for i=2:length(H)
        plot(H(i),[tstart]/60, [0], 'g^',[tstop]/60, [0], 'gv')
        %xlim(H(i),[tset/60 - 1, tstop/60+1])
    end
    
    fprintf('Close figure window to continue\n');
    while (ishandle(H(1)))
        pause(0.1)  % twiddle thumbs until figure window closed
    end
end

catch me
    me.getReport()
    keyboard
end
    