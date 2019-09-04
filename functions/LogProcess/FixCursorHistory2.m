% FixCursorHistory.m
%
% Due to bug in LabVIEW code, recorded cursor trajectory is meaningless.
% This code reconstructs it by integrating from the recorded starting
% position of each trial. Cursor position before the first trial is
% unknown. We can dump out integration error at the beginning of each trial
% to see if there is an inconsistency with the reconstruction

% 2013-02-03 This version is for the post- 2012-12-23 code fix that records
% actual cursor position correctly.

function dataFixed = FixCursorHistory2(data, dbFlag)

if (nargin < 2)
    dbFlag = 0;
end


% post-processed flag
data.stateHist.valid = logical(data.stateHist.cursor(:,1));  
data.trials.valid = false(size(data.trials.time));


if (dbFlag > 1)
    h1 = figure(1);
    clf
    subplot(2,1,1), title('Cursor trajectory in screen coordinates');
    %subplot(3,1,2), title('Cursor trajectory in error coordinates');
    subplot(2,1,2), title('Norm of cursor error');
    h2 = figure(2);
    clf
else
    h1 = NaN;
    h2 = NaN;
end

for i = 1:(length(data.trials.time)-1)
    % (re)-initialize cursor position with trial start position
    % (LabVIEW has clocks everywhere and is not hard-realtime, so
    % recorded event times may not be identical to sample times)
    % NOTE: trial.time is *END* of trial, when trial data recorded
    
    % Time interval is from beginning of current trial to 
    % time of trial end
    
    idxRange = find(data.stateHist.time < data.trials.time(i) - data.trials.duration(i),1,'last'):...
                find(data.stateHist.time < data.trials.time(i),1,'last'); % - data.trials.duration(i+1),1,'last');
    %idx = find(data.stateHist.time <= data.trials.time(i),1,'last');
    idx = idxRange(1);
    
    %data.stateHist.cursor(max(1,idx-1),:) = data.trials.startPos(i,:);     % note that the -1 causes integrator error to be reset
    
    % find out whether integrator was enabled on each output channel
    pvmapIdx = find(data.pvmap.time <= data.trials.time(i),1,'last');
    paramIdx = find(data.params.time <= data.trials.time(i),1,'last');
    intDisable = data.pvmap.intDisable(pvmapIdx,:);
    
    hitWall = false;
    
        % LabVIEW bugfix 2012-12-23 writes cursor position to log where
        % error used to go (of course error was wrong)
%        data.stateHist.cursor(idxRange,:) =
%        data.stateHist.error(idxRange,:); (done in LogRead.m already)

        data.stateHist.error(idxRange,:) = -data.stateHist.cursor(idxRange,:)...
        +repmat(data.trials.targetPos(i,:),length(idxRange),1);
    
    % tweak 2013-11-20 flip sign of error if it starts off negative
%     if (data.stateHist.error(idxRange(2),1) < 0)
%         data.stateHist.error(idxRange,1) = -data.stateHist.error(idxRange,1);
%     end
    
    % the condition for a valid trial is that one of the first three
    % points in the current trajectory agrees with startPos.
    % Typically there can be a one time sample error in when each was
    % recorded, which is not a problem.
    
    warn = false;
    if (min(sqrt(sum((data.stateHist.cursor(idxRange(1:min(end,3)),:)-repmat(data.trials.startPos(i,:),min(3,length(idxRange)),1)).^2,2))) > 0.02)
        fprintf('Start position disagrees by more than 0.02\n');
        warn = true;
    else
        data.stateHist.valid(idxRange) = true;
        data.trials.valid(i) = true;
    end
    
    if (dbFlag > 1 && warn)
        if (~ishandle(h1))
            h1 = figure;
        end
        figure(h1)
        subplot(3,1,1), hold on
        dim = 1;
        l0 = plot(data.stateHist.time(idxRange),data.stateHist.cursor(idxRange,dim));
        plot(data.stateHist.time(idxRange(1)),data.stateHist.cursor(idxRange(1),dim),'bo')
        l1 = plot(data.stateHist.time(idxRange(end)),data.stateHist.cursor(idxRange(end),dim),'og');
        l2 = plot(data.trials.time(i)-data.trials.duration(i),data.trials.startPos(i,dim),'+r');
        l3 = plot(data.trials.time(i),data.trials.targetPos(i,dim),'xr');
        legend([l0 l1 l2 l3], 'Trajectory','Traj End','startPos','Target Loc')
        title('X trajectory')
        
        subplot(3,1,2), hold on
        dim = 2;
        plot(data.stateHist.time(idxRange),data.stateHist.cursor(idxRange,dim))
        plot(data.stateHist.time(idxRange(1)),data.stateHist.cursor(idxRange(1),dim),'bo')
        plot(data.stateHist.time(idxRange(end)),data.stateHist.cursor(idxRange(end),dim),'og')
        plot(data.trials.time(i)-data.trials.duration(i),data.trials.startPos(i,dim),'+r')
        plot(data.trials.time(i),data.trials.targetPos(i,dim),'xr')
        title('Y trajectory')
        
        subplot(3,1,3), hold on
        plot(data.stateHist.time(idxRange),...
            sqrt(data.stateHist.error(idxRange,1).^2 + data.stateHist.error(idxRange,2).^2))
        plot(data.stateHist.time(idxRange(1))+[0 data.trials.duration(i)], ones(1,2)*data.params.successRadius(paramIdx),':')
        if (data.trials.success(i))
            plt = '^-';
        else
            plt = 'v-';
        end
        plot(data.trials.time(i)-[data.params.dwellTime(paramIdx) 0], ones(1,2)*data.params.successRadius(paramIdx),plt)
        title('Norm of position error')
 
        
        % display other exciting signals in another figure window
        if (~ishandle(h2))
            h2 = figure;
        end
        figure(h2)
        subplot(4,1,1)
        hold on
        plot(data.stateHist.time(idxRange),data.stateHist.spikes(idxRange,:))
        ylabel('Spikes/bin')
        xlabel('Time (s)')
        title(sprintf('Neuron count = %d\n',data.neurons));
        
        subplot(4,1,2)
        hold on
        plot(data.stateHist.time(idxRange),data.stateHist.rates(idxRange,:))
        ylabel('Spike rates (pps)')
        
        subplot(4,1,3)
        hold on
        plot(data.stateHist.time(idxRange),data.stateHist.velocity(idxRange,:))
        xlabel('Velocity Command')
        
        subplot(4,1,4)
        hold on
        plot(data.stateHist.time(idxRange),data.stateHist.error(idxRange,:))
        ylabel('Output Error')
        drawnow
        keyboard
        
    end

 
    
end % trial loop

dataFixed = data;

end
