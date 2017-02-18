% FixCursorHistory.m
%
% Due to bug in LabVIEW code, recorded cursor trajectory is meaningless.
% This code reconstructs it by integrating from the recorded starting
% position of each trial. Cursor position before the first trial is
% unknown. We can dump out integration error at the beginning of each trial
% to see if there is an inconsistency with the reconstruction

function dataFixed = FixCursorHistory(data, dbFlag)

if (nargin < 2)
    dbFlag = 0;
end

% set cursor trajectory back to 0
data.stateHist.cursor = zeros(length(data.stateHist.time),2);
% post-processed flag
data.stateHist.valid = logical(data.stateHist.cursor(:,1));  
data.trials.valid = false(size(data.trials.time));


if (dbFlag > 1)
    h1 = figure;
    clf
    subplot(2,1,1), title('Cursor trajectory in screen coordinates');
    %subplot(3,1,2), title('Cursor trajectory in error coordinates');
    subplot(2,1,2), title('Norm of cursor error');
    h2 = figure;
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
    
    % Time interval is from beginning of current trial to beginning of
    % next.
    idxRange = find(data.stateHist.time > data.trials.time(i) - data.trials.duration(i) & ...
                data.stateHist.time < data.trials.time(i+1) - data.trials.duration(i+1));
    %idx = find(data.stateHist.time <= data.trials.time(i),1,'last');
    idx = idxRange(1);
    data.stateHist.cursor(max(1,idx-1),:) = data.trials.startPos(i,:);     % note that the -1 causes integrator error to be reset
    
    % find out whether integrator was enabled on each output channel
    pvmapIdx = find(data.pvmap.time <= data.trials.time(i),1,'last');
    paramIdx = find(data.params.time <= data.trials.time(i),1,'last');
    intDisable = data.pvmap.intDisable(pvmapIdx,:);
    
    hitWall = false;
    
    % integrate velocity until start of next trial
    %while (data.stateHist.time(idx) < data.trials.time(i+1))
    for idx = idxRange(1):idxRange(end)     % fucking stupid syntax
        
        % iterate over dimensions
        for dim = 1:2
            % integrate one step, conditional on integrator enable
            if (intDisable(dim))
                val = data.stateHist.velocity(idx,dim);
            else
                val = data.stateHist.velocity(idx,dim)/data.sampleRate + ...
                    data.stateHist.cursor(max(1,idx-1),dim);
            end
               
            % enforce screen boundaries
            val = max(-0.5,val);
            val = min(0.5,val);
            if (val == 0.5 || val == -0.5)
                hitWall = true;
            end
            data.stateHist.cursor(idx,dim) = val;
        end
        
        %idx = idx+1;
    end % of samples in trial
    
        data.stateHist.error(idxRange,:) = data.stateHist.cursor(idxRange,:)-...
        repmat(data.trials.targetPos(i,:),length(idxRange),1);
    
        warn = false;
        if (norm(data.stateHist.cursor(idx,:)-data.trials.startPos(i+1,:)) > 0.1)
            warn = true;
        else
            data.stateHist.valid(idxRange) = true;
            data.trials.valid(i) = true;
        end
        
        if (dbFlag > 0 && warn)
            % calculate accumulated integrator error
            fprintf('Trial #%d, %2.1fs, int err (%1.2f, %1.2f), wall %d, intDis = (%d,%d), win = %d\n',...
            i, data.trials.duration(i), data.stateHist.cursor(idx,1)-data.trials.startPos(i+1,1),...
            data.stateHist.cursor(idx,2)-data.trials.startPos(i+1,2),hitWall,...
            intDisable(1),intDisable(2),data.trials.success(i));
        end
        

    
    if (dbFlag > 1)
        if (~ishandle(h1))
            h1 = figure;
        end
        figure(h1)
        subplot(3,1,1), hold on
        dim = 1;
        plot(data.stateHist.time(idxRange),data.stateHist.cursor(idxRange,dim))
        plot(data.stateHist.time(idxRange(end)),data.stateHist.cursor(idxRange(end),dim),'og')
        plot(data.trials.time(i+1)-data.trials.duration(i+1),data.trials.startPos(i+1,dim),'+r')
        title('X trajectory')
        
        subplot(3,1,2), hold on
        dim = 2;
        plot(data.stateHist.time(idxRange),data.stateHist.cursor(idxRange,dim))
        plot(data.stateHist.time(idxRange(end)),data.stateHist.cursor(idxRange(end),dim),'og')
        plot(data.trials.time(i+1)-data.trials.duration(i+1),data.trials.startPos(i+1,dim),'+r')
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

fprintf('%2.1f%% of all trajectory data successfully reconstructed.\n',...
    100*length(find(data.stateHist.valid))/length(data.stateHist.valid));
fprintf('~%2.1f%% of in-trial trajectories reconstructed.\n',...
    100*length(find(data.stateHist.valid))/(sum(data.trials.duration)*data.sampleRate));

dataFixed = data;

end