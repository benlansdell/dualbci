% Overview.m
%
% General analyses of an experiment log file 
%
% INPUT: fname  filename of .mat file, no extension needed
%        path   optional path to .mat file, default is current dir
%
% OUTPUT: h     handle to figure window
%         h1-h4 handles to subplots

function H = MATOverview(fname, path)
if (nargin == 1)
    path = pwd;
end

oldpath = pwd;

if (nargin == 0)
% Select file here...
[fname,path] = uigetfile('*.mat');    
end

cd(path) 
% Open file, read only
load(fname)
cd(oldpath)

% assumes the existence of data structure called data with member struct
% trials.
%   time        Session time in seconds
%   index       global (for that day) trial index
%   startPos    starting cursor position (2D)
%   duration    trial time is sec = time limit if failed
%   targetPos   target position (2D)
%   threshold   effective target radius
%   success     1 if trial successful

% convert time to minutes since experiment start
%times = (data.trials.time - data.trials.time(1))/60;    

wins = find(data.trials.success == 1);
losses = find(data.trials.success == 0);

dwell = data.params(end,end);
if (data.params.time(end) > data.trials.time(1))
    disp('Warning: task parameters changed during session');
end

%% Iterate over continuous sessions (delineated by run/pause)
%for tStartIdx = data.

%% Calculate moving average performance in trials/min
tavg = 120;  % seconds over which to calculate performance average

[tavg, ravg, idx] = TrialsPerMin(data, data.messages.time(1), data.messages.time(end), tavg);

hfig = figure(10);
clf
set(gcf,'Name',fname)

h1 = subplot(4,1,1);

runTimes = data.messages.time([data.messages.type] == 1 ...
                                & [data.messages.data(:,1)] == 1);
%runTimes = [runTimes; data.messages.time(end)]; % fake final 'run event'

pauseTimes = data.messages.time([data.messages.type] == 0 ...
                                & [data.messages.data(:,1)] == 1);
%pauseTimes = pauseTimes(2:end); % system is initialized into pause state


% TODO create & plot shaded polygons when system is paused
% pairs = [pauseTimes, runTimes]/60;
% for i = 1:size(pairs, 1)
%     fill([pairs(i,:) fliplr(pairs(i,:))], [0 0 10 10], [.75 .75 .75])
%     hold on
% end

plot(data.trials.time(idx)/60, tavg,'b',data.trials.time(idx)/60, ravg,'g')
hold on
stem(data.trials.time/60, data.trials.success,'k');
xlabel('Time (min)');
title('Trials summary');
legend('Avg trials/min','Avg rewards/min', 'Trial events');
limits = axis();
%limits = [min(data.pvmap.time(end)/60,limits(1)) limits(2)];
xlim(limits(1:2));

figure(11)
clf
plot(tavg,ravg,'o')
hold on
h = plot([0 max(ravg)],[0 max(ravg)],'k-');
legend(h,'Rewards = Trials (prob = 1)')
ylabel('Rewards')
xlabel('Correct Trials')
title('Cheater plot: slope should equal reward probability')
figure(10)

%% Retrieve samples of baseline rates recorded when new parameters applied
h2 = subplot(4,1,2);
hn = plot(repmat(data.pvmap.time/60,1,data.neurons), -data.pvmap.offsets,'.-');
hold on

% calculate moving average baseline rates (note this is *not* the way they
% are calculated online)
tavg = 120; % 120s moving window average
window = ones(tavg*data.sampleRate,1);
window = window/sum(window);
avgRates = conv2(data.stateHist.rates,window,'same');
ha = plot(repmat(data.stateHist.time/60,1,size(data.stateHist.rates,2)),avgRates,':');
time = data.stateHist.time; % for return argument

xlim(limits(1:2))
title('Baseline firing rates & torques')
xlabel('Time (min)')
ylabel('Spikes/sec')
legend(hn,'1 avg','2 avg','3 avg','4 avg')
    
% TODO: gray shading on plot when trials are paused

%% DO it better, by recalculating actual spike averages


%% Plot graphical representation of pop vector mapping over time
h3 = subplot(4,1,3);
hold on
for i=1:size(data.pvmap.angles,2)
    enabled = find(data.pvmap.enabled(:,i));
    quiver(data.pvmap.time(enabled)/60,ones(length(enabled),1)+i-1,...
        cos(data.pvmap.angles(enabled,i)*pi/180),...
        sin(data.pvmap.angles(enabled,i)*pi/180),...
        .1)
end
grid on
% loop is repeated so that color cycling works nicely for quiver plot
for i=1:size(data.pvmap.angles,2)
h=plot(data.pvmap.time(data.pvmap.enabled(:,i)==1)/60, ones(length(find(data.pvmap.enabled(:,i))),1)*i,'ok',...
    data.messages.time/60,zeros(size([data.messages.time])),'^k');
end
xlim(limits(1:2))
ylim([0 7])
xlabel('Time (min)')
ylabel('Neuron index')
title('Population vector mapping over time')
legend(h,'Enabled neurons','Log Events')

%%
h4 = subplot(4,1,4);
stem(data.params.time/60,data.params.yrange(:,2))
hold on
stem(data.params.time/60,data.params.yrange(:,1))
xlim(limits(1:2))
xlabel('Time (min)')
ylabel('Y target range')
title('Game parameters over time')

H = [hfig h1 h2 h3 h4];    % returned array of figure handles

%figure
%FiringRateStats(data,3);