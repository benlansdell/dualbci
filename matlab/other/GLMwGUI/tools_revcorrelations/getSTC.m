function [Mode1,Mode2,covdiff] = getSTC(Stim,spike_time,t,ppms,dt)
% This function generates the spike triggered covariance(STC) with given 
% input(stimulus) and output(spike train) and ppms(timestep).
%
% 1st August. 2013 
% Author : SangWook Lee
Stim = Stim(:);
sta = getSTA(Stim,spike_time,t,ppms,dt); % generate STA
Stim = Stim - mean(Stim);
Flength = t*dt*ppms;
dsamp = 1/dt; % sample rate
Stim = Stim(dsamp/2:dsamp:end);
spike_time = round(spike_time*dt);

% generate a stc
stc = zeros(Flength,Flength); % Flength := length of the filter.
pc = zeros(Flength,Flength); % Prior matrix

randInd = randsample(Flength+1:length(Stim),length(spike_time));
for i = 1:length(spike_time)
    chunk = Stim(randInd(i)-Flength+1:randInd(i));
    pc = pc + chunk*chunk';
end

for i = 1:length(spike_time) % spike_time vector contains the time of spike happens
    if spike_time(i) > Flength
        piece = Stim(spike_time(i)-Flength+1:spike_time(i)); % piece component for the each spike
        stc = stc + (piece-sta)*(piece-sta)'; % STC calculation
    end
end

stc = stc./length(spike_time); % spike triggered covariance matrix
pc = pc./length(spike_time); % covariance prior matrix
covdiff = stc-pc; % spike triggered difference matrix

if mean(mean(covdiff(round(length(covdiff)*(3/4)):length(covdiff),...
        round(length(covdiff)*(3/4)):length(covdiff)))) < 0
else
    covdiff = -covdiff;
end

[v,E] = eig(covdiff);

Mode1 = v(:,1);
Mode2 = v(:,2);

% choose the right sign of the feature
if mean(Mode1(round(length(Mode1)*0.2)+1:length(Mode1))) < 0
    Mode1 = -Mode1;
end

if mean(Mode2) > 0
    Mode2 = -Mode2;
end
covdiff = -covdiff;