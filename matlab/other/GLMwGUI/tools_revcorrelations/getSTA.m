function sta = getSTA(Stim,spike_time,t,ppms,dt)
% This function generates the spike triggered average(STA) with given 
% input(stimulus) and output(spike time) and ppms(timestep).

% 31st July. 2013 
% Author : SangWook Lee
Stim = Stim(:);
Stim = Stim - mean(Stim); % subtracting mean
Flength = t*dt*ppms;
dsamp = 1/dt; % sample rate
Stim = Stim(dsamp/2:dsamp:end);
spike_time = round(spike_time*dt);

% generate a sta
sta = 0;
for i = 1:length(spike_time)
    if spike_time(i) > Flength
        sta = sta + Stim(spike_time(i)-Flength+1:spike_time(i));
    end
end
sta = sta./length(spike_time);