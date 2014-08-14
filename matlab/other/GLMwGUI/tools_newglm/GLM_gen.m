function [Stimfilt,Histfilt,negloglival] = GLM_gen(DataStruct)
% This function perform the glm fitting using maximum likelihood method.
% The input 'DataStruct' should contain as it is specified in the
% instructor document.
%
% 2nd August 2013
% Author : SangWook Lee.
timeb4spk = DataStruct.timeb4spk;
timeafterspk = DataStruct.timeafterspk;
dt = DataStruct.dt;
spike_time = DataStruct.spike_time;
Stim = DataStruct.stim;
ppms = DataStruct.ppms;
Stim = Stim(:);

Stim = Stim - mean(Stim);
Flength = timeb4spk*ppms*dt;
dsamp = 1/dt; % sample rate
Stim = Stim(dsamp/2:dsamp:end);
spike_time = round(spike_time*dt);

sta = 0;
for i = 1:length(spike_time)
    if spike_time(i) > Flength
        sta = sta + Stim(spike_time(i)-Flength+1:spike_time(i));
    end
end
sta = sta./length(spike_time);

GlmStatBar = waitbar(0,'Calculating GLM filters');
gg0 = makeFittingStruct_GLM(sta,dt,timeafterspk,ppms);  % projects sta into basis for fitting k
waitbar(0.25,GlmStatBar,'Calculating GLM filters');

gg0.tsp = spike_time;  % Insert spikes into fitting struct
gg0.tspi = 1;   % First spike to use (you can ask it to ignore the first "n" spikes)
                
[logli0,rr0,tt] = neglogli_GLM(gg0,Stim); % Compute logli of initial params (if desired)
waitbar(0.35,GlmStatBar,'Calculating GLM filters');

% Do ML estimation of model params
opts = {'display', 'iter', 'maxiter', 100};
[gg,negloglival] = MLfit_GLM(gg0,Stim,opts,GlmStatBar); % do ML (requires optimization toolbox)
waitbar(1,GlmStatBar,'Calculating GLM filters');

Histfilt = gg.ihbas*gg.ih;
Stimfilt = gg.k;

Histfilt = Histfilt(1:end-1);
close(GlmStatBar)