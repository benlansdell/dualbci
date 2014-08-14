function [FiltteredResponse,PriorDist,FiltteredDist] = FeatSelect(stim,...
    spike_time,feature)
% This function generates the probability distribution of prior, and the
% filtered responses that are involved with the spike time.
%
% 31st July. 2013 
% Author : SangWook Lee
% reorganizing the data so that they are in column vector.
stim = stim(:);
feature = feature(:);

dt = 0.5;
len = 100;

% sampling
dsamp = round(1/dt); % sample rate
stim = stim(1:dsamp:end);
spike_time = round(spike_time*dt);

proj_stim = crosscorr_mex(stim',feature');
priorresp = [zeros(1,length(feature)),proj_stim(1:length(stim)-length(feature))];

spikeresp = priorresp(spike_time); % extract the part of response that are involved with the spike.

% edge data for generating probability distribution
edge = linspace(min(priorresp),max(priorresp),len);

% calculating the prior distribution
pprior = (histc(priorresp,edge)./norm(histc(priorresp,edge)))';

for i = 1:length(pprior)
    if pprior(i) == 0
        pprior(i) = 1e-5;
    end
end

% probability distribution of a response of each features
pspikes = (histc(spikeresp,edge)./norm(histc(spikeresp,edge)));
pspikes = pspikes./norm(pspikes);

FiltteredResponse = priorresp;
PriorDist = pprior';
PriorDist = PriorDist./norm(PriorDist);
FiltteredDist = pspikes;
FiltteredDist = FiltteredDist./norm(FiltteredDist);