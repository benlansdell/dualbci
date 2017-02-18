function spikes = SptimesToSpikes(sptimes, tt)
% by C. Matlack 2012-07-26
%
% USAGE:    spikes = SptimesToSpikes(sptimes, tt)
% Converts array of spike times to binned samples in bins ending at tt.
%
% All spikes prior to tt(1) assumed to be in first bin.
%
% Output is row vector.

%tt = (floor(min(sptimes)*SR):1:ceil(max(sptimes)*SR))/SR;
spikes = 0*tt;

for ti=1:length(tt)
    if (ti == 1)
        spikes(ti) = sum(sptimes <= tt(ti));
    else
        spikes(ti) = sum((sptimes > tt(ti-1)) & (sptimes <= tt(ti)));
    end
end

if (sum(spikes) < length(sptimes))
    fprintf('WARNING (SptimesToSpikes): %d spikes not captured by sample time vector\n',length(sptimes)-sum(spikes))
    keyboard
end