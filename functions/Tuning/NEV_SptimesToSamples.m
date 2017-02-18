% NEV_SptimesToSamples.m
% by C. Matlack 2013-11-25
%
% USAGE:    nev = NEV_SptimesToSamples(nev, SR, catalog, toffset, dbFlag)
%
% INPUTS:   nev         NEV file data object
%           SR          desired sample rate
%           catalog     N x 2 matrix of [electrode, unit] listings 
%           toffset = 0 offset time from file start at which first sample
%                       period should start
% OUTPUTS:  nev         updated nev object
%           catIdx      mapping from input catalog to sampled catalog. If
%                       no existing sampled channels, this is just 1:N
%
% Converts spike times to sampled spike train. 
% Catalog = [chans sorts] specifies which channels to convert
% Toffset is an optional time offset to apply before sampling, e.g. in
% order to achieve perfect match-up of spikes per sample period with a
% Labview data set
%
% CREATES THE FOLLOWING DATA STRUCTURES WITHIN THE NEV STRUCTURE:
%
%   nev.Sampled.Spikes       binned spike counts
%              .ToffsetSec   time offset of first sample time
%              .Channels     [Electrode Unit] identities for spiken trains
%              .SampleRateHz sample rate for bins

function [nev,catIdx] = NEV_SptimesToSamples(nev, SR, catalog, toffset, dbFlag)
global NEVCACHE
tic
if (nargin < 4)
    toffset = 0;
else
    if (toffset > 1/SR)
        fprintf('Warning: time offset %f greater than one sample period %f.\n',...
            toffset, 1/SR);
    end
end

assert(sum(sum(catalog == round(catalog))) == numel(catalog),...
    'Catalog contains non-integer values, please fix.');

if (nargin < 5)
    dbFlag = false;
end

try
        % (Assume until exception thrown) 
        % at least some sampled channels available. Determine those that
        % must be generated new.
        assert(toffset == nev.Sampled.ToffsetSec && ...
            SR == nev.Sampled.SampleRateHz, 'Must regenerate all sampled spikes');
        newChans = setdiff(catalog, nev.Sampled.Channels, 'rows');
        if isempty(newChans)
            if (dbFlag)
                fprintf('Sampled spike trains already calculated, returning.\n')
            end
            [dummy,catIdx] = ismember(catalog, nev.Sampled.Channels,'rows');
            return
        end
catch me
    %me.getReport
    nev.Sampled.Spikes = [];
    nev.Sampled.ToffsetSec = toffset;
    nev.Sampled.TimeSec = (toffset:1/SR:ceil(nev.MetaTags.DataDurationSec*SR)/SR)';
    nev.Sampled.Channels = [];
    nev.Sampled.SampleRateHz = SR;
    newChans = catalog;
    % continue with function
    fprintf('(Re)generating sampled spike trains...\n')
    
end

%%% Preallocate new space for new catalog entries
try
nev.Sampled.Spikes = [nev.Sampled.Spikes, zeros(length(nev.Sampled.TimeSec),size(newChans,1))];
nev.Sampled.Channels = [nev.Sampled.Channels; newChans];





[dummy,catIdx] = ismember(newChans, nev.Sampled.Channels,'rows');
% ci indexes into newChans, and catIdx maps ci to existing catalog indices
for ci = 1:size(newChans,1)
    
    elec = newChans(ci,1);
    unit = newChans(ci,2);
    if (dbFlag)
        fprintf('%d-%d in slot %d, %d spikes...\n',...
            elec,unit,catIdx(ci), sum(nev.Data.Spikes.Electrode == elec & nev.Data.Spikes.Unit == unit));
    end
    nev.Sampled.Spikes(:,catIdx(ci)) = histc(...
        double(nev.Data.Spikes.TimeStamp((nev.Data.Spikes.Electrode == elec) & (nev.Data.Spikes.Unit == unit)))/...
        nev.MetaTags.TimeRes,...
        nev.Sampled.TimeSec)*SR;
end

% cut off last histogram value which is anything falling on last sample bin
% edge
%nev.Sampled.Spikes = nev.Sampled.Spikes(1:end-1,:)*SR;

% We've just counted spikes in bins, need to convert this to a sampled
% value of firing rate. To do this multiply bin count by number of bins/sec
% to get number of spikes/sec.

% mapping from requested channels to all sampled channels
[dummy,catIdx] = ismember(catalog, nev.Sampled.Channels,'rows');

if (dbFlag)
    fprintf('done in %.2fs\n',toc);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Save sampled data to NEV file cache to avoid duplication
NEV = nev;  % variable name used by openNEV()

[tf, idx] = ismember([NEV.MetaTags.FilePath filesep NEV.MetaTags.Filename], NEVCACHE.files);
if tf
    if (dbFlag)
        fprintf('Saving sampled spikes to NEV cache.\n')
    end
    NEVCACHE.nevs{idx} = NEV;
    NEVCACHE.lastRead{idx} = clock();
    openNEV();      % see if old cache items should be dumped to free mem
else
    fprintf('Cache miss after NEV loaded, this should not be possible...\n')

end

%save([nev.MetaTags.FilePath filesep nev.MetaTags.Filename '.mat'], 'NEV', '-v7.3');
catch me
    me.getReport
    keyboard
end

