% NEV_ScreenForMinRate.m
% By C. Matlack 2013-11-25
%
% USAGE:    catalog = NEV_ScreenForMinRate(nev, minRateHz)
%
% Creates a data structure catalog = [ chans sorts ] 
% which contains the channel numbers and sorts codes (col vectors) of units
% with a firing rate of at least minRateHz, sorted descending by firing
% rate.
%
% This catalog is meant to be dropped into nev.UnitsByRate by
% NEV_SptimesToSamples, after verifying that sampled versions of each
% channel in the catalog exist.


function catalog = NEV_ScreenForMinRate(nev, minRateHz, dbFlag)

if (nargin < 3)
    dbFlag = false;
end

if (dbFlag)
fprintf('Finding good channels + sorts with at least %.1f Hz spike rate...\n', minRateHz)
tic
end

% try
%     if (~isempty(nev.UnitsByRate))
%         if (dbFlag)
%             fprintf('Unit catalog already exists, returning...\n')
%         end
%         return
%     end
% catch me
%     % isempty will throw exception if nev.UnitsByRate doesn't exist
% end

% Combine sort codes and electrode numbers into a single unique channel
% number. Note that an entry will exist for every spike
channelData = [10*nev.Data.Spikes.Electrode + uint16(nev.Data.Spikes.Unit)];

% Generate list of channels with spikes
channels = unique(channelData)';

% Count the number of spikes on each channel to calculate average rates
spikeAvg = zeros(size(channels));
for i=1:length(spikeAvg)
    spikeAvg(i) = length(find(channelData == channels(i)))/nev.MetaTags.DataDurationSec;
end

% Pull out index set of cells firing above minimum rate
fastCellsIdx = find(spikeAvg >= minRateHz);
Nunits = length(fastCellsIdx);

% Sort descending by firing rate
[dummy,sortIdx] = sort(spikeAvg(fastCellsIdx),'descend');
fastCellsIdx = fastCellsIdx(sortIdx);

% write channel & sort code list to data structure
% typecast is because native type is uint16, can't do multiplication
catalog = single([channels(fastCellsIdx)/10, mod(channels(fastCellsIdx),10)]);

if (dbFlag)
fprintf('Elpased time: %3.fms\n',toc*1000);
fprintf('Elec\tSort\tAvg rate\n')
disp([channels(fastCellsIdx)/10, mod(channels(fastCellsIdx),10), spikeAvg(fastCellsIdx)])
end

% fprintf('Saving back to NEV MAT\n')
% NEV = nev;
% save([nev.MetaTags.FilePath filesep nev.MetaTags.Filename '.mat'],'NEV');
end