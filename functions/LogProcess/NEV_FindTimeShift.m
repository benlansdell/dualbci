% NEV_FindTimeShift.m
% By C. Matlack 2013-12-02
%
% USAGE: ToffsetSec = NEV_FindTimeShift(samplesData, SRData, samplesNEV, SRnev, tRangeSec)
%
% This function calculates the squared error of samplesData and samplesNEV,
% with the latter histogrammed at SRData, and uses this to determine the
% exact time offset ToffsetSec when samplesNEV begins within samplesData.
% If the error is too high, the function returns an empty set.

% TODO: Parallelize this to take matrix for samplesData, column per
% separate channel
function ToffsetSec = NEV_FindTimeShift(samplesData, SRData, sptimesNEV, SRnev, tRangeSec, dbFlag)
try
if (nargin < 6)
    dbFlag = 0;
end

if isempty(sptimesNEV) || sum(isnan(tRangeSec))
    fprintf('No NEV spikes/time range provided, can''t help you!\n')
    ToffsetSec = -1;
    return
end

nCh = size(samplesData,2);
binWidthNEV = SRnev/SRData;
relTol = 0.5;              % relative tolerance for match

% will store errors. Initialize to big since we'll be finding min.
errs = zeros(length(tRangeSec),nCh);

ti = 1;
for tOffsetSec = reshape(tRangeSec,1,numel(tRangeSec))
    
    % time between Toffset and next sample bin @ SR, in units of nev
    % samples
    iNevOffset = SRnev*(1/SRData - mod(tOffsetSec,1/SRData));   

    % same time as above except with 1 sample period added, units of data
    iDataOffset = 1 + ceil(tOffsetSec*SRData); 
    
    % now we bin the nev data using this offset index
    nevBinned = histc(sptimesNEV,iNevOffset: SRnev/SRData : sptimesNEV(end));
    
    % discard last value which is any spike right at bin boundary; make col
    nevBinned = nevBinned(1:end-1)';
    
if (iDataOffset + length(nevBinned) > length(samplesData) || iDataOffset < 1)
    fprintf('Off the end of the samples, fail.\n')
    ToffsetSec = -1;
    return
end

    % now we can compare the NEV bin counts to the matched-up data samples
    errs(ti,:) = ...
        sum( (repmat(nevBinned,1,nCh) - ...
                samplesData(iDataOffset + (0:length(nevBinned)-1),:) ).^2 ...
           );
    
    % Plot progress
    if (dbFlag >= 3 )
        figure(1),clf
        subplot(2,1,1)
        h1 = plot(conv(nevBinned,gausswin(7),'same'),'k-.');
        hold on
        h2 = plot(conv2(samplesData(iDataOffset + (0:length(nevBinned)-1),:),gausswin(7),'same'));
        %legend(h2, 'NEV bin counts')
        xlim([0 250])
        ylim([0 4])
        
        subplot(2,1,2)
        h3 = plot(errs);
        drawnow
        pause(0.001)
    end
   ti = ti + 1; 
end
if (dbFlag >= 3)
    keyboard
end
% Return time offsets minimizing error only when error is driven exactly to
% 0.
[minErrs,iMinErrs] = min(errs);
tol = relTol*max(errs);
ToffsetSec = tRangeSec(iMinErrs).*(abs(minErrs) < tol) -1*(abs(minErrs) > tol);

if (dbFlag == 1 && ToffsetSec ~= -1) || dbFlag > 1
        figure(1),clf
        subplot(2,1,1)
        h1 = plot(conv(nevBinned,gausswin(7),'same'),'k-.');
        hold on
        h2 = plot(conv2(samplesData(iDataOffset + (0:length(nevBinned)-1),:),gausswin(7),'same'));
        %legend(h2, 'NEV bin counts')
        xlim([0 250])
        ylim([0 4])
        
        subplot(2,1,2)
        if (ToffsetSec == -1)
            plt = 'r-';
        else
            plt = 'g-';
        end
        plot([1 1]*iMinErrs, [0 1000],plt)

        hold on
        h3 = plot(errs);
        ylim([0 max(errs)])
        drawnow
        pause(0.001)
end
    
% if (dbFlag >= 2 )
%     max(errs)
%     min(errs)
% end

catch me
    me.getReport
    keyboard
end
