% NEVLoadForMAT.m
% by C. Matlack 2013-11-03
%
% USAGE: [chans, dT, fnames] = NEVFindTimesForMAT1b(data, matpath, nevpath, metaData, <usemat>, <dbFlag>)
%
% DEPENDS:  NPMK 2.3.0  (package of NEV/NS3 utility functions from
%                        Blackrock)
%
% <usemat> = false      tells openNEV() to save a MAT file for each NEV
%                       loaded if none already exist, speeding the loading 
%                       process the next time.
%
% <dbFlag> = false      prints extra debugging information
%
% This function identifies the NEV files recorded simultaneously with the
% provided LabVIEW data structure (contents of a .MAT file).
%
% For each NEV file recorded on the same day with the same animal name, it
% finds the fastest neurons and runs a cross-covariance between them and the
% spike channels in the LabVIEW data. Using the results it determines
% channel assignments and timestamp mappings between the two data
% structures. 
%
% Algorithm idiosyncracy: covariance is fast but silly, used to identify a
% candidate time offset. Squared error between smoothed spikes trains is
% used to confirm it's real (i.e. we're even on the right channel).
%
% Note that the Blackrock system clock is grossly incorrect,
% and timestamps within LabVIEW files have been known to disagree by an
% hour with the filename timestamp, the latter possibly due to Daylight
% Savings. 
%
% What is returned is an array in which each row is the NEV channels on
% that LabVIEW channel, the time offset, and the NEV filename. The array is
% split into an N x 4 array, a column vector, and a cell array of
% filenames.


% TODO: preallocate output data structure, start from last NEV file and
% work back since less channel shuffling at end of day. Set a time marker
% anytime all 4 chans agree within 1 sec, and use to trim size of stateHist
% that is checked.

% VERSION 1b: Based on 1a, which uses FindTimeShift to finely adjust and
% discover time offset using error cost function. New in this version, as
% soon as a channel is matched for a particular NEV file, the channel loop
% is re-run *only* doing the fine adjustment to detect remaining channels.


function data = NEV_FindTimesForMAT1b(data, matpath, nevpath, metaData, Tguess, dbFlag)
try
Nspikes = 500;     % How many spikes to use when calculating correlations
                   % (Just has to be good enough for unique
                   % identification of a channel matchup)
minRateHz = .5;     % minimum overall spike rate to consider a channel

Tgiveup = 15;      % sec w/o progress during channel search to give up, move on
Ngiveup = 50;       % number of elec/unit combos to check before moving on
NstartGiveup = 15;  % number to try if nothing yet found
nMaxBins = 100;      % maximum number of bins to check
nMaxSpread = 30;    % maximum range of sample point candidates allowed into brute force search
nConverge = 3;      % Number of consective bin counts with only one match before calling it good.
convergeWidth = 10;
tTrimSec = 5;       % Number of seconds @ start & end of NEV to trim off in
                    % order to ignore poorly-timed channel switches
ToffsetSpread = 2;  % maximum distribution of time offsets for 'ok' match (IN SAMPLES)

if (nargin < 5)
    Tguess = 15.5*60;    % guess for clock offset between data & nev
end

if (nargin < 6)
    dbFlag = false;
end

correctIDcount = 0;     % number of cells ID'd previously
skew = Inf;             % integrated time uncertainty


%%% Find out if this data already exits in the MAT file %%%%%%%%%%%%%%%%%%%
%if ( exist('data.nev','var') )  %% <--- this always returns false
nChanExisting = 0;
try
    %%% Check whether the MAT file was updated more recently than this file
    codeInfo = dir([mfilename('fullpath') '.m']);
    if (data.nevUpdateNum > codeInfo.datenum && ~isempty(data.nev))
        fprintf('NEV mapping already calculated, returning...\n');
        return
    end
    fprintf('MAT %s out of date with script %s,\nneed to rebuild NEV mapping...\n',...
        datestr(data.nevUpdateNum), datestr(codeInfo.datenum));
    correctIDcount = sum(sum(cell2mat({data.nev.chans}') > 0));     % number of cells successfully identified previously
    Toff = cell2mat({data.nev.Toffset}');
    skew = sum(diff(Toff(:,2:3)'));         % total uncertainty in offset times
    nChanExisting = sum(sum(cell2mat({data.nev.chans}) ~= 0));
catch me
    if (nargin < 5)
        dbFlag = 1;
    end
end

%%% Find metaData entry matching this data object %%%%%%%%%%%%%%%%%%%%%%%%%
% This will be used to figure out the predominant mapping during an nev
% file
midx = find(ismember(metaData.fileList, data.fname), 1);

if (isempty(midx))
    fprintf('No metaData entry for %s\n',data.fname)
    return
end


%%% Extract animal name and date from filename %%%%%%%%%%%%%%%%%%%%%%%%%%%%

name = regexp(data.fname,'([a-zA-Z]+)','match');
name = name{1};
date = data.fname(strfind(data.fname,'_')+1:end);
year = date(1:4);
month = date(6:7);
day = date(9:10);
pattern = [year month day name];

SR = data.sampleRate;
Ndata = length(data.stateHist.time);


%%% Find all NEV files matching animal name and date %%%%%%%%%%%%%%%%%%%%%%

fprintf('\nFor MAT file %s,\n looking for NEV matching %s\n in %s\n',...
    data.fname, pattern, nevpath);

dirData = dir([nevpath filesep pattern '*.nev']);
%assert(~isempty(dirData),'No NEVs match %s\n\n',[nevpath filesep pattern '*.nev'] )


matchIdx = false(length(dirData),1);
for fi = 1:length(matchIdx)
    matchIdx(fi) = strncmp(pattern,dirData(fi).name,length(pattern));
end

fileList = {dirData(matchIdx).name};

if (isempty(fileList))
    fprintf('I found no matching NEV files; check nevpath\n\n')
    return
end

%%% Preallocate entire data structure


%% Now iterate over the matching files, looking for matches with data %%%%%
fnames = {};
chans = [];
Toffsets = [];
tOffsetGuess = NaN;     % time offset (known for sure) for nev relative to data

for file = fliplr(fileList)                     % work backwards
    file = file{1};     % make not cell array
        
    nev = openNEV(strcat(nevpath,filesep,file));

    fi = find(ismember(fileList,{file}));        % index of file in array
    fprintf('\n===%s\t%.2fmin\n',file,nev.MetaTags.DataDurationSec/60);
  
    SRnev = single(nev.MetaTags.SampleRes);

    data.nev(fi).Toffset = [0 0 0]; % pre-set value in case total match fail
    data.nev(fi).chans = [0 0 0 0];
    data.nev(fi).nevfile = [];
    data.nev(fi).DurationSec = 0;
    
    
    %%% Sort all neural channels by firing rate %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is specifically written to minimize memory usage since there are
    % a LOT of spikes.
    try
        spikeData = single(10*nev.Data.Spikes.Electrode + uint16(nev.Data.Spikes.Unit));
    catch me
        me.getReport
        keyboard
    end
    activeChans = unique(spikeData);
    spikeCounts = hist(spikeData,activeChans);
    [spikeCounts,sortIdx] = sort(spikeCounts,'descend');
    catalog = sortIdx(spikeCounts/nev.MetaTags.DataDurationSec > minRateHz);
    if ~isnan(tOffsetGuess)
        activeChans = activeChans(catalog);
    else
        activeChans = activeChans(catalog(1:min(end,NstartGiveup)));
    end
        % activeChans now contains each channel+sort code, sorted descending by
    % total number of spikes, only channels above minRateHz
        
    % if we got any matches the last go-around
    if (~isempty(chans) && sum(chans(end,:) ~= 0))
        % search the channels we matched in the last file first
        activeChans = [unique(chans(end,chans(end,:) ~= 0)), activeChans(~ismember(activeChans,chans(end,:)))];
        
        if (dbFlag >= 1)
            fprintf('\n***Expecting: %s\n',num2str(chans(end,:)))
        end
        % Use the time offset previously detected, by correcting for offset
        % between NEV file start times
        dDateTimeRaw = nev.MetaTags.DateTimeRaw - lastDateTimeRaw;
        %tOffsetLast = data.nev(fi-1).Toffset(1)/SR;  (this assumes order
        % of visiting files and so is fragile)
        
        dTimeSec = dDateTimeRaw*[0 0 0 0 60*60 60 1 .001]';
        tOffsetGuess = dTimeSec + tOffsetLast;
        
        fprintf('Last NEV offset to data is %dm %2.4fs; NEV diff is %.0fs;\n',...
            floor(tOffsetLast/60), mod(tOffsetLast,60), dTimeSec);
        fprintf(' so current offset to data should be %dm %2.4fs...\n',...
            floor(tOffsetGuess/60), mod(tOffsetGuess,60));
        if tOffsetGuess < 0
            fprintf('Next file appears to be before data record, aborting.\n')
            continue
        end
    else
        tOffsetGuess = NaN;
    end
    lastDateTimeRaw = nev.MetaTags.DateTimeRaw;

    %     if (dbFlag > 1)
    %         activeChans/10'
    %     end
    
    
    %%% For each channel, see if it correlates with %%%%%%%%%%%%%%%%%%%%%%%
    %   one of the data channels. We can quit as soon as we've figured
    %   out each labview channel.
    chanSet = zeros(1,4);
    Tset = zeros(1,4);
    tiStart = 1;                % data sample index to start search
    tic
    
%%%% OPTIMISTIC FAST PASS - Assume a guessed offset time is correct
secondPass = true;
if ~isnan(tOffsetGuess) && tOffsetGuess > 0     % catches condition of multiple data records per day
    fprintf('First fast pass\n')
    for ach = activeChans
        elec = round(ach/10);
        unit = ach - 10*elec;
        
        
        
        % get indices of first Nspikes matching this channel
        spikeIdx = find(spikeData == ach, Nspikes);
        if (length(spikeIdx) < Nspikes && dbFlag)
            %fprintf('Warning: only found %d spikes on channel %2.1f\n',...
            %    length(spikeIdx), ach/10);
        end
        
        % retrieve the spike times from the nev data, uint32 -> float32
        sptimesNEV = single(nev.Data.Spikes.TimeStamp(spikeIdx));
        
        
        for ch = find(chanSet == 0)
            if isnan(tOffsetGuess)
                tGuessSec = mean(Tset(Tset > 0))/SR;
            else
                tGuessSec = tOffsetGuess;
            end
if (dbFlag >= 1)
            fprintf('Checking %d-%d -> %d\n',elec,unit,ch);
end
            tOffsetRangeSec = (tGuessSec-0.1):0.0005:(tGuessSec+0.1);
            if ach == chans(end,ch) && dbFlag
                Toffsets = NEV_FindTimeShift(data.stateHist.spikes(:,ch),SR,sptimesNEV,SRnev,tOffsetRangeSec,1);
                if Toffsets == -1 && (dbFlag >= 2)
                    NEV_FindTimeShift(data.stateHist.spikes(:,ch),SR,sptimesNEV,SRnev,tOffsetRangeSec,2);
                    
                end
            else
                Toffsets = NEV_FindTimeShift(data.stateHist.spikes(:,ch),SR,sptimesNEV,SRnev,tOffsetRangeSec,dbFlag);
            end
            
            if (Toffsets ~= -1)     % winner
                secondPass = false;
                iOffsetData = Toffsets*SR;
                fprintf('Elec %d Unit %d <==> LabVIEW ch %d ',...
                            elec, unit, ch);
                        fprintf('Toffset %dm %.4fs\n',...
                            floor(iOffsetData/SR/60), mod(iOffsetData/SR,60));
                        chanSet(ch) = ach;
                        Tset(ch) = iOffsetData;
                        if sum(isnan(Tset))
                            fprintf('Getting nonsense times; bailing on whole day\n\n')

                            return
                        end
            end
        end

    end
    else
        fprintf('Skipping fast pass, no offset guess!\n')
end

    
%%%%%% MAIN PASS USING BIN COUNT MATCHING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    firstPass = true;
    fprintf('Main slow pass using bin count matching...\n')

    for ach = activeChans
        % bail out of this as soon as we've ID'd one channel
    if ~sum(chanSet) && ~(tOffsetGuess > 0)     % covers NaN case for tOffsetGuess     
        elec = round(ach/10);
        unit = ach - 10*elec;
 if (dbFlag >= 1)      
        fprintf('Checking %d-%d\n',elec,unit);
 end       
        % get indices of first Nspikes matching this channel
        spikeIdx = find(spikeData == ach, Nspikes);
        if (length(spikeIdx) < Nspikes && dbFlag)
            %fprintf('Warning: only found %d spikes on channel %2.1f\n',...
            %    length(spikeIdx), ach/10);
        end
        
        % retrieve the spike times from the nev data, uint32 -> float32
        sptimesNEV = single(nev.Data.Spikes.TimeStamp(spikeIdx));
        
        % Trim off spikes at beginning & end in case channels were changed
        sptimesNEV = sptimesNEV(sptimesNEV > tTrimSec*SRnev & ...
                   sptimesNEV < (sptimesNEV(end) - tTrimSec*SRnev));
        
        % find all gaps longer than 2 labview sample periods
        iGapSetNev = find(diff(sptimesNEV) > 2*SRnev/SR);    
        
        if (length(iGapSetNev) < 2)
            fprintf('Not enough spike gaps for me to work with')
            continue
        end
        
        % time, in nev samples, at end of each gap
        iGapEndSetNev = sptimesNEV(iGapSetNev + 1);  
        
        % time, relative to NEV file, of the *end* of the *first* gap
        tFirstGapEndSec = iGapEndSetNev*SRnev;
        
        % time offset, within NEV, of splits, in units of data samples
        iBinOffsetData = iGapEndSetNev(1)*SR/SRnev;

        % number of data sample periods in the biggest gap
        gapIdx = length(iGapSetNev);
        
        % 'splits' are continuous stretches of spikes we'll look at
        bins = [];
        binIdx = 1;

        while (~isempty(gapIdx) && gapIdx > 1 )
            
            % TODO: This could be tweaked to get one number that always
            % works
            % number of data sample periods in this bin
            
            %number of NEV sample periods in this bin
            binTime = (iGapEndSetNev(gapIdx) - 2*SRnev/SR - iGapEndSetNev(1));    
            bins(binIdx).periodCount = ceil( binTime*SR/SRnev);
            %bins(binIdx).periodCount = ceil((gapTimes(gapIdx) - 2*SRnev/SR - gapTimes(1))...
            %                    / (SRnev/SR));
            % number of spikes in this split (bounded by gaps)
            
            
            % the -1 time shift is to capture exactly the first spike on
            % the boundary, but not the spike that terminates the second
            % gap
            bins(binIdx).spikeCount = sum(histc(sptimesNEV,iGapEndSetNev([1 gapIdx])-1 ));
            
            
%             stem(sptimes(1:50),ones(50,1),'k')
%             hold on
%             plot(gapTimes([1 1 gapIdx gapIdx])-1,[0 2 2 0],'r')
%             drawnow,pause(0.05)
%             xlim([0 max(sptimes(50))])
  
            
            % binary partitioning of the entire spike sequence
            %gapIdx = round(gapIdx*2/3);
            gapIdx = gapIdx - 1;
            %gapIdx = find(gapTimes <= gapTimes(gapIdx)/2,1,'last');
            binIdx = binIdx + 1;
            if (binIdx > nMaxBins)
                break
            end
        end

        %bins = fliplr(bins);    % first bin is biggest, 
                                % lets us literally narrow the field as we go
        bins = bins(randperm(length(bins)));
                                
        %if (min(abs(data.nev(fi).chans - ach/10)) < 0.1)
        if dbFlag
            % convert spike times from nev samples to data samples
            %sptimes = sptimesNEV(1:bins(end).spikeCount)*SR/single(nev.MetaTags.SampleRes);
            %samples = double(hist(sptimes,sptimes(end)-sptimes(1)));
        end
        if  dbFlag
            skip = false;       % set to true via keyboard to skip channel
            figure(3),clf
            %stem(1.1*samples,'b')
            %stem(sptimes,0.9*ones(size(sptimes)),'b')
            hold on
        end
        
        
        for ch = find(chanSet == 0)
            try
                tiSet = 1:length(data.stateHist.time);      % TODO prune if possible
                si = length(bins);
                matches = tiSet;
                
%                 if (abs(ach/10 - data.nev(fi).chans(ch)) < 0.1)
                if  dbFlag && (sum(chanSet > 0) || ~skip)
%                     try
%                         Tguess = data.nev(fi).Toffset(1);
%                     catch me
                        Tguess = mean(Tset(chanSet>0));
%                     end
                                    
                hData = [];     % reference to plot of spikes from data    
                hMatch = [];
                hScreen = [];
                end
                
                breakFlag = false;
                convergeCount = nConverge;
                while ~isempty(matches) %&& (convergeCount > 0))     % we expect to see matches @ same places
                    
                    % If all sum checks have been passed
                    if si == 0 ||( max(matches)-min(matches) < convergeWidth && convergeCount == 0)
                        iOffsetData = median(matches)-iBinOffsetData;
%                         if length(matches) > 2
%                             fprintf('Warning: %d matches,',length(matches));
%                             fprintf(' spread %d samples\n',max(matches)-min(matches));
%                         end
                        
                        fprintf('Counts %d iOffsetData %dm %.4fs\n',...
                            length(bins)-si, floor(iOffsetData/SR/60), mod(iOffsetData/SR,60));
                        
                        if (max(matches) - min(matches)) > nMaxSpread
                            fprintf('Too large a spread, can''t deal with this, moving on...\n')
                            break
                        end
                        %%%%%% Run brute-force error minimization to nail
                        %%%%%% exact offset time
                        iOffsetRange = matches-iBinOffsetData;
                        iOffsetRange = [min(iOffsetRange)-10 max(iOffsetRange)+5];
                        tOffsetRangeSec = iOffsetRange/SR;
                        tOffsetRangeSec = (tOffsetRangeSec(1)-1/SR):1/SRnev:(tOffsetRangeSec(end)+1/SR);
                        if ~isempty(chans) && (ach == chans(end,ch) && dbFlag)
                            Toffsets = NEV_FindTimeShift(data.stateHist.spikes(:,ch),SR,sptimesNEV,SRnev,tOffsetRangeSec,1);
                        else
                            Toffsets = NEV_FindTimeShift(data.stateHist.spikes(:,ch),SR,sptimesNEV,SRnev,tOffsetRangeSec,dbFlag);
                        end
                        if (Toffsets ~= -1)
                            iOffsetData = Toffsets*SR;
                            
                            if chanSet
                                % update search start pos to sooner of existing,
                                % or current offset minus 1 second tolerance
                                tiStart = min(tiStart, iOffsetData-SR);
                            else
                                tiStart = iOffsetData-SR;
                            end
                            
                            fprintf('Elec %d Unit %d <==> LabVIEW ch %d ',...
                                elec, unit, ch);
                            fprintf('Counts %d Toffset %dm %.4fs\n',...
                                length(bins)-si, floor(iOffsetData/SR/60), mod(iOffsetData/SR,60));
                            if (sum(isnan(Tset)))
                                fprintf('Getting nonsense results, bailing on this whole day\n')
                            end
                            chanSet(ch) = ach;
                            Tset(ch) = iOffsetData;
                        else
                            fprintf('Failed to find minimum\n')
                        end
                        breakFlag = true;
                        %firstPass = false;
                        if (~dbFlag)
                            break
                        end
                    end
                    
                    if (~breakFlag)
                        
                        % TODO: can prune down stateHist here
                        sIdx = min(matches):max(matches);
                        screenA = conv(data.stateHist.spikes(:,ch),ones(bins(si).periodCount,1),'valid');
                        screenB = conv(data.stateHist.spikes(:,ch),ones(bins(si).periodCount+1,1),'valid');
%                        screenC = conv(data.stateHist.spikes(:,ch),ones(bins(si).periodCount+2,1),'valid');
                        newMatches =  find(screenA(1:end-2) == bins(si).spikeCount...
                            | screenB(1:end-1) == bins(si).spikeCount); % ...
%                            | screenC == bins(si).spikeCount);
                        
                        matches = intersect(newMatches ,matches);
                        if (max(matches) - min(matches) < convergeWidth)
                            convergeCount = convergeCount -1;
                        end
                    end
                    
                   % if (abs(ach/10 - data.nev(fi).chans(ch)) < 0.1)
                    if dbFlag && ~skip && ~isempty(matches)
%                         breakFlag || ...
%                             (sum(chanSet) && dbFlag && ~ismember(ach,chanSet) && ~isempty(matches))...
%                             || ~skip
%                         figure(2)
% 
%                         plot(matches+round(TsplitStartSR),'.')
%                         hold on
%                         plot([0 length(matches)],[1 1]*data.nev(fi).Toffset(1),'r-')
                        
                        % get matches closest to mean of already-found
                        % channels

                        [~,bestMatchIdx] = sort(abs(matches + iBinOffsetData - Tguess),'ascend');
                        bestMatch = matches(bestMatchIdx(1:min(10,end)));
                        % matches are global references to data indices 
                        
                        if breakFlag
                            si = 1; % to plot last iteration 
                        end
                        
                        % everything in this plot uses current best guess Toffset
                        figure(3) % everything in real time relative to Data start
                        
              
                     if (ishandle(hData))
                            delete(hData)
                            delete(hMatch)
                            delete(hScreen)
                     end
                        pCount = bins(si).periodCount;
                        % the -1 here is so the stairs 'enclose' the blue
                        % NEV plotted samples, since the stairs behavior is
                        % ZOH and we want the instant value to represent
                        % the sum of spikes over the current time sample
                        hData = [];
%                         for m = reshape(bestMatch,1,numel(bestMatch))
%                             binCumSum = cumsum(data.stateHist.spikes(round(m + iBinOffsetData) +(0:pCount),ch));
%                             hData = [hData, stairs([zeros(round(iBinOffsetData)-1,1);...
%                             ],'g.-')];
%                             hold on
%                         end
                        
%%% Hmmmmmm
%                        hScreen = [ plot(screenA(round(bestMatch(1))+(0:bins(si).periodCount))-bins(si).spikeCount,'m.'),...
%                         plot(screenB(round(bestMatch(1))+(0:length(samples)))-bins(si).spikeCount,'m.')];%,...
                         %plot(screenC(round(bestMatch(1))+(0:length(samples)))-bins(si).spikeCount,'m.')];

                        
                        hold on
                        plot([0 0 1 1]*pCount+iBinOffsetData,[0 1 1 0]*bins(si).spikeCount,'r');
                        hold on
                        %stem(TsplitStartSR +data.stateHist.spikes(data.nev(fi).Toffset(1)+1+(0:splits(si).periodCount+5),ch),'k')
                       
                        hMatch = plot(bestMatch-bestMatch(1)+iBinOffsetData,zeros(size(bestMatch)),'rp','MarkerSize',10);
                        ylim([-5 5])
                        drawnow,pause(0.05)

                        
                    end
if (breakFlag)
                            break
                        end

                    si = si -1;
                end
                
               
                
            catch me
                me.getReport()
                keyboard
                
            end
            end
            
        
        
        
        

        
        end
%         if (toc > Tgiveup || Ntried > Ngiveup)
%             fprintf('**Moving on; have searched %d of %d electrode/unit combinations.\n',...
%                 find(ach == activeChans),length(activeChans));
%             break
%         end
        if ~firstPass
            break
        end
        
    end
    if (true)
        fprintf('\nSECOND PASS - brute force error minimization\n')
    end
    
    
%%%%%%%%% SECOND PASS - using error minimization over short guess windows
    if sum(Tset > 0) && secondPass
    for ach = activeChans
        elec = round(ach/10);
        unit = ach - 10*elec;
        
        
        
        % get indices of first Nspikes matching this channel
        spikeIdx = find(spikeData == ach, Nspikes);
        if (length(spikeIdx) < Nspikes && dbFlag)
            %fprintf('Warning: only found %d spikes on channel %2.1f\n',...
            %    length(spikeIdx), ach/10);
        end
        
        % retrieve the spike times from the nev data, uint32 -> float32
        sptimesNEV = single(nev.Data.Spikes.TimeStamp(spikeIdx));
        
        
        for ch = find(chanSet == 0)
            if isnan(tOffsetGuess)
                tGuessSec = mean(Tset(Tset > 0))/SR;
            else
                tGuessSec = tOffsetGuess;
            end
if (dbFlag >= 1)

            fprintf('Checking %d-%d -> %d\n',elec,unit,ch);
end
            tOffsetRangeSec = (tGuessSec-0.1):0.0005:(tGuessSec+0.1);
            if ~isempty(chans) && (ach == chans(end,ch) && dbFlag)
                Toffsets = NEV_FindTimeShift(data.stateHist.spikes(:,ch),SR,sptimesNEV,SRnev,tOffsetRangeSec,1);
            else
                Toffsets = NEV_FindTimeShift(data.stateHist.spikes(:,ch),SR,sptimesNEV,SRnev,tOffsetRangeSec,dbFlag);
            end
            
            if (Toffsets ~= -1)     % winner
                iOffsetData = Toffsets*SR;
                fprintf('Elec %d Unit %d <==> LabVIEW ch %d ',...
                            elec, unit, ch);
                        fprintf('Toffset %dm %.4fs\n',...
                            floor(iOffsetData/SR/60), mod(iOffsetData/SR,60));
                        chanSet(ch) = ach;
                        Tset(ch) = iOffsetData;
                        if sum(isnan(Tset))
                            fprintf('Getting nonsense times; bailing on whole day\n\n')

                            return
                        end
            end
        end

    end
    else
        fprintf('Skipping, no previous matches to work with!\n')
    end
    
        
    data.nev(fi).chans = chanSet/10 - (chanSet == 0);
    T = round(mean(Tset(Tset ~= 0)));
    tOffsetLast = T/SR;
    if isnan(T)
        data.nev(fi).Toffset = [0 0 0];
    else
        data.nev(fi).Toffset = [T min(Tset(Tset ~= 0))-T max(Tset(Tset ~= 0))-T];
    end
    data.nev(fi).nevfile = file;
    data.nev(fi).DurationSec = nev.MetaTags.DataDurationSec;
    fi = fi +1;
    chans = [chans; chanSet];
    Toffsets = [Toffsets; ]; % identified channels vote on time offset
    fnames = [fnames; file];
end


%%% Final step: figure out which mapping was used during each NEV%%%%%%%%%%
%
%
%
trials = metaData.trials{midx};
Toffset = cell2mat({data.nev.Toffset}');
Toffset = Toffset(:,1);

for fi = 1:length(data.nev)
    try
        data.nev(fi).ok = false;
        data.nev(fi).OverlapSec = 0;
    if find(data.nev(fi).chans)
    % figure out LabVIEW time indices corresponding to start & end of nev
    tiStart = data.nev(fi).Toffset(1);
    tiStop = tiStart + round(nev.MetaTags.DataDurationSec*SR);

    
%     si = find(data.nev(fi).Toffset(1)/SR < [trials.segInfo.tpause] & ...
%               data.nev(fi).Toffset(1)/SR + data.nev(fi).DurationSec > ...
%                                            [trials.segInfo.trun]);
    % This method ensures 51% of the nev overlaps with the data segment
    % avoiding duplicate matches between data segments and nevs
    si = find(min([trials.segInfo.tpause] - data.nev(fi).Toffset(1)/SR,...
                  data.nev(fi).Toffset(1)/SR + data.nev(fi).DurationSec - ...
                    [trials.segInfo.trun]) >= ...
                        0.51*data.nev(fi).DurationSec);
                 % 0.50*([trials.segInfo.tpause] - [trials.segInfo.trun]));

    if ~isempty(si)
        data.nev(fi).segInfo = trials.segInfo(si);
        data.nev(fi).map = trials.segInfo(si).map;          % ASSUMPTION
        data.nev(fi).OverlapSec = ...       % overlap between trial block & NEV file
            min(trials.segInfo(si).tpause, data.nev(fi).Toffset(1)/60+data.nev(fi).DurationSec) - ...
            max(trials.segInfo(si).trun, data.nev(fi).Toffset(1)/60);
        
        % Set 'ok' flag if no significant timing issues & all channels ID'd
        if (diff(data.nev(fi).Toffset(2:3)) <= ToffsetSpread) && ...
                ~sum(data.nev(fi).chans(logical(data.nev(fi).segInfo.enabled(1:4))) == -1)
            data.nev(fi).ok = true;
        end
            
    else
        data.nev(fi).segInfo = [];
    end
    
    if isempty(data.nev(fi).segInfo)
        
        if (~isnan(tiStart) && ~isnan(tiStop))
            % convert to time in seconds during LabVIEW recording
            tStart = data.stateHist.time(max(tiStart,1));
            tStop = data.stateHist.time(min(tiStop, Ndata));
            
            % find trials that ended within that interval
            idxSet = trials.time > tStart & trials.time < tStop;
            
            if (sum(idxSet) < 10)
                data.nev(fi).map = '(waitfile)';    % no trials during this nev (e.g. wait file)
            else
                % find mappings used by those trials, hopefully just one!
                maps = unique([trials.Map(idxSet)]);
                maps = setdiff(maps,{''});
                
                if (length(maps) > 1)
                    fprintf('Multiple trial types during %s\n', data.nev(fi).nevfile)
                    for m = maps
                    fprintf('%d of %s\n',sum(ismember(trials.Map(idxSet),m{1})),m{1})
                    end
                    data.nev(fi).map = 'Multiple';
                    
                else if (isempty(maps))
                        data.nev(fi).map = '';
                    else
                        data.nev(fi).map = maps{1};
                    end
                end
                
            end
        else
            data.nev(fi).map = 'Error';
            if (dbFlag >= 2)
                keyboard
            end
        end
    end
    
    else
        data.nev(fi).map = '';
        data.nev(fi).ok = false;
    end
  catch me
            me.getReport()
            keyboard
        end  
end

if (dbFlag)
    fprintf('\nSummary:\n')
    Toffsets = cell2mat({data.nev.Toffset}');
    num2str([cell2mat({data.nev.chans}') floor(Toffsets(:,1)/SR/60) mod(Toffsets(:,1)/SR,60)],'%2.1f\t%2.1f\t%2.1f\t%2.1f\t%d:%2.3f\n')
end

    Toff = cell2mat({data.nev.Toffset}');   % recalc for new results     
    
    if (sum(sum(cell2mat({data.nev.chans}) ~= 0)) >= nChanExisting || sum(diff(Toff(:,2:3)')) < skew/2)
    data.nevUpdateNum = now;
    fprintf('Saving channel & timing info to %s.mat\n',data.fname);
    save([matpath filesep data.fname '.mat'],'data');
else
    fprintf('Dumping results, existing mapping looks better...\n\n')
    end

    catch me
        me.getReport
        keyboard
end
