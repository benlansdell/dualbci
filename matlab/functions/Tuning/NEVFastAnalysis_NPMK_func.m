function tuningmatrix = NEVFastAnalysis_NPMK_func(fFullPath)
    % By C. Matlack 2012-12-15
%
% Work in progress
% Goal: cross-correlation analysis of top-firing neurons

% DEPENDS:      SptimesToSpikes()
%               NPMK library from BlackRock

%clear all

%% PARAMETERS: 
SR = 50;        % sampling rate for sampled reconstruction of spike trains (Hz)
minRate = 5;   % Minimum firing rate to consider cell (Hz)
wKern = 0.1;    % width of Gaussian kernel to use for filtering
clip = 4;       % factor of w at which to truncate Gaussian kernel
maxLag = 90;   % maximum time shift to check for in cross-correlation
maxPeak = 1; % maximum time shift for admissible peak in covariance

plotFlag = false;

maxLag = round(maxLag*SR)/SR;   % make maxLag even number of samples

%% DEFINITIONS
Gauss = @(s,w) 1/sqrt(2*pi)/w*exp(-s.^2/2/w^2);   % define Gaussian kernel

%% Load NEV file
%[fname, fpath] = uigetfile('.nev','Select NEV file');
%
wd = pwd;

%fFullPath = strcat(fpath,fname)
nev = openNEV(fFullPath,'report','nosave','nomat');
fprintf('Please wait, loading NSx file...\n')
nsx = openNSx(strcat(fFullPath(1:end-3),'ns3'),'read','report');


%% TO-DO: Quickly identify good neurons by firing rate
fprintf('Finding good neurons...\n')

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
fastCellsIdx = find(spikeAvg >= minRate);
Nunits = length(fastCellsIdx);

nev.catalog = [channels(fastCellsIdx)/10, mod(channels(fastCellsIdx),10)];

fprintf('Elec\tSort\tAvg rate\n')
disp([channels(fastCellsIdx)/10, mod(channels(fastCellsIdx),10), spikeAvg(fastCellsIdx)])

%% OPTIONAL TO-DO: Calculate optimal Gaussian kernel width for each electrode



%% Convert to sampled data THIS PART IS SLOW

fprintf('Preallocating sampled data array....\n')

nev.sampleTimes = (0:1/SR:(nev.MetaTags.DataDurationSec+1))';
nev.sampledSpikes = repmat(0*nev.sampleTimes,1, Nunits);
nev.sampleRate = SR;

fprintf('Converting to sampled spikes...\n');

for i = 1:Nunits
        elec = channels(fastCellsIdx(i))/10;
        unit = mod(channels(fastCellsIdx(i)),10);
        fprintf('Elec %d, unit %d, %5.0f spikes...\n',...
            elec,unit,spikeAvg((fastCellsIdx(i)))*nev.MetaTags.DataDurationSec);
        nev.sampledSpikes(:,i) = SptimesToSpikes(...
        single(nev.Data.Spikes.TimeStamp((nev.Data.Spikes.Electrode == elec) & (nev.Data.Spikes.Unit == unit)))/...
            nev.MetaTags.TimeRes,...
        nev.sampleTimes);
end
fprintf('done\n');

%% Convolve with Gaussian kernel of fixed width

fprintf('Smoothing & calculating stats...');
gaussKern = Gauss(-clip*wKern:1/SR:clip*wKern,wKern);

nev.smoothedSpikes = zeros(size(nev.sampledSpikes));
nev.modulation = zeros(Nunits, 2);

for i = 1:size(nev.smoothedSpikes,2)
    nev.smoothedSpikes(:,i) = conv(nev.sampledSpikes(:,i),gaussKern,'same');
    %nev.smoothedSpikes(:,i) = nev.smoothedSpikes(:,i) - mean(nev.smoothedSpikes(:,i));
    
    % Calculate basline firing rate and modulation depth as mean & moment of rate
    nev.modulation(i,:) = [mean(nev.smoothedSpikes(:,i)), std(nev.smoothedSpikes(:,i))];
end

fprintf('done\n');

%%
%fprintf('\tElec\tSort  Baseline Std\n')
[junk, sortIndex] = sort(nev.modulation(:,2),1,'descend');
%nev.catalog = nev.catalog(sortIndex,:);
%nev.modulation = nev.modulation(sortIndex,:);

%disp([sortIndex, nev.catalog(sortIndex,:), nev.modulation(sortIndex,:)])
%disp([nev.catalog, nev.modulation])

%% Load torque data from NSx 
% Down-sample to 60Hz
nChannels = size(nsx.Data,1);
nsx.torque = nsx.Data((nChannels-2):(nChannels-1),:)';  % pull out torque channels, make col vectors
% resample at SR/SamplingFreq * SamplingFreq = SR
nsx.torque = resample(double(nsx.torque), SR, nsx.MetaTags.SamplingFreq);
% run torque signal through Gaussian filters
nsx.torque(:,1) = conv(nsx.torque(:,1),gaussKern,'same');
nsx.torque(:,2) = conv(nsx.torque(:,2),gaussKern,'same');
nsx.jerk = diff(nsx.torque);

Nsamp = min(length(nsx.jerk),length(nev.smoothedSpikes));




%% Torque covariance analysis
fprintf('Finding torque covariances...\n')
%figure(2)
tt = -maxLag:1/SR:maxLag;
peakRange = find(tt > -maxPeak & tt < maxPeak);    % area to search for peak

fprintf('Sampled NEV has %d samples, downsampled NSX has %d, using %d.\n\n',...
    length(nev.smoothedSpikes),length(nsx.torque),Nsamp);

peakFE = zeros(Nunits,1);
peakRU = peakFE;
lagFE = peakFE;
lagRU = lagFE;
stdFE = peakFE;
stdRU = peakRU;

for i = 1:Nunits
    
    %subplot(Nunits,2,i*2-1)
    covFE = xcov(nev.smoothedSpikes(1:Nsamp,i), single(nsx.torque(1:Nsamp,1)),SR*maxLag,'unbiased');
    % normalize against spikes auto-covariance
    covFE = covFE / sqrt(xcov(nev.smoothedSpikes(1:Nsamp,i),0));
    covFE = covFE / sqrt(xcov(single(nsx.torque(1:Nsamp,1)),0));
    peakFE(i) = covFE(abs(covFE) == max(abs(covFE(peakRange))));
    lagFE(i) = tt((covFE == peakFE(i)) | (covFE == -peakFE(i)));
    stdFE(i) = std(covFE);
    avg = mean(covFE);
    peakFE = peakFE - avg;
    
    if (plotFlag)
    figure
    subplot(1,3,1)
    plot(tt,covFE,[-maxLag maxLag],[avg avg],...
                  [lagFE(i) lagFE(i)],[avg avg+peakFE(i)],...
                  [-maxLag maxLag],avg+[stdFE(i) stdFE(i)]*sign(peakFE(i)));
    title(sprintf('FE, Ch%d, Elec%d-%d, Score %2.1f',i,nev.catalog(i,1),nev.catalog(i,2),peakFE(i)/stdFE(i)))
    xlim([-maxLag maxLag])
    end
    
    covRU = xcov(nev.smoothedSpikes(1:Nsamp,i), single(nsx.torque(1:Nsamp,2)),SR*maxLag,'unbiased');
    covRU = covRU / sqrt(xcov(nev.smoothedSpikes(1:Nsamp,i),0));
    covRU = covRU / sqrt(xcov(single(nsx.torque(1:Nsamp,2)),0));
    peakRU(i) = covRU(abs(covRU) == max(abs(covRU(peakRange))));
    stdRU(i) = std(covRU);
    lagRU(i) = tt((covRU == peakRU(i)) | (covRU == -peakRU(i)));
    if (plotFlag)
    subplot(1,3,2)
    plot(tt,covRU,[-maxLag maxLag],[avg avg],...
                  [lagRU(i) lagRU(i)],[avg avg+peakRU(i)],...
                  [-maxLag maxLag],avg+[stdRU(i) stdRU(i)]*sign(peakRU(i)));
    
    title(sprintf('RU Cov, Ch %d, Elec %d, Sort %d, Score %2.1f',i,nev.catalog(i,1),nev.catalog(i,2),peakRU(i)/stdRU(i)))
    xlim([-maxLag maxLag])
    
    subplot(1,3,3)
    scatter(-covFE,covRU,1,tt);
    axis equal
    axis square
    %fprintf('Press enter to continue...\n')
    %pause %(0.5)
    end
end

if (length(nev.catalog) >= 7)
fprintf('\n********************************************************\n');
[lscores,srtIdx] = sort(peakFE./stdFE,'ascend');
fprintf('\n Best RIGHT (0) neurons and scores: \n')
fprintf('%3.1f\t',((single(nev.catalog(srtIdx(1:7),1:2))*[1; 0.1])' ))
fprintf('\n')
fprintf('%3.1f\t\t',-lscores(1:7))
fprintf('\n\n Best LEFT (180) neurons and scores: \n')
[rscores,srtIdx] = sort(peakFE./stdFE,'descend');
fprintf('%3.1f\t',((single(nev.catalog(srtIdx(1:7),1:2))*[1; 0.1])' ))
fprintf('\n')
fprintf('%3.1f\t\t',rscores(1:7))
fprintf('\n\n If doing Push/Push, go with ')
if (sum(-lscores(1:4)) > sum(rscores(1:4)))
    fprintf('LEFT (180).\n');
else
    fprintf('RIGHT (0).\n');
end

end

%%
fprintf('\n********************************************************\n');
fprintf('Detailed neuron information (pos score good left, neg good right):\n')
fprintf('Elec\tSort\tAvg\t\tStd\t\tFE Score\tRU Score\n')
[~,srtIdx] = sort(nev.modulation(:,2),'descend');
fprintf('%d\t\t%d\t\t%g\t\t%g\t\t%2.1f\t%2.1f\n',[single(nev.catalog(srtIdx,:)), round(nev.modulation(srtIdx,:)),peakFE(srtIdx)./stdFE(srtIdx),peakRU(srtIdx)./stdRU(srtIdx)]')
fprintf('\n********************************************************\n');
tuningmatrix = [single(nev.catalog(srtIdx,:)), round(nev.modulation(srtIdx,:)),peakFE(srtIdx)./stdFE(srtIdx),peakRU(srtIdx)./stdRU(srtIdx)];

return


%%% Tuning function analysis - MEMORY KiLLER
%fprintf('\n\nFinding OLS torque tuning model\n')
%Nsamp = min(Nsamp, 60*60*3);
%
%X = [ones(Nsamp,1), nsx.torque(1:Nsamp,:)-repmat(mean(nsx.torque,1), Nsamp,1), nsx.jerk(1:Nsamp,:)]; 
%regression = (X'*X)\X';
%P = X*regression;
%nev.tuningModel = zeros(Nunits,6);
%delete X
%% tuning model: b, cx, cy, cdx, cdy, R^2
%onevec = ones(Nsamp, 1);
%L = eye(Nsamp) - onevec*onevec'/Nsamp;
%for i=1:Nunits
%    y = nev.smoothedSpikes(1:Nsamp,i);
%    nev.tuningModel(i,1:5) = regression*y;
%    nev.tuningModel(i,6) = 1 - (y'*L*P*y)/(y'*L*y);
%end
%
%fprintf('Elec\tSort\tAvg\t\tStd\tBaseline\t Cx\t\t\tCy\t\t\tR^2\n')
%fprintf('%d\t\t%d\t\t%g\t\t%g\t\t%2.0f\t\t%1.4f\t\t%1.4f\t\t%1.4f\t\t%1.4f\t\t%1.3f\n',[single(nev.catalog), round(nev.modulation), nev.tuningModel]')
%%disp([single(nev.catalog), nev.modulation, nev.tuningModel(:,1),nev.tuningModel(:,2:end)*1000])
%
%fprintf('\n\nBest LEFT neurons in decreasing order:\n')
%[~,srtIdx] = sort(nev.tuningModel(:,2),'ascend');
%fprintf('%3.1f\t',((single(nev.catalog(srtIdx(1:5),1:2))*[1; 0.1])' ))
%fprintf('<----\n\nBest RIGHT neurons in decreasing order:\n')
%[~,srtIdx] = sort(nev.tuningModel(:,2),'descend');
%fprintf('%3.1f\t',((single(nev.catalog(srtIdx(1:5),1:2))*[1; 0.1])' ))
%fprintf('---->\n')
%
%%% Pairwise cross-correlation analysis
%
%figure(1)
%tt = -maxLag:1/SR:maxLag;
%
%Nsamp = length(nev.smoothedSpikes);
%
%peakCov = zeros(Nunits);  % will only use upper triangle of this
%peakLag = peakCov;         % lag corresponding to max covariance
%
%fprintf('Finding auto-variances\n')
%% first find auto-correlations for normalization purposes
%
%for i = 1:Nunits
%    [cov, idx] = xcov(nev.smoothedSpikes(:,i),SR*maxLag,'unbiased');
%    peakCov(i,i) = max(abs(cov));
%    peakLag(i,i) = idx((cov == peakCov(i,i)) | (cov == -peakCov(i,i)))/SR;
%    
%end
%%%
%fprintf('Finding cross-covariances\n')
%for i = 1:Nunits
%    for j = i:Nunits
%        
%        [xcr, idx] = xcov(nev.smoothedSpikes(:,i),...
%                    nev.smoothedSpikes(:,j),...
%                    SR*maxLag,'unbiased');        
%        
%        norm = sqrt(peakCov(i,i))*sqrt(peakCov(j,j));
%        %[i,j,norm]
%        % normalize against auto-variance of each neuron
%        xcr = xcr / norm;
%        % This should make every auto-covariance maximum = 1
%        
%        xcrPerm = xcov(nev.smoothedSpikes(randperm(Nsamp),i),...
%                        nev.smoothedSpikes(randperm(Nsamp),j),...
%                        SR*maxLag,'unbiased');
%        xcrPerm = xcrPerm/norm;
%        
%        if (i ~= j)
%            peakCov(i,j) = max(abs(xcr));
%        
%            if (isempty(find(xcr == peakCov(i,j), 1)))
%                peakCov(i,j) = -peakCov(i,j);
%            end
%            
%            peakLag(i,j) = idx(xcr == peakCov(i,j))/SR;
%        end
%                    
%        %subplot(neurons,neurons,(i-1)*neurons + j)
%        %plot(idx/SR,xcr,'b-',...
%        %     idx/SR,repmat(mean(xcrPerm) - 2*moment(xcrPerm,2),1,length(idx)),'g-.',...
%        %     idx/SR,repmat(mean(xcrPerm) + 2*moment(xcrPerm,2),1,length(idx)),'g-.');
%        % hold on
%        % stem(peakLag(i,j),peakCov(i,j))
%         
%        %title(sprintf('(%d, %d)',nev.catalog(i),nev.catalog(j)))
%        %ylim([-.5 .5])
%    end
%end
%
%
%for i = 1:Nunits
%    peakCov(i,i) = 1;
%end
%%%
%
%HeatMap(flipdim(peakCov,1),'Colormap',redbluecmap)
%title('Covariance strength')
%HeatMap(flipdim(peakLag,1),'Colormap',redbluecmap)
%title('Lag at peak covariance')
%
%fprintf('Lags:\n')
%
%%%
%sorted = sort(reshape(abs(peakCov),1,Nunits^2),'descend');
%sorted = sorted((Nunits+1):end);    % clear out 1 values of auto-variance
%sorted = sorted(sorted ~= 0);        % clear out 0's in lower triangle
%

% TO-DO: Run cross-correlation on all pairs, compare to time-shifted
% average