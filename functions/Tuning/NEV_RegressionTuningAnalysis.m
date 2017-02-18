% NEV_RegressionTuningAnalysis.m
% by C. Matlack 2013-12-10
%
% USAGE: fits = NEV_RegressionTuningAnalysis(nev, nsx, delaySec, dbFlag)
%
% INPUTS:   nev         nev data object (spike recordings)
%           nsx         nsx data object (analog recordings)
%           delaySec    delay from torque to neural firing in seconds
%
%   If no args given, will prompt for selection via GUI.



function [fits, cumR2score, tuningmatrix] = NEV_RegressionTuningAnalysis(nev, nsx, delaySec, wFilter, dbFlag)

binWidthSec = 1/25;
minRateHz = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load NEV & NS3 files if not already provided
if (nargin <=4)
    % These values were optimal for one tested data set of neurons above
    % 10Hz
    wFilter = 0.3;      % Std dev of gaussian filter to apply
    delaySec = -0.075;
    dbFlag = 1;
end

if (nargin == 0)
    
    % Load NEV file
    [fname, fpath] = uigetfile('.nev','Select NEV file');
    fFullPath = strcat(fpath,fname);
    nev = openNEV(fFullPath,'nosave');
    
    % Load NS3 file, assumed to be at same location
    fprintf('Please wait, loading NSx file...\n')
    nsx = openNSx(strcat(fFullPath(1:end-3),'ns3'),'read');
    fprintf('done.\n');
   
    dbFlag = 1;
end




% Things I don't like to type over and over
TendSec = nev.MetaTags.DataDurationSec;
SRnev = double(nev.MetaTags.SampleRes);
SRnsx = double(nsx.MetaTags.SamplingFreq);
Gauss = @(s,w) 1/sqrt(2*pi)/w*exp(-s.^2/2/w^2);   % define Gaussian kernel



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Extract and bin neuron spike information
% Screen for neurons
fprintf('Finding neurons firing at least %.1f Hz...\n',minRateHz)
catalog = NEV_ScreenForMinRate(nev,minRateHz);
Nunits = length(catalog);

% Use histogram method to group spikes into bins of specified width
% sampled = zeros(floor(TendSec/binWidthSec)+1,Nunits);
% for i = 1:Nunits
%     elec = catalog(i,1);
%     unit = catalog(i,2);
%     
%     spIdx = find((nev.Data.Spikes.Electrode == uint16(elec)) & ...
%             (nev.Data.Spikes.Unit == uint8(unit)));
%     
%     sampled(:,i) = histc(double(nev.Data.Spikes.TimeStamp(spIdx))/SRnev,0:binWidthSec:TendSec);
% end
nev = NEV_SptimesToSamples(nev,SRnsx, catalog, 0, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter sampled spike trains
fprintf('Smoothing spike trains (this takes a minute)... ');
tt = 0:1/SRnsx:4*wFilter;
tt = [-fliplr(tt) tt(2:end)];
gaussKern = Gauss(tt,wFilter);
gaussKern = gaussKern/sum(gaussKern);
smoothed = conv2(nev.Sampled.Spikes',gaussKern,'same')';

if (dbFlag > 1)
    figure(1),clf
    plot(1:SRnsx,nev.Sampled.Spikes(1:SRnsx,1),1:SRnsx,smoothed(1:SRnsx,1))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use Matlab's resampling function to get down-sampled & smoothed copy of
% torque signals at matching sample rate
nChannels = size(nsx.Data,1);

% pull out torque channels, make col vectors, flip sign of FE for
% coordinate system consistency
torqueRaw = [nsx.Data(nChannels-1,:); -nsx.Data(nChannels-2,:)]'; 

% normalize to full-scale value of 1 (raw is int16)
torqueRescaled = double(torqueRaw)/(2^15);

% subtract mean
torqueCentered = torqueRescaled - repmat(mean(torqueRescaled),length(torqueRescaled),1);

% Gaussian filter
torqueFiltered = conv2(torqueCentered',gaussKern,'same')';

%%% Both signals should be same length at this point. Trim as necessary.
nSamp = min(size(torqueFiltered,1),size(smoothed,1));
torqueFiltered = torqueFiltered(1:nSamp,:);
smoothed = smoothed(1:nSamp,:);
fprintf('done.\n')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply delay, rounding off to even number of samples
delaySamples = round(delaySec*SRnsx);
if (delaySamples > 0)
    smoothed = smoothed(1+delaySamples:end,:);
    torqueFiltered = torqueFiltered(1:end-delaySamples,:);
else if (delaySamples < 0)
    smoothed = smoothed(1:end+delaySamples,:);
    torqueFiltered = torqueFiltered(1-delaySamples:end,:);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Resample at low sampling rate to avoid overfitting
%
fprintf('Resampling at %.1f Hz...',1/binWidthSec);
% resample at SR/SamplingFreq * SamplingFreq = SR
torqueResamp = resample(torqueFiltered, 1/binWidthSec, SRnsx);
smoothed = resample(smoothed,1/binWidthSec,SRnsx);
fprintf('done.\n');

% Make sure all of the above looks clean, and timing isn't skewed
if (dbFlag >= 2)
    figure(1),clf
    plot((1:length(torqueRescaled))*1/SRnsx, torqueRescaled,'k')
    hold on
    plot((0:length(torqueResamp)-1)*binWidthSec, torqueResamp,'b')
    plot((0:length(torqueCentered)-1)*1/SRnsx, torqueCentered,'g')
    plot((0:length(torqueFiltered)-1)*1/SRnsx, torqueFiltered,'m')
    xlim([0 10])
    pause(.5)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we can run the model fits that predict neural firing based on torque,
% with both signals in binWidthSec bins.
fprintf('Fittng... ');
% Prune signals: knock off time 0 value from torque.
%  Hopefully they're the same length!
torque = torqueResamp;
sampled = smoothed;

%%% Model fit:  Firing rate = baseline + 2D torque
TuningModel = fittype(@(a,ru,fe,Tru,Tfe) a + ru*Tru + fe*Tfe,...
                    'independent',{'Tru','Tfe'});
cumR2score = 0;        
if (dbFlag >= 1)
    fprintf('\nChan\tFE\tRU\tR^2\n')
end

try
    tuningmatrix = [];
for i=1:Nunits
    [fitObj, gof, output] = fit(torque,sampled(:,i),TuningModel,'StartPoint',[1 1 1]);
    fits(i).gof = gof;
    fits(i).fo = fitObj;
    fits(i).neuron = catalog(i,:);
    cumR2score = cumR2score + gof.rsquare;

    if (dbFlag >= 1)
        fprintf('%d-%d\t%.1f\t%.1f\t%.2f\n',...
            catalog(i,1),catalog(i,2),fitObj.fe,fitObj.ru,gof.rsquare);
            tuningmatrix = [tuningmatrix; catalog(i,1) catalog(i,2) fitObj.fe fitObj.ru gof.rsquare];
    else
            fprintf('.');
    end
    
    if (dbFlag >= 2)
        figure(1),clf
        plot(torque)
        hold on
        plot(sampled(:,i),'r')
        xlim([0 100])
        pause(.25)
    end
        
end
fprintf('done.\n')

catch me
    me.getReport
    keyboard
end

%%%

if (dbFlag)
    figure(1),clf
    labels = {};
    for fi = 1:length(fits)
        rsq = fits(fi).gof.rsquare;
        fe = fits(fi).fo.fe;
        ru = fits(fi).fo.ru;
        subplot(10,1,1:7)
        plot([0 fe], [0 ru],'LineWidth',round(rsq*10)+1)
        hold all
        labels = {labels{:},num2str(catalog(fi,1)+catalog(fi,2)/10)};
        text(fe,ru,labels{fi})
        
        subplot(10,1,9:10)
        stem(fi,fits(fi).fo.a,'LineWidth',round(rsq*10)+1)
        hold all
    end
    set(gca,'XtickL',labels(:))
    ylabel('Baseline rates')
    xlabel('Channel + Sort')
    
    subplot(10,1,1:7)
    title('Population vector mapping')
    xlabel('Flexion/Extension tuning depth')
    ylabel('Radial/Ulnar tuning depth')
end

end


