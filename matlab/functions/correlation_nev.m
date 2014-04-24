function correlation_nev(nevfile, fn_out, threshold, binsize)
	%correlation_nev	Function to fit the following model to spike and torque data:
	%			lambda_i(t) = \lambda_0 + \sum_j^N k_j^1 . x_i^1(t+\tau+jh) + \sum_j^N k_j^2 . x_i^2(t+\tau+jh)
	%		That is, it fits a linear filter to the torque data. In the above formula, time-step size h
	%		is given by parameter binsize; max value of N is given by kernellength; \tau is given by offset.
	%		Program fits all models with kernels of length between 1 and N, meaning that for each single-unit N models
	%		will be fit.
	%
	% 		Input:
	%			nevfile = file to process
	%			fn_out = base name for output plots
	%			threshold = (optional, default = 5) threshold firing rate below which unit is ignored
	%			binsize = (optional, default = 0.05) time bin size for spikes (torque resmapled to this rate0)
	%		
	%		Output:
	%			(none) produces plots of the filter for each single-unit channel, along with summary of quality of fits for 
	%			each channel
	%
	%		Test code:
	%			nevfile = './testdata/20130117SpankyUtah001.nev';
	%			threshold = 5;
	%			fn = './worksheets/tuning/plots/corr_20130117SpankyUtah001';
	%			correlation_nev(nevfile, fn, threshold);
	
	%Optional arguments
	if (nargin < 3)	threshold = 5; end
	if (nargin < 4) binsize = 0.05; end
	%Set this to above 0 if want debug info/plots
	verbosity = 1;
	%Total number of possible units recorded from
	nE = 128;
	nunits = 5; 
	nU = nE*nunits;
	%Offset to apply
	offset = 0;
	%Std of gaussian filter to apply
	sigma = 5;
	%Size of gaussian filter to apply
	sz = 30;
	samplerate = 1/binsize;
	%Max lag for correlations
	maxlag = 90;
	maxlag = round(maxlag*samplerate)/samplerate;
	maxpeak = 1;
	%Make sure we can perform the sample rate conversion easily
	assert(rem(samplerate,1) == 0, 'Select a binsize corresponding to an integer sample rate.');
	ns3file = [nevfile(1:end-3) 'ns3'];

	%%%%%%%%%%%%%%%%%%%%%%
	%Process spiking data%
	%%%%%%%%%%%%%%%%%%%%%%
	NEV = openNEV(nevfile);
	%Find the duration and sample rate of the nev file recording
	nevsamplerate = NEV.MetaTags.TimeRes;
	dur = NEV.MetaTags.DataDuration/nevsamplerate;
	%Convert spike times into array of binned spikes, one for each spike sorted channel
	spiketimes = double(NEV.Data.Spikes.TimeStamp)/nevsamplerate;
    elecs = cell(1,nU);
    spikemuas = struct('times', elecs);
    unitnames = cell(1,nU);
    averate = zeros(1,nU);
    for idx=1:nU
        spikemuas(idx).times = [0];    
    end
    for i=1:length(spiketimes)
       	E = NEV.Data.Spikes.Electrode(i);
       	unit = NEV.Data.Spikes.Unit(i);
       	U = single((E-1)*nunits)+single(unit)+1;
       	spikemuas(U).times = [spikemuas(U).times; spiketimes(i)];
       	unitnames{U} = [num2str(E) '.' num2str(unit)];
    end
    %Check which channels are doing stuff
    for idx=1:nU
    	averate(idx) = (length(spikemuas(idx).times)-1)/dur;
    	if (verbosity) display(['Electrode.Unit: ' unitnames{idx} ' Spike count: ' num2str(length(spikemuas(idx).times)-1) ' Mean firing rate (Hz): ' num2str(averate(idx))]); end
    end
    %Set a threshold firing rate, below which we ignore that unit
    abovethresh = averate > threshold;
    %Update nU
    nU = sum(abovethresh);
    unitnames = unitnames(abovethresh);
    spikemuas = spikemuas(abovethresh);
    averate = averate(abovethresh);
   	%Bin spikes (chronux function)
    binnedspikes = binspikes(spikemuas, samplerate);
    %From this apply gaussian filter to spike train for each electrode
    x = linspace(-sz/2, sz/2, sz);
    gaussFilter = exp(-x.^2/(2*sigma^2));
    gaussFilter = gaussFilter/sum(gaussFilter);
    for idx=1:nU
        gf = conv(binnedspikes(:,idx), gaussFilter, 'same');
        rates(:,idx)=gf*samplerate;
    end

	%%%%%%%%%%%%%%%%%%%%%
	%Process torque data%
	%%%%%%%%%%%%%%%%%%%%%
	NS3 = openNSx(ns3file, 'read', 'c:138:139');
	nsxtorque = double(NS3.Data);
	nsxsamplerate = double(NS3.MetaTags.SamplingFreq);
	%Switch sign of FE axis for coordinate consistency
	nsxtorque(2,:)=-nsxtorque(2,:);
    for j=1:2
        %Scale from uint16 value to proportion
        nsxtorque(j,:) = nsxtorque(j,:)/(2^15);
        %Subtract mean
        nsxtorque(j,:) = nsxtorque(j,:)-mean(nsxtorque(j,:));
        %Smooth
        nsxtorque(j,:) = conv(nsxtorque(j,:),gaussFilter,'same');
    end
	%Resample at rate of binsize
	torque = resample(nsxtorque', samplerate, nsxsamplerate);		
	%Check they're the same length, and trim
	nsamp = min(size(torque,1), size(rates,1));
	torque=torque(1:nsamp,:);
	rates = rates(1:nsamp,:);
	%Apply offset to data
	delaysamples = round(offset*samplerate);
	if (delaysamples > 0)
    	rates = rates(1+delaysamples:end,:);
    	torque = torque(1:end-delaysamples,:);
	elseif (delaysamples < 0)
	    rates = rates(1:end+delaysamples,:);
    	torque = torque(1-delaysamples:end,:);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Cross- and auto-correlation%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	tt = -maxlag:binsize:maxlag;
	peakrange = find(tt > -maxpeak & tt < maxpeak);    % area to search for peak
    peakRU = zeros(1,nU);
    peakFE = zeros(1,nU);
    lagRU = zeros(1,nU);
    lagFE = zeros(1,nU);
    stdRU = zeros(1,nU);
    stdFE = zeros(1,nU);
    scoreRU = zeros(1,nU);
    scoreFE = zeros(1,nU);
    for i=1:nU
    	%Compute cross correlation
	    covFE = xcov(rates(:,i), torque(:,1),samplerate*maxlag,'unbiased');
    	% normalize against spikes auto-covariance
    	covFE = covFE / sqrt(xcov(rates(:,i),0));
    	covFE = covFE / sqrt(xcov(torque(:,1),0));
    	peakFE(i) = covFE(abs(covFE) == max(abs(covFE(peakrange))));
	    stdFE(i) = std(covFE);
    	lagFE(i) = tt((covFE == peakFE(i)) | (covFE == -peakFE(i)));
    	avg = mean(covFE);
    	peakFE(i) = peakFE(i) - avg;
    	scoreFE(i) = peakFE(i)/stdFE(i);
	    covRU = xcov(rates(:,i), torque(:,2),samplerate*maxlag,'unbiased');
    	% normalize against spikes auto-covariance
    	covRU = covRU / sqrt(xcov(rates(:,i),0));
    	covRU = covRU / sqrt(xcov(torque(:,2),0));
    	peakRU(i) = covRU(abs(covRU) == max(abs(covRU(peakrange))));
    	stdRU(i) = std(covRU);
    	lagRU(i) = tt((covRU == peakRU(i)) | (covRU == -peakRU(i)));
    	avg = mean(covRU);
    	peakRU(i) = peakRU(i) - avg;
    	scoreRU(i) = peakRU(i)/stdRU(i);

    	%Plot cross- and auto-correlations
    	figure
    	subplot(2,2,1);
    	plot(tt,covFE,[-maxlag maxlag],[avg avg],...
               [lagFE(i) lagFE(i)],[avg avg+peakFE(i)],...
                  [-maxlag maxlag],avg+[stdFE(i) stdFE(i)]*sign(peakFE(i)));
		title(['FE Unit: ' unitnames{i} ' Score: ' num2str(scoreFE(i))])
		xlim([-maxlag maxlag])
    	subplot(2,2,2);
    	plot(tt,covRU,[-maxlag maxlag],[avg avg],...
               [lagRU(i) lagRU(i)],[avg avg+peakRU(i)],...
                  [-maxlag maxlag],avg+[stdRU(i) stdRU(i)]*sign(peakRU(i)));
		title(['RU Score: ' num2str(scoreRU(i))])
		xlim([-maxlag maxlag])
    	subplot(2,2,3);
    	plot(tt,covFE,[-maxlag maxlag],[avg avg],...
               [lagFE(i) lagFE(i)],[avg avg+peakFE(i)],...
                  [-maxlag maxlag],avg+[stdFE(i) stdFE(i)]*sign(peakFE(i)));
		xlim([-maxpeak maxpeak])
    	subplot(2,2,4);
    	plot(tt,covRU,[-maxlag maxlag],[avg avg],...
               [lagRU(i) lagRU(i)],[avg avg+peakRU(i)],...
                  [-maxlag maxlag],avg+[stdRU(i) stdRU(i)]*sign(peakRU(i)));
		xlim([-maxpeak maxpeak])
		saveplot(gcf, [fn_out '_unit_' unitnames{i} '_cross_maxscore_' num2str(max(scoreFE(i), scoreRU(i))) '.eps'], 'eps', [6 6]);
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Heat map of tau and xcov for all units%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	figure
	gridx = floor(sqrt(nU));
	if (gridx < sqrt(nU))
		gridx = gridx+1;
		gridy = gridx; 
		peakFE = [peakFE, zeros(1, gridx*gridy-nU)];
		peakRU = [peakRU, zeros(1, gridx*gridy-nU)];
		lagFE = [lagFE, zeros(1, gridx*gridy-nU)];
		lagRU = [lagRU, zeros(1, gridx*gridy-nU)];
		scoreFE = [scoreFE, zeros(1, gridx*gridy-nU)];
		scoreRU = [scoreRU, zeros(1, gridx*gridy-nU)];
	end
	subplot(3,2,1)
	image(reshape(abs(peakFE),gridx,gridy), 'CDataMapping', 'scaled');
	zaxis = [0 max(abs(peakFE))];
	caxis(zaxis);
	set(gca,'Zlim',zaxis,'Ztick',zaxis);
	xlabel('channel');
	title(['Peak FE cross-correlation'])
	colorbar
	subplot(3,2,2)
	image(reshape(abs(peakRU),gridx,gridy), 'CDataMapping', 'scaled');
	zaxis = [0 max(abs(peakRU))];
	caxis(zaxis);
	set(gca,'Zlim',zaxis,'Ztick',zaxis);
	xlabel('channel');
	title(['Peak RU cross-correlation'])
	colorbar
	subplot(3,2,3)
	image(reshape(lagFE,gridx,gridy), 'CDataMapping', 'scaled');
	zaxis = [min(lagFE) max(lagFE)];
	caxis(zaxis);
	set(gca,'Zlim',zaxis,'Ztick',zaxis);
	xlabel('channel');
	title(['Lag FE'])
	colorbar
	subplot(3,2,4)
	image(reshape(lagRU,gridx,gridy), 'CDataMapping', 'scaled');
	zaxis = [min(lagRU) max(lagRU)];
	caxis(zaxis);
	set(gca,'Zlim',zaxis,'Ztick',zaxis);
	xlabel('channel');
	title(['Lag RU'])
	colorbar
	subplot(3,2,5)
	image(reshape(abs(scoreFE),gridx,gridy), 'CDataMapping', 'scaled');
	zaxis = [0 max(abs(scoreFE))];
	caxis(zaxis);
	set(gca,'Zlim',zaxis,'Ztick',zaxis);
	xlabel('channel');
	title(['Score FE cross-correlation' num2str(i)])
	colorbar
	subplot(3,2,6)
	image(reshape(abs(scoreRU),gridx,gridy), 'CDataMapping', 'scaled');
	zaxis = [0 max(abs(scoreRU))];
	caxis(zaxis);
	set(gca,'Zlim',zaxis,'Ztick',zaxis);
	xlabel('channel');
	title(['Score RU cross-correlation' num2str(i)])
	colorbar
	saveplot(gcf, [fn_out '_summary.eps'], 'eps', [6 4]);	
end
