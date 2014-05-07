function [binnedspikes rates torque unitnames] = preprocess_nev(nevfile, binsize, sigma_fr, sigma_trq, offset, verbosity)
	%preprocess_nev	Function to fit the following model to spike and torque data:
	%			lambda_i(t) = \lambda_0 + \sum_j^N k_j^1 . x_i^1(t+\tau+jh) + \sum_j^N k_j^2 . x_i^2(t+\tau+jh)
	%
	% 		Input:
	%			nevfile = file to process
	%			fn_out = base name for output plots
	%			binsize = (optional, default = 0.05) size of timebins over which to compute regression
	%			sigma_fr = (optional, default = 0.25) width of gaussian filter to apply to spikes for firing rate. If 0 then no filter applied
	%			sigma_trq = (optional, default = 0.25) width of gaussian filter to apply to torque. If 0 then no filter applied
	%				Note: Note: both sigmas are in units of seconds, and then are scaled according to binsize
	%			offset = (optional, default = 0) number of seconds to add to spike data before comparing with torque
	%			versbose = (optional, default = 0) verbosity level. If above 0 then plot/print extra info
	%		
	%		Output:
	%			binnedspikes = [nB x nU] array with spikes from all channels binned according to binsize. nB = no. bins, 
	%			rates = filtered (smoothed) binnedspikes data according to sigma_fr filter width
	%			torque = 
	%			unitnames = 
	%
	%		Test code:
	%			nevfile = './testdata/20130117SpankyUtah001.nev';
	%			binsize = 0.05;
	%			kernellength = 6;
	%			sigma_fr = 0.25;
	%			sigma_trq = 0.25;
	%			offset = 0;
	%			verbosity = 1;
	%			fn_out = './worksheets/tuning/01172013/20130117SpankyUtah001';
	%			regression_nev(nevfile, fn_out, kernellength, binsize, sigma_fr, sigma_trq, offset, verbosity);
	
	%Optional arguments
	if (nargin < 3)	kernellength = 6; end
	if (nargin < 4) binsize = 0.05; end
	if (nargin < 5) sigma_fr = 0.25; end
	if (nargin < 6) sigma_trq = 0.25; end
	if (nargin < 7) offset = 0; end
	if (nargin < 8) verbosity = 0; end
	%Set this to above 0 if want debug info/plots
	%verbosity = 1;
	%Total number of possible units recorded from
	nE = 128;
	nunits = 5; 
	nU = nE*nunits;
	%Threshold firing reate below which we ignore that unit
	threshold = 5;
	samplerate = 1/binsize;
	%Make sure we can perform the sample rate conversion easily
	assert(rem(samplerate,1) == 0, 'Select a binsize corresponding to an integer sample rate.');
	ns3file = [nevfile(1:end-3) 'ns3'];

	%Filters
    if sigma_fr > 0
    	sigma_fr = sigma_fr*samplerate;
   		sz = sigma_fr*3*2;
	    x = linspace(-sz/2, sz/2, sz);
	    gaussFilter_fr = exp(-x.^2/(2*sigma_fr^2));
    	gaussFilter_fr = gaussFilter_fr/sum(gaussFilter_fr);
    end
    if sigma_trq > 0
    	sigma_trq = sigma_trq*samplerate;
   		sz = sigma_trq*3*2;
	    x = linspace(-sz/2, sz/2, sz);
	    gaussFilter_trq = exp(-x.^2/(2*sigma_trq^2));
    	gaussFilter_trq = gaussFilter_trq/sum(gaussFilter_trq);
    end
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
    for idx=1:nU
    	if sigma_fr > 0
    		gf = conv(binnedspikes(:,idx), gaussFilter_fr, 'same');
	        rates(:,idx)=gf*samplerate;
    	else
    		rates(:,idx) = binnedspikes(:,idx)*samplerate;
    	end
    end

	%%%%%%%%%%%%%%%%%%%%%
	%Process torque data%
	%%%%%%%%%%%%%%%%%%%%%
	sigma_trq = 0.5*samplerate;
	clear torque;
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
        torque(:,j) = conv(nsxtorque(j,:),gaussFilter_trq,'same');
    end
    %Resample at rate of binsize
    torque=resample(torque,samplerate,nsxsamplerate);
    plot(torque(1:100,1), torque(1:100,2));
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

    if (verbosity == 1)
    	%Plot a bunch of preprocessing diagnostics
    	figure
    	subplot(3,2,1)
    	%Plot the kernel used
    	plot((1:sz)/nsxsamplerate,gaussFilter_trq)
    	title('Gaussian filter used torque')
    	%And make a plot of smoothed data compared with binned spikes
    	subplot(3,2,2);
    	t = 50; unit = 18;
    	times = (1:(t*samplerate))*binsize;
		plot(times, rates(1:(t*samplerate), unit)*binsize, times, binnedspikes(1:(t*samplerate),unit))
		title('Smoothed rate vs binned spikes');
		subplot(3,2,3)
    	t = 10;
    	times = (1:(t*samplerate))*binsize;
    	plot(times, torque(1:(t*samplerate),1),times, torque(1:(t*samplerate),2))
		title('Smoothed torque');		
		subplot(3,2,4)
		%Compute auto- and cross-correlation in torque and example firing rate
		maxlag = 90;
		autotorqueFE = xcov(torque(:,1),samplerate*maxlag);%, 'coeff');
   		autotorqueRU = xcov(torque(:,2),samplerate*maxlag);%, 'coeff');
		covFE = xcov(rates(:,unit), torque(:,1),samplerate*maxlag,'unbiased');
    	% normalize against spikes auto-covariance
    	autorate = xcov(rates(:,unit),samplerate*maxlag);%, 'coeff');
    	covFE = covFE / sqrt(xcov(rates(:,unit),0));
    	covFE = covFE / sqrt(xcov(torque(:,1),0));
   		tt = -maxlag:binsize:maxlag;
    	plot(tt, covFE);
		title(['cross-corr FE, unit ' num2str(unitnames{unit})]);		
		subplot(3,2,5)
		plot(tt, autotorqueFE)
		title('auto-corr torque FE');
		subplot(3,2,6)
		plot(tt, autorate);
		title(['auto-corr rate, unit ' num2str(unitnames{unit})])
		saveplot(gcf, [fn_out '_preprocess.eps'], 'eps', [3 6]);

	end
	display('Printed preprocessing diagnostics.');

    
