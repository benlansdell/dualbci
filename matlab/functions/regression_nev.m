function regression_nev(nevfile, fn_out, kernellength, binsize, sigma, offset)
	%regression_nev	Function to fit the following model to spike and torque data:
	%			lambda_i(t) = \lambda_0 + \sum_j^N k_j^1 . x_i^1(t+\tau+jh) + \sum_j^N k_j^2 . x_i^2(t+\tau+jh)
	%		That is, it fits a linear filter to the torque data. In the above formula, time-step size h
	%		is given by parameter binsize; max value of N is given by kernellength; \tau is given by offset.
	%		Program fits all models with kernels of length between 1 and N, meaning that for each single-unit N models
	%		will be fit.
	%
	% 		Input:
	%			nevfile = file to process
	%			fn_out = base name for output plots
	%			kernellength = (optional, default = 6) max 
	%			binsize = (optional, default = 0.05) size of timebins over which to compute regression
	%			sigma = (optional, default = 5) width of gaussian filter to apply to spikes for firing rate
	%			offset = (optional, default = 0) number of seconds to add to spike data before comparing with torque
	%		
	%		Output:
	%			(none) produces plots of the filter for each single-unit channel, along with summary of quality of fits for 
	%			each channel
	%
	%		Test code:
	%			nevfile = './testdata/20130117SpankyUtah001.nev';
	%			bs = 0.05;
	%			kl = 6;
	%			fn = './worksheets/diagnostics/plots/test_regression_nev_20130117SpankyUtah001';
	%			regression_nev(nevfiles, fn, kl, bs);
	
	%Optional arguments
	if (nargin < 3)	kernellength = 6; end
	if (nargin < 4) binsize = 0.05; end
	if (nargin < 5) sigma = 5; end
	if (nargin < 6) offset = 0; end
	%Set this to above 0 if want debug info/plots
	verbosity = 1;
	%Total number of possible units recorded from
	nE = 128;
	nunits = 5; 
	nU = nE*nunits;
	%Threshold firing reate below which we ignore that unit
	threshold = 1;
	%Size of gaussian filter to apply
	sz = 30;
	samplerate = 1/binsize;
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
        rates(:,idx)=gf;
    end
    if (verbosity)
    	%Plot the kernel used
    	plot((1:sz)/samplerate,gaussFilter)
    	%And make a plot of smoothed data compared with binned spikes
    	figure
    	t = 50; unit = 18;
    	times = (1:(t*samplerate))*binsize;
		plot(times, rates(1:(t*samplerate), unit), times, binnedspikes(1:(t*samplerate),unit))
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
    if (verbosity)
    	plot(torque(1:100,1),torque(1:100,2))
    end

    %%%%%%%%%%%%
    %Fit models%
    %%%%%%%%%%%%
    %Produce shifted torque data to fit model to
    laggedtorque = [];
    for i=1:2
		for j=1:kernellength
	    	laggedtorque = [laggedtorque, torque(j:end-kernellength+j,i)];
    	end
    end
    rates = rates(1:end-kernellength+1,:);
    szRates = size(rates,1);
    fittedrates = zeros(kernellength, nU, szRates);
    %Do fitting for a range of kernel lengths
    models = {};
    lambda_0 = zeros(kernellength, nU);
    kernelRU = zeros(kernellength, nU, kernellength);
    kernelFE = zeros(kernellength, nU, kernellength);
    r2 = zeros(kernellength, nU);
    sumr2 = zeros(kernellength, 1);
    for i=1:kernellength
		display(['Fitting model with kernel length: ' num2str(i)])
		for j=1:nU
			%Estimate parameters as a linear model
			X = [laggedtorque(:,1:i), laggedtorque(:,kernellength+1:kernellength+i)];
			Y = rates(:,j);
			model = LinearModel.fit(X,Y);
			models{i,j} = model;
			r2(i,j) = model.Rsquared.Ordinary;
			coeffs = model.Coefficients.Estimate;
			lambda_0(i,j) = coeffs(1);
			kernelRU(i,j,1:i) = coeffs(2:i+1);
			kernelFE(i,j,1:i) = coeffs(i+2:end);
			%Do the same thing but do the least squares fit ourselves...
			%%Add a constant intercept term
			X = [ones(szRates, 1), X];
			%%Solve normal equations
			beta = (X'*X)\(X'*Y);
	    	%%...these give the same result, but it's easy to get the predicted Y values from here
	    	fittedrates(i,j,:) = X*beta;
	    end
    end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   	%Plot kernels for each neuron%
   	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   	colors = get(gca, 'colororder');
    for i=1:nU
    	figure
    	keys = {};
    	subplot(1,2,1)
		for j=1:kernellength
			hold on
			plot((1:j)*binsize, squeeze(kernelRU(j,i,1:j)), 'Color', colors(j,:));
			keys = {keys{:}, ['R^2=' num2str(r2(j,i)) ', lambda_0=' num2str(lambda_0(j,i))]}
		end
		title(['Unit: ' unitnames(i)])
		xlabel('time (s)')
		ylabel('k [Radial-Ulnar]')
		subplot(1,2,2)
		for j=1:kernellength
			hold on
			plot((1:j)*binsize, squeeze(kernelFE(j,i,1:j)), 'Color', colors(j,:));
		end
		legend(keys)
		title(['Fitted lambda_0: ' num2str(lambda_0(j,i))])
		xlabel('time (s)')
		ylabel('k [Flexion-Extension]')
		saveplot(gcf, [fn_out '_kernels_unit_' unitnames(i) '.eps'])
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Heat map of goodness-of-fit for all units%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for i=1:kernellength
		subplot(kernellength,1,i)
		gridx = floor(sqrt(nU));
		if (gridx < sqrt(nU)) 
			gridy = gridx+1; 
			r2_imputed = [squeeze(r2(i,:)), ones(1,gridx*(gridx+1)-nU)];
		else
			r2_imputed = [squeeze(r2(i,:))];
		end
		image(reshape(r2_imputed,gridx,gridy), 'CDataMapping', 'scaled');
		zaxis = [0 1];
		caxis(zaxis);
		set(gca,'Zlim',zaxis,'Ztick',zaxis);
		xlabel('channel');
  		title(['R^2 values for kernels of length ' num2str(i)])
   		colorbar
	end
	saveplot(gcf, [fn_out '_R2_heatmap.eps'], 'eps', [3 10]);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Plot filtered data%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for i=1:kernellength
		for j=1:nU
			%Plot histogram of filtered data
			%hist(squeeze(fittedrates(i,j,:));
			%Plot scatter plot of filtered data vs firing rate
			%plot(fittedrates(i,j,:), rates(:,j));
			%smoothhist2D([reshape(squeeze(fittedrates(i,j,:)),szRates,1), rates(:,j)], 5, [100 100], 0.05);
		end
	end

end
