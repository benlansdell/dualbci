function processed = preprocess(nevfile, binsize, threshold, offset, fn_out)
	%Preprocess both torque data and firing rate data from an .nev file and a corresponding .ns3 file.
	%Will do the following:
	%- resample spikes and torque data into units of binsize (seconds)
	%- apply a temporal offset between the two
	%- apply a threshold on average firing rate, below which, unit is not returned
	%
	%Usage:
	%		processed = preprocess(nevfile, binsize, threshold, offset, fn_out)
	%	
	%Input:
	%		nevfile = .nev file to process. For loading torque data, assumes that an .nsx file of the same name and location exists.
	%		binsize = (optional, default = 0.05) size of timebins over which to compute regression
	%		threshold = (optional, default = 5) threshold firing rate below which unit is ignored
	%		offset = (optional, default = 0) number of seconds to add to spike data before comparing with torque
	%		fn_out = (optional) file to output diagnostic plots to, if provided
	%	
	%Output:
	%		processed is a structure containing the following fields:
	%			binnedspikes = [nB x nU] array with spikes from all channels binned according to binsize. nB = no. bins, nU = no. units.
	%			rates = [nB x nU] array of binnedspikes multiplied by samplerate to give a per second firing rate approximation
	%			torque = [nB x 2] array of torque inputs 
	%			dtorque = [nB x 2] array of diff of torque inputs (velocity)
	%			ddtorque = [nB x 2] array of diff of diff of torque inputs (acceleration)
	%			unitnames = String of the format Electrode.SortCode used to distinguish individual units from multi-unit electrode activity
	%			tspks = cell array containing spike times for each active channel
	%			binsize = binsize used
	%
	%Test code:
	%		nevfile = './testdata/20130117SpankyUtah001.nev';
	%		binsize = 0.002;
	%		offset = 0.0;
	%		threshold = 5;
	%		fn_out = './worksheets/diagnostics/test_pre.eps';
	%		processed = preprocess(nevfile, binsize, threshold, offset, fn_out);
	
	%Optional arguments
	if (nargin < 2) binsize = 0.05; end
	if (nargin < 3) threshold = 5; end		
	if (nargin < 4) offset = 0; end
	if (nargin < 5) fn_out = 0; end
	%Total number of possible units recorded from
	nE = 128;
	%Total number of sort codes per channel used
	nunits = 5; 
	%Total number of possible units, later truncated
	nU = nE*nunits;
	samplerate = 1/binsize;
	%Make sure we can perform the sample rate conversion easily
	assert(rem(samplerate,1) == 0, 'Select a binsize corresponding to an integer sample rate.');
	%Torque data
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
	isvalid = zeros(1,nU);
	%Prepare spike time structure to bin spikes
	for idx=1:nU
		spikemuas(idx).times = [0];    
	end
	%Bin spikes according to unit identifier
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
		if length(spikemuas(idx).times)>1
			if (spikemuas(idx).times(2)<20) & (spikemuas(idx).times(end)>(dur-20))
				isvalid(idx)=1;
			end
		end
		display(['Electrode.Unit: ' unitnames{idx} ' Spike count: ' num2str(length(spikemuas(idx).times)-1) ' Mean firing rate (Hz): ' num2str(averate(idx))]);
	end
	%Set a threshold firing rate, below which we ignore that unit
	abovethresh = (averate > threshold) & isvalid;
	%Update nU
	nU = sum(abovethresh);
	display(['preprocess: Found ' num2str(nU) ' units above ' num2str(threshold) 'Hz']);
	unitnames = unitnames(abovethresh);
	spikemuas = spikemuas(abovethresh);
	averate = averate(abovethresh);
	%Bin spikes (chronux function)
	binnedspikes = binspikes(spikemuas, samplerate);
	%Compute firing rate simply as spike count divided by binsize
	for idx=1:nU
		rates(:,idx) = binnedspikes(:,idx)*samplerate;
	end
	%%%%%%%%%%%%%%%%%%%%%
	%Process torque data%
	%%%%%%%%%%%%%%%%%%%%%
	clear torque;
	NS3 = openNSx(ns3file, 'read', 'c:138:139');
	nsxtorque = double(NS3.Data);
	nsxsamplerate = double(NS3.MetaTags.SamplingFreq);
	%Switch sign of FE axis for coordinate consistency
	nsxtorque(2,:)=-nsxtorque(2,:);
	nsxpts = ((1:size(nsxtorque,2))-1)/nsxsamplerate;
	pts = ((1:size(binnedspikes,1))-1)/samplerate;
	%Translate, scale and resample torque data
	for j=1:2
		%Scale from uint16 value to proportion
		nsxtorque(j,:) = nsxtorque(j,:)/(2^15);
		%Subtract mean
		nsxtorque(j,:) = nsxtorque(j,:)-mean(nsxtorque(j,:));
		torque(:,j) = resample(nsxtorque(j,:)', samplerate, nsxsamplerate);
	end

	%Check they're the same length, and trim
	nsamp = min(size(torque,1), size(rates,1));
	torque=torque(1:nsamp,:);
	rates = rates(1:nsamp,:);
	binnedspikes = binnedspikes(1:nsamp,:);
	%Apply offset to data
	delaysamples = round(offset*samplerate);
	if (delaysamples > 0)
		binnedspikes = binnedspikes(1+delaysamples:end,:);
		rates = rates(1+delaysamples:end,:);
		torque = torque(1:end-delaysamples,:);
	elseif (delaysamples < 0)
		binnedspikes = binnedspikes(1:end+delaysamples,:);
		rates = rates(1:end+delaysamples,:);
		torque = torque(1-delaysamples:end,:);
	end
	%Compute dtorque and ddtorque
	dtorque = [diff(torque); 0 0];
	ddtorque = [diff(diff(torque)); 0 0; 0 0];

	%Plot a bunch of preprocessing diagnostics, if desired
	if isstr(fn_out)
		figure
		subplot(2,2,1)
		t = 10;
		times = (1:(t*samplerate))*binsize;
		plot(times, torque(1:(t*samplerate),1),times, torque(1:(t*samplerate),2))
		title('Smoothed torque');		
		subplot(2,2,2)
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
		subplot(2,2,3)
		plot(tt, autotorqueFE)
		title('auto-corr torque FE');
		subplot(2,2,4)
		plot(tt, autorate);
		title(['auto-corr rate, unit ' num2str(unitnames{unit})])
		saveplot(gcf, fn_out, 'eps', [6 6]);
	end

	%Return data
	processed.binnedspikes = binnedspikes;
	processed.rates = rates;
	processed.torque = torque;
	processed.dtorque = dtorque; 
	processed.ddtorque = ddtorque;
	processed.unitnames = unitnames;
	processed.tspks = spikemuas;
	processed.binsize = binsize;