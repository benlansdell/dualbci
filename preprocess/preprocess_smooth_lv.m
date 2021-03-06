function processed = preprocess_smooth_lv(nevfile, matfile, binsize, sigma_fr, sigma_trq, threshold, offset, fn_out, verbose, units)
	%Preprocess both torque data and firing rate data from an .nev file and a corresponding .ns3 file.
	%Will do the following:
	%- resample spikes and torque data into units of binsize (seconds)
	%- apply a Guassian filter of width sigma_fr to the binned spikes to compute a firing rate
	%- apply a Gaussian filter of width sigma_trq to the torque channels to smooth torque input
	%- apply a threshold on average firing rate, below which, unit is not returned
	%- apply a temporal offset between the two
	%- plot some diagnostics if desired
	%
	%Usage:
	%		processed = preprocess_smooth_lv(nevfile, matfile, binsize, sigma_fr, sigma_trq, threshold, offset, fn_out, verbose, units)
	%
	%Input:
	%		nevfile = file to process. For loading torque data, assumes that an .nsx file of the same name and location exists.
	%		matfile = matfile containing cursor information for nev file
	%		binsize = (optional, default = 0.05) size of timebins over which to compute regression
	%		sigma_fr = (optional, default = 0.25) width of gaussian filter to apply to spikes for firing rate. If 0 then no filter applied
	%		sigma_trq = (optional, default = 0.25) width of gaussian filter to apply to torque. If 0 then no filter applied
	%			Note: Note: both sigmas are in units of seconds, and then are scaled according to binsize
	%		threshold = (optional, default = 5) threshold firing rate below which unit is ignored
	%		offset = (optional, default = 0) number of seconds to add to spike data before comparing with torque
	%		fn_out = (optional) If provided then plot/print extra info
	%		verbose = (optional) if provided then will output info about firing of units
	%		units = (optional) if a cell array of unit names is provided then will load only these, and ignore the threshold constraint
	%	
	%Output:
	%	processed is a structure containing the following fields:
	%		binnedspikes = [nB x nU] array with spikes from all channels binned according to binsize. nB = no. bins, nU = no. units.
	%		rates = [nB x nU] array of binnedspikes multiplied by samplerate to give a per second firing rate approximation
	%		torque = [nB x 2] array of torque inputs 
	%		dtorque = [nB x 2] array of diff of torque inputs (approximation of velocity)
	%		ddtorque = [nB x 2] array of diff of diff of torque inputs (approximation of acceleration)
	%		cursor = [nB x 2] array of cursor outputs 
	%		dcursor = [nB x 2] array of diff of cursor outputs (approximation of velocity)
	%		ddcursor = [nB x 2] array of diff of diff of cursor outputs (approximation of acceleration)
	%		unitnames = String of the format Electrode.SortCode used to distinguish individual units from multi-unit electrode activity
	%		tspks = cell array containing spike times for each active channel
	%		binsize = binsize used
	%
	%Test code:
	%		nevfile = './testdata/20130117SpankyUtah001.nev';
	%		matfile = './testdata/Spanky_2013-01-17-1325.mat';
	%		binsize = 0.05;
	%		sigma_fr = 0.25;
	%		sigma_trq = 0.25;
	%		offset = 0;
	%		threshold = 5;
	%		fn_out = './worksheets/diagnostics/plots/test_smooth_pre.eps';
	%		verbose = 0;
	%		processed = preprocess_smooth_lv(nevfile, matfile, binsize, sigma_fr, sigma_trq, threshold, offset, fn_out, verbose);

	%Optional arguments
	if (nargin < 3) binsize = 0.05; end
	if (nargin < 4) sigma_fr = 0.25; end
	if (nargin < 5) sigma_trq = 0.25; end
	if (nargin < 6) threshold = 5; end		
	if (nargin < 7) offset = 0; end
	if (nargin < 8) fn_out = 0; end
	if (nargin < 9) verbose = 0; end
	if (nargin < 10) units = {}; end

	nE = 128;
	nunits = 5; 
	%Total number of possible units recorded from
	nU = nE*nunits;
	samplerate = 1/binsize;
	%Make sure we can perform the sample rate conversion easily
	assert(rem(samplerate,1) == 0, 'Select a binsize corresponding to an integer sample rate.');
	ns3file = [nevfile(1:end-3) 'ns3'];
	%Create Filters
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
	isvalid = zeros(1,nU);
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
		if length(spikemuas(idx).times)>1
			if (spikemuas(idx).times(2)<20) & (spikemuas(idx).times(end)>(dur-20))
				isvalid(idx)=1;
			end
		end
		if verbose ~= 0
			display(['Electrode.Unit: ' unitnames{idx} ' Spike count: ' num2str(length(spikemuas(idx).times)-1) ' Mean firing rate (Hz): ' num2str(averate(idx))]);
		end
	end
	if length(units) == 0
		%Set a threshold firing rate, below which we ignore that unit
		abovethresh = (averate > threshold) & isvalid;
		%Update nU
		nU = sum(abovethresh);
		display(['Found ' num2str(nU) ' units above ' num2str(threshold) 'Hz']);
		unitnames = unitnames(abovethresh);
		spikemuas = spikemuas(abovethresh);
		averate = averate(abovethresh);
	else
		indices = cellfun(@(x)ismember(num2str(x), units), unitnames);
		if ~any(indices)
			display(['Cannot find units in nev file. Please provide cell array of strings between 1-128 '...
				'followed by a sort code of 0-3.'])
		end
		unitnames = unitnames(indices);
		spikemuas = spikemuas(indices);
		averate = averate(indices);
		nU = sum(indices);
	end

	%Bin spikes (chronux function)
	binnedspikes = binspikes(spikemuas, samplerate);
	%From this apply gaussian filter to spike train for each electrode
	for idx=1:nU
		if sigma_fr > 0
			gf = conv(binnedspikes(:,idx), gaussFilter_fr, 'same');
			rates(:,idx)=gf*samplerate;
			binnedspikes(:,idx) = rates(:,idx)/samplerate;
		else
			rates(:,idx) = binnedspikes(:,idx)*samplerate;
		end
	end
	tspks = spikemuas;

	%%%%%%%%%%%%%%%%%%%%%
	%Process torque data%
	%%%%%%%%%%%%%%%%%%%%%
	clear torque;

	%Load from labview, also. Needed for brain control datasets
	load(matfile);
	labviewsamplerate = data.sampleRate;
	foundnev = 0;
	[pathstr,nevfilename,ext] = fileparts(nevfile);
	for idx = 1:length(data.nev)
		%Find nev structure
		if length(data.nev(idx).nevfile)==0
			continue
		end
		[pathstr,curr_nevfilename,ext] = fileparts(data.nev(idx).nevfile);
		if strcmp(curr_nevfilename, nevfilename)
			%Find offset and duration in seconds
			nevoffset = single(data.nev(idx).Toffset(1))/60;
			nevdur = data.nev(idx).DurationSec;
			foundnev = 1;
			break
		end 
	end
	assert(foundnev == 1, 'Cannot find nev file in labview mat file, check your input files');
	%extract from stateHist.cursor the cursor position during the nev recording
	withinnev = (data.stateHist.time < nevoffset+nevdur) & (data.stateHist.time > nevoffset);
	lvtorque = data.stateHist.cursor(withinnev,:);
	lvdtorque = [diff(lvtorque); 0, 0]*labviewsamplerate;
	lvddtorque = [diff(lvdtorque); 0, 0]*labviewsamplerate;

	%Resample at desired rate
	lvtorque = resample(lvtorque, samplerate, labviewsamplerate);
	lvdtorque = resample(lvdtorque, samplerate, labviewsamplerate);
	lvddtorque = resample(lvddtorque, samplerate, labviewsamplerate);

	%Smooth with same filter...
	%for j = 1:2
	%	lvtorque(:,j) = conv(lvtorque(:,j),gaussFilter_trq,'same');
	%	lvdtorque(:,j) = conv(lvdtorque(:,j),gaussFilter_trq,'same');
	%	lvddtorque(:,j) = conv(lvddtorque(:,j),gaussFilter_trq,'same');
	%end

	%%%%%%%%%%%%%%%%%%
	%Load torque data%
	%%%%%%%%%%%%%%%%%%

	%Torque data are in last two channels
	chans = findTorqueChannels(nevfile);
	NS3 = openNSx(ns3file, 'read', ['c:' num2str(chans(1)) ':' num2str(chans(2))]);
	nsxtorque = double(NS3.Data);
	nsxsamplerate = double(NS3.MetaTags.SamplingFreq);
	%Switch sign of FE axis for coordinate consistency
	nsxtorque(1,:)=-nsxtorque(1,:);
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
	%Compute dtorque and ddtorque
	dtorque = [diff(torque); 0 0];
	ddtorque = [diff(diff(torque)); 0 0; 0 0];
	if isstr(fn_out)
		plot(torque(1:100,1), torque(1:100,2));
	end

	%Check they're the same length, and trim
	nsamp = min([size(torque,1), size(rates,1), size(lvtorque,1)]);
	torque=torque(1:nsamp,:);
	dtorque=dtorque(1:nsamp,:);
	ddtorque=ddtorque(1:nsamp,:);
	lvtorque=lvtorque(1:nsamp,:);
	lvdtorque=lvdtorque(1:nsamp,:);
	lvddtorque=lvddtorque(1:nsamp,:);
	rates = rates(1:nsamp,:);
	binnedspikes = binnedspikes(1:nsamp,:);
	%Apply offset to data
	delaysamples = round(offset*samplerate);
	if (delaysamples > 0)
		binnedspikes = binnedspikes(1+delaysamples:end,:);
		rates = rates(1+delaysamples:end,:);
		torque = torque(1:end-delaysamples,:);
		dtorque = dtorque(1:end-delaysamples,:);
		ddtorque = ddtorque(1:end-delaysamples,:);
		lvtorque = lvtorque(1:end-delaysamples,:);
		lvdtorque = lvdtorque(1:end-delaysamples,:);
		lvddtorque = lvddtorque(1:end-delaysamples,:);
	elseif (delaysamples < 0)
		binnedspikes = binnedspikes(1:end+delaysamples,:);
		rates = rates(1:end+delaysamples,:);
		torque = torque(1-delaysamples:end,:);
		dtorque = dtorque(1-delaysamples:end,:);
		ddtorque = ddtorque(1-delaysamples:end,:);
		lvtorque = lvtorque(1-delaysamples:end,:);
		lvdtorque = lvdtorque(1-delaysamples:end,:);
		lvddtorque = lvddtorque(1-delaysamples:end,:);
	end
	
	if isstr(fn_out)
		%Plot a bunch of preprocessing diagnostics
		figure
		subplot(3,2,1)
		%Plot the kernel used
		plot((1:sz)/nsxsamplerate,gaussFilter_trq)
		title('Gaussian filter used torque')
		%And make a plot of smoothed data compared with binned spikes
		subplot(3,2,2);
		t = 50; 
		unit = min(18, nU);
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
		saveplot(gcf, fn_out, 'eps', [3 6]);
	end

	%Return data
	processed.binnedspikes = binnedspikes;
	processed.rates = rates;
	processed.torque = torque;
	processed.dtorque = dtorque; 
	processed.ddtorque = ddtorque;
	processed.cursor = lvtorque;
	processed.dcursor = lvdtorque; 
	processed.ddcursor = lvddtorque;
	processed.unitnames = unitnames;
	processed.tspks = tspks;
	processed.binsize = binsize;
	processed.nevfile = nevfile;
	processed.labviewfile = matfile;