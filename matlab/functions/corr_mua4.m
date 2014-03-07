function corr4 = corr_mua4(nevfiles, binsize, fn_out)
	%corr_mua4	Function to compute correlation between spike trains (firing rates) for a set of 
	%		Blackrock recording (nev and ns3) files for each electrode (multi unit) and the 4 cardinal directions. 
	%		Function takes list of nev files, concatenates their spike trains and torque data. Torque data
	%		assumed to reside in file of same name, but with extension .ns3.
	%		(B Lansdell)
	%
	% 		Input:
	%			nevfiles = cell array of nev files to process
	%			binsize = (optional, default = 0.1) size of time bin (in seconds) over which to 
	%				compute firing rate.
	%			fn_out = (optional, default = '') if provided will plot kernel density estimates of distribution
	%		
	%		Output:
	%			corr4 = [nE x 4] array listing correlation between each of nE electrodes and torque motion
	%				in each of 4 cardinal directions
	%
	%		Test code:
	%			nevfiles = {'./blackrock/20130821SpankyUtah017.nev'};
	%			binsize = 0.1;
	%			fn_out = './worksheets/diagnostics/plots/corr_mua4_20130821SpankyUtah017';
	%			corr4 = corr_mua4(nevfiles, binsize, fn_out);
	if (nargin < 3)
		binsize = 0.1;
		fn_out = '';
	elseif (nargin < 2)
		fn_out = '';
	end

	%Number of electrodes recorded from Utah arrays
	nE = 128;

	totaltime = 0;
	corr4 = zeros(nE, 4);
	samplerate = 1/binsize;
	spikes = []; torques = [];

	%Make sure we can perform the sample rate conversion easily
	assert(rem(samplerate,1) == 0, 'Select a binsize corresponding to an integer sample rate.');

	%Load data from each file
	for idx = 1:length(nevfiles)
		nevfile = nevfiles{idx};
		ns3file = [nevfiles{idx}(1:end-3) 'ns3'];
		NEV = openNEV(nevfile);
		%Find the duration and sample rate of the nev file recording
		nevsamplerate = NEV.MetaTags.TimeRes;
		duration = NEV.MetaTags.DataDuration/nevsamplerate;

		%Place the spikes into bins
	        spiketimes = single(NEV.Data.Spikes.TimeStamp)/nevsamplerate;
        	nT = ceil(duration/binsize);
        	nevspikes = zeros(nE, nT);
        	for i=1:length(spiketimes)
                	T = ceil(spiketimes(i)/binsize);
                        E = NEV.Data.Spikes.Electrode(i);
                        nevspikes(E,T) = nevspikes(E,T) + 1;
        	end
		
		%Load torque data from NS3 file
		NS3 = openNSx(ns3file, 'read', 'c:138:139');
		nsxtorque = double(NS3.Data);
		nsxsamplerate = double(NS3.MetaTags.SamplingFreq);
		%Resample at rate of binsize
		torque = resample(nsxtorque', samplerate, nsxsamplerate)';		

		%Concatenate to previously loaded files
		totaltime = totaltime + duration;
		torques = [torques, torque];
		spikes = [spikes, nevspikes];
	end

	%For each electrode
	for idx = 1:nE
		%Compute correlation for each electrode
		corr4(idx, 1) = corr(subplus(torques(1,:)'), spikes(idx,:)');
		corr4(idx, 2) = corr(-subplus(-torques(1,:)'), spikes(idx,:)');
		corr4(idx, 3) = corr(subplus(torques(2,:)'), spikes(idx,:)');
		corr4(idx, 4) = corr(-subplus(-torques(2,:)'), spikes(idx,:)');
	
		%If filename provided, plot density estimates of distributions
		if exist('fn_out', 'var')
			plot(subplus(torques(1,:)), spikes(idx,:), '.');
			title(['correlation = ' num2str(corr4(idx,1))]);
			saveplot(gcf, [fn_out '_channel_' num2str(idx) '_torque1+.eps']);
			plot(-subplus(torques(1,:)), spikes(idx,:), '.');
			title(['correlation = ' num2str(corr4(idx,2))]);
			saveplot(gcf, [fn_out '_channel_' num2str(idx) '_torque1-.eps']);
			plot(subplus(torques(2,:)), spikes(idx,:), '.');
			title(['correlation = ' num2str(corr4(idx,3))]);
			saveplot(gcf, [fn_out '_channel_' num2str(idx) '_torque2+.eps']);
			plot(-subplus(torques(2,:)), spikes(idx,:), '.');
			title(['correlation = ' num2str(corr4(idx,4))]);
			saveplot(gcf, [fn_out '_channel_' num2str(idx) '_torque2-.eps']);
		end
	end
end
