function corrbands = corr_lfp_spd(nevfiles, binsize, fn_out, movingwin, bands, offset)
	%corr_lfp_spd	Function to compute correlation between LFP bands for a set of 
	%		Blackrock recording (ns3) files for each electrode and the torque velocity.
	%		Function takes list of nev files, concatenates their LFP and torque data. Torque data
	%		assumed to reside in file of same name, but with extension .ns3.
	%		(B Lansdell)
	%
	% 		Input:
	%			nevfiles = cell array of nev files to process
	%			binsize = (optional, default = 0.1) size of window over which to compute correlation between torque and FRs
	%			fn_out = (optional, default = '') if provided will plot kernel density estimates of distribution
	%			movingwin = (optional, default = [1 0.1]) window size and window-step size used to compute spectrogram
	%			bands = (optional, default = {[0, 5],[10, 40],[45, 65],[70, 200],[200, 400]}) the frequency bands used
	%			offset = (optional, default = 0) number of seconds to add to spike data before computing corr with torque
	%		
	%		Output:
	%			corrbands = [nE x 4 x 5] array listing correlation between each of nE electrodes and torque motion
	%				in each axis for each of 5 frequency bands of interest
	%
	%		Test code:
	%			nevfiles = {'./testdata/20130117SpankyUtah001.nev'};
	%			bs = 0.1;
	%			fn = './worksheets/diagnostics/plots/test_corr_lfp_vel_20130117SpankyUtah001';
	%			corrbands = corr_lfp_spd(nevfiles, bs, fn);
	
	if (nargin < 2)
                movingwin = [1 0.1];
		bands = {[0, 5],[10, 40],[45, 65],[70, 200],[200, 400]};
                fn_out = '';
                binsize = 0.1;
                offset = 0;		
	elseif (nargin < 3)
		movingwin = [1 0.1];
 		bands = {[0, 5],[10, 40],[45, 65],[70, 200],[200, 400]};
		fn_out = '';
		offset = 0;
	elseif (nargin < 4)
		bands = {[0, 5],[10, 40],[45, 65],[70, 200],[200, 400]};
		movingwin = [1 0.1];
		offset = 0;
	elseif (nargin < 5)
		bands = {[0, 5],[10, 40],[45, 65],[70, 200],[200, 400]};
		offset = 0;
	elseif (nargin < 6)
		offset = 0;
	end

	nbands = length(bands);

	%Number of electrodes recorded from Utah arrays
	nE = 128;
	labviewsamplerate = 60;
	totaltime = 0;
	corrbands = zeros(nE,nbands);
	p = zeros(nE,nbands);
	samplerate = 1/binsize;
	
	params.fpass = [0 1000];
	params.tapers = [5 9];
	params.trialave = 0;
	params.err = 0;

	spds = [];
	lfps = [];
	%Make sure we can perform the sample rate conversion easily
	assert(rem(samplerate,1) == 0, 'Select a binsize corresponding to an integer sample rate.');

	%Load data from each file
	for idx = 1:length(nevfiles)

		nevfile = nevfiles{idx};
		ns3file = [nevfiles{idx}(1:end-3) 'ns3'];
		%Load data from NS3 file
		NS3 = openNSx(ns3file, 'read');
		NS3Data = double(NS3.Data);
		%Test run
		%nsxelectrodes = NS3Data(1:2,:);
		nsxelectrodes = NS3Data(1:128,:);
		nsxtorque = NS3Data(138:139,:);
		nsxsamplerate = double(NS3.MetaTags.SamplingFreq);

		%Params for spectrogram
		params.Fs = nsxsamplerate;
		%Compute LFP data. Memory intensive... needs to be run on americano.amath or similar
		[S,t,f] =  mtspecgramc(nsxelectrodes', movingwin, params);	
		size(S);

		%Split into frequency bands
		for j = 1:nbands
			idxb = (f>bands{idx}(1)) & (f<bands{idx}(2));
			Sb = S(:,idxb,:);
			%Average spectral content over those bands to get one feature per timestep per electrode
			lfpb = squeeze(mean(Sb,2));
			lfp(j,:,:) = lfpb;
			%Normalize this by its mean and std
		end

		%Normalize torque data
		for j=1:2
			nsxtorque(j,:) = (nsxtorque(j,:)-mean(nsxtorque(j,:)))/std(nsxtorque(j,:));
		end

		%Resample at rate of binsize
		torque = resample(nsxtorque', samplerate, nsxsamplerate);		

		%Concatenate to previously loaded files
		%totaltime = totaltime + duration;
		vel = [torque(2:end,:) - torque(1:(end-1),:); 0,0];
		spd = sqrt(vel(:,1).^2 + vel(:,2).^2);
		%Truncate on either side to match with specgram data
		truncby = movingwin(1)/movingwin(2)/2
		spd = spd((truncby+1):(end-truncby),:);
                spds = [spds; spd];
		for j = 1:nbands
			if idx == 1
				lfps(j,:,:) = lfp(j,:,:);
			else
				lfps(j,:,:) = [lfps(j,:,:); lfp(j,:,:)];
			end
		end
	end

	%For each electrode compute correlation and make scatter plot
	for i = 1:nE
	%for i = 1:2
		%Compute correlation for each electrode for each band
		for j = 1:nbands
			[corrbands(i,j), p(i,j)] = corr(spds(:,1), lfps(j,:,i)');
		end
		%If filename provided, plot density estimates of distributions
		%if (length(fn_out) > 0)
		%	plot(vels(:,1), rates(:,i), '.');
		%	title(['correlation = ' num2str(corr2(i,1)) ' p-value ' num2str(p2(i,1))]);
		%	saveplot(gcf, [fn_out '_channel_' num2str(i) '_torque1.eps']);
		%	plot(vels(:,2), rates(:,i), '.');
		%	title(['correlation = ' num2str(corr2(i,2)) 'p-value ' num2str(p2(i,2))]);
		%	saveplot(gcf, [fn_out '_channel_' num2str(i) '_torque2.eps']);
		%end
	end

	%Heat map of correlations and p values
	if length(fn_out)
		for idx = 1:nbands
			image(reshape(corrbands(:,idx), 8,16), 'CDataMapping', 'scaled');
			zaxis = [-1 1];
        		%xaxis = [0 112];
   			caxis(zaxis);
   			%set(gca,'Zlim',zaxis,'Ztick',zaxis, 'NextPlot', 'replacechildren');
   			set(gca,'Zlim',zaxis,'Ztick',zaxis);
   			xlabel('channel');
 	   		title(['correlation with speed per channel. band ' num2str(bands{idx}(1)) ' to ' num2str(bands{idx}(2)) 'Hz'])
    			colorbar
			saveplot(gcf, [fn_out '_band_' num2str(idx) '_heatmap.eps']);
        		image(reshape(p(:,idx), 8,16), 'CDataMapping', 'scaled');
        		zaxis = [0 1];
        		caxis(zaxis);
        		%set(gca,'Zlim',zaxis,'Ztick',zaxis, 'NextPlot', 'replacechildren');
        		set(gca,'Zlim',zaxis,'Ztick',zaxis);
       			xlabel('channel');
        		title(['p-value of non-zero correlation with speed per channel. band ' num2str(bands{idx}(1)) ' to ' num2str(bands{idx}(2)) 'Hz'])
        		colorbar
        		saveplot(gcf, [fn_out '_band_' num2str(idx) '_heatmap_p-val.eps']);
		end
	end
end
