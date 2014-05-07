function corrbands = corr_lfp_vel(nevfiles, binsize, fn_out, movingwin, bands, offset)
	%corr_lfp_vel	Function to compute correlation between LFP bands for a set of 
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
	%			corrbands = corr_lfp_vel(nevfiles, bs, fn);
	
	if (nargin < 2)
		movingwin = [1 0.1];
		fn_out = '';
		binsize = 0.1;
		offset = 0;
	elseif (nargin < 3)
		fn_out = '';
		movingwin = [1 0.1];
		offset = 0;
	elseif (nargin < 4)
		movingwin = [1 0.1];
		offset = 0;
	elseif (nargin < 5)
		offset = 0;
	end

	%Number of electrodes recorded from Utah arrays
	nE = 128;
	labviewsamplerate = 60;
	totaltime = 0;
	corrbands4 = zeros(nE,4,5);
	p4 = zeros(nE,4,5);
	samplerate = 1/binsize;
	
	params.fpass = [0 1000];
	params.tapers = [5 9];
	params.trialave = 0;
	params.err = 0;

	%Set up frequency bands of interest:
	falpha = [0 5];
	fbeta = [10 40];
	fgamma1 = [45 65];
	fgamma2 = [70 200];
	feps = [200 400];

	lfps_alpha = [];
	lfps_beta = [];
	lfps_gamma1 = [];
	lfps_gamma2 = [];
	lfps_eps = [];
	vels = [];
	%Make sure we can perform the sample rate conversion easily
	assert(rem(samplerate,1) == 0, 'Select a binsize corresponding to an integer sample rate.');

	%Load data from each file
	for idx = 1:length(nevfiles)

		nevfile = nevfiles{idx};
		ns3file = [nevfiles{idx}(1:end-3) 'ns3'];
		%Load data from NS3 file
		NS3 = openNSx(ns3file, 'read');
		NS3Data = double(NS3.Data);
		nsxelectrodes = NS3Data(1:2,:);
		%nsxelectrodes = NS3Data(1:128,:);
		nsxtorque = NS3Data(138:139,:);
		nsxsamplerate = double(NS3.MetaTags.SamplingFreq);

		%Params for spectrogram
		params.Fs = nsxsamplerate;
		%Compute LFP data. Memory intensive... needs to be run on americano.amath or similar
		[S,t,f] =  mtspecgramc(nsxelectrodes', movingwin, params);	
		size(S);

		%Split into frequency bands
		idxalpha = (f>falpha(1)) & (f<falpha(2));
		idxbeta = (f>fbeta(1)) & (f<fbeta(2));
		idxgamma1 = (f>fgamma1(1)) & (f<fgamma1(2));
		idxgamma2 = (f>fgamma2(1)) & (f<fgamma2(2));
		idxeps = (f>feps(1)) & (f<feps(2));

		Salpha = S(:,idxalpha,:);
		Sbeta = S(:,idxbeta,:);
		Sgamma1 = S(:,idxgamma1,:);
		Sgamma2 = S(:,idxgamma2,:);
		Seps = S(:,idxeps,:);

		%Average spectral content over those bands to get one feature per timestep per electrode
		lfp_alpha = squeeze(mean(Salpha,2));
		lfp_beta = squeeze(mean(Sbeta,2));
		lfp_gamma1 = squeeze(mean(Sgamma1,2));
		lfp_gamma2 = squeeze(mean(Sgamma2,2));
		lfp_eps = squeeze(mean(Seps,2));

		%Normalize this by its mean and std
		

		%Normalize torque data
		for j=1:2
			nsxtorque(j,:) = (nsxtorque(j,:)-mean(nsxtorque(j,:)))/std(nsxtorque(j,:));
		end

		%Resample at rate of binsize
		torque = resample(nsxtorque', samplerate, nsxsamplerate);		

		%Concatenate to previously loaded files
		totaltime = totaltime + duration;
		vel = [torque(2:end,:) - torque(1:(end-1),:); 0,0];
		vels = [vels; vel];
		lfps_alpha = [lfps_alpha; lfp_alpha];
		lfps_beta = [lfps_beta; lfp_beta];
		lfps_gamma1 = [lfps_gamma1; lfp_gamma1];
		lfps_gamma2 = [lfps_gamma2; lfp_gamma2];
		lfps_eps = [lfps_eps; lfp_eps];
	end

	%For each electrode compute correlation and make scatter plot
	for i = 1:nE
		%Compute correlation for each electrode for each band
		for j = 1:nbands
			[corr2(i,1,j), p2(i,1,j)] = corr(vels(:,1), lfps(:,i,j));
			[corr2(i,2,j), p2(i,2,j)] = corr(vels(:,2), lfps(:,i,j));
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
		for idx = 1:2
			image(reshape(corr2(:,idx), 8,16), 'CDataMapping', 'scaled');
			zaxis = [-1 1];
            %xaxis = [0 112];
   			caxis(zaxis);
   			%set(gca,'Zlim',zaxis,'Ztick',zaxis, 'NextPlot', 'replacechildren');
   			set(gca,'Zlim',zaxis,'Ztick',zaxis);
   			xlabel('channel');
 	   		title('correlation with torque axis per channel')
    		colorbar
			saveplot(gcf, [fn_out '_axis_' num2str(idx) '_heatmap.eps']);
            image(reshape(p2(:,idx), 8,16), 'CDataMapping', 'scaled');
            zaxis = [0 1];
            caxis(zaxis);
            %set(gca,'Zlim',zaxis,'Ztick',zaxis, 'NextPlot', 'replacechildren');
            set(gca,'Zlim',zaxis,'Ztick',zaxis);
            xlabel('channel');
            title('p-value of non-zero correlation with torque axis per channel')
            colorbar
            saveplot(gcf, [fn_out '_axis_' num2str(idx) '_heatmap_p-val.eps']);

		end
	end
end
