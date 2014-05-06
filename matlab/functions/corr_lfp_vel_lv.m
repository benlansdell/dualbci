function corrbands = corr_lfp_vel(nevfiles, binsize, fn_out, movingwin, offset)
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
	%			offset = (optional, default = 0) number of seconds to add to spike data before computing corr with torque
	%		
	%		Output:
	%			corrbands = [nE x 4 x 4] array listing correlation between each of nE electrodes and torque motion
	%				in each axis for each of 4 frequency bands
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
	%Size of gaussian filter to apply
	sz = 30;
	totaltime = 0;
	corrbands = zeros(nE,4,4);
	p4 = zeros(nE,4,4);
	samplerate = 1/binsize;
	lfp_alpha = [];
	lfp_beta = [];
	lfp_gamma = [];
	lfp_beta = [];
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
		nsxelectrodes = NS3Data(1:128,:);
		nsxtorque = NS3Data(138:139,:);
		nsxsamplerate = double(NS3.MetaTags.SamplingFreq);

		%Params for spectrogram
		params.Fs = nsxsamplerate;
    	params.fpass = [0 1000];
    	params.tapers = [5 9];
    	params.trialave = 0;
    	params.err = 0;
		%Compute LFP data
		[S,t,f] =  mtspecgramc(nsxelectrodes', movingwin, params);	
		size(S)
		
		%Normalize torque data
		for j=1:2
			nsxtorque(j,:) = (nsxtorque(j,:)-mean(nsxtorque(j,:)))/std(nsxtorque(j,:));
		end

		%Resample at rate of binsize
		torque = resample(nsxtorque', samplerate, nsxsamplerate);		
		%t = size(torque);
		%n = size(nevrates);
		%if t(1) < n(1)
		%	torque = [torque; 0 0];
		%end

		%Concatenate to previously loaded files
		totaltime = totaltime + duration;
		vel = [torque(2:end,:) - torque(1:(end-1),:); 0,0];
		vels = [vels; vel];
		lfps_alpha = [lfps_alpha; lfp_alpha];
		lfps_beta = [lfps_beta; lfp_beta];
		lfps_gamma = [lfps_gamma; lfp_gamma];
		lfps_delta = [lfps_delta; lfp_delta];
	end

	%For each electrode compute correlation and make scatter plot
	for idx = 1:nE
		%Compute correlation for each electrode
		[corr2(idx, 1), p2(idx,1)] = corr(vels(:,1), rates(:,idx));
		[corr2(idx, 2), p2(idx,2)] = corr(vels(:,2), rates(:,idx));
	
		%If filename provided, plot density estimates of distributions
		if (length(fn_out) > 0)
			plot(vels(:,1), rates(:,idx), '.');
			title(['correlation = ' num2str(corr2(idx,1)) ' p-value ' num2str(p2(idx,1))]);
			saveplot(gcf, [fn_out '_channel_' num2str(idx) '_torque1.eps']);
			plot(vels(:,2), rates(:,idx), '.');
			title(['correlation = ' num2str(corr2(idx,2)) 'p-value ' num2str(p2(idx,2))]);
			saveplot(gcf, [fn_out '_channel_' num2str(idx) '_torque2.eps']);
		end
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
