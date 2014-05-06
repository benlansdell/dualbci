function corr4 = corr_mua4(nevfiles, binsize, fn_out, sigma, offset)
	%corr_mua4	Function to compute correlation between spike trains (firing rates) for a set of 
	%		Blackrock recording (nev and ns3) files for each electrode (multi unit) and the 4 cardinal directions. 
	%		Function takes list of nev files, concatenates their spike trains and torque data. Torque data
	%		assumed to reside in file of same name, but with extension .ns3.
	%		(B Lansdell)
	%
	% 		Input:
	%			nevfiles = cell array of nev files to process
	%			binsize = (optional, default = 0.1) size of window over which to compute correlation between torque and FRs
	%			fn_out = (optional, default = '') if provided will plot kernel density estimates of distribution
	%			sigma = (optional, default = 5) width of gaussian filter to apply to spikes for firing rate
	%			offset = (optional, default = 0) number of seconds to add to spike data before computing corr with torque
	%		
	%		Output:
	%			corr4 = [nE x 4] array listing correlation between each of nE electrodes and torque motion
	%				in each of 4 cardinal directions
	%
	%		Test code:
	%			nevfiles = {'./testdata/20130117SpankyUtah001.nev'};
	%			bs = 0.1;
	%			fn = './worksheets/diagnostics/plots/test_corr_mua4_20130117SpankyUtah001';
	%			corr4 = corr_mua4(nevfiles, bs, fn);
	
	if (nargin < 2)
		sigma = 5;
		fn_out = '';
		binsize = 0.1;
		offset = 0;
	elseif (nargin < 3)
		fn_out = '';
		sigma = 5;
		offset = 0;
	elseif (nargin < 4)
		sigma = 5;
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
	corr4 = zeros(nE, 4);
	p4 = zeros(nE,4);
	samplerate = 1/binsize;
	rates = []; torques = [];
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

		spiketimes = double(NEV.Data.Spikes.TimeStamp)/nevsamplerate;
        elecs = cell(1,nE);
        spikemuas = struct('times', elecs);
        for idx=1:nE
                spikemuas(idx).times = [0];
        end
        for i=1:length(spiketimes)
               	E = NEV.Data.Spikes.Electrode(i);
               	spikemuas(E).times = [spikemuas(E).times; spiketimes(i)];
        end
        	%Bin spikes
        binnedspikes = binspikes(spikemuas, labviewsamplerate);

	    %From this apply gaussian filter to spike train for each electrode
	    x = linspace(-sz/2, sz/2, sz);
	    gaussFilter = exp(-x.^2/(2*sigma^2));
	    gaussFilter = gaussFilter/sum(gaussFilter);
	    for idx=1:nE
	            gf = conv(binnedspikes(:,idx), gaussFilter, 'same');
                nevrates(:,idx) = resample(gf, samplerate, labviewsamplerate)*samplerate;
        end
		size(nevrates)

		%Load torque data from NS3 file
		NS3 = openNSx(ns3file, 'read', 'c:138:139');
		nsxtorque = double(NS3.Data);
		nsxsamplerate = double(NS3.MetaTags.SamplingFreq);
		%Smooth and subtract mean from data
                for j=1:2
                        %nsxtorque(j,:) = conv(nsxtorque(j,:), gaussFilter, 'same')-mean(nsxtorque(j,:));
                        nsxtorque(j,:) = nsxtorque(j,:)-mean(nsxtorque(j,:));
                end
		%Resample at rate of binsize
		torque = resample(nsxtorque', samplerate, nsxsamplerate);		
                t = size(torque)
                n = size(nevrates)
                if t(1) < n(1)
                        torque = [torque; 0 0];
                end

		size(torque)

		%Concatenate to previously loaded files
		totaltime = totaltime + duration;
		torques = [torques; torque];
		rates = [rates; nevrates];
	end

	%For each electrode compute correlation and make scatter plot
	for idx = 1:nE
		%Compute correlation for each electrode
		%Find those values which are positive, those which are negative
		pindices1 = torques(:,1)>0;
		nindices1 = torques(:,1)<0;
		pindices2 = torques(:,2)>0;
		nindices2 = torques(:,2)<0;
		%[corr4(idx, 1), p4(idx,1)] = corr(subplus(torques(:,1)), rates(:,idx));
		%[corr4(idx, 2), p4(idx,2)] = corr(-subplus(-torques(:,1)), rates(:,idx));
		%[corr4(idx, 3), p4(idx,3)] = corr(subplus(torques(:,2)), rates(:,idx));
		%[corr4(idx, 4), p4(idx,4)] = corr(-subplus(-torques(:,2)), rates(:,idx));
		[corr4(idx, 1), p4(idx,1)] = corr(torques(pindices1,1), rates(pindices1,idx));
		[corr4(idx, 2), p4(idx,2)] = corr(torques(nindices1,1), rates(nindices1,idx));
		[corr4(idx, 3), p4(idx,3)] = corr(torques(pindices2,2), rates(pindices2,idx));
		[corr4(idx, 4), p4(idx,4)] = corr(torques(nindices2,2), rates(nindices2,idx));

	
		%If filename provided, plot density estimates of distributions
		if (length(fn_out) > 0)
			plot(torques(pindices1,1), rates(pindices1,idx), '.');
			title(['correlation = ' num2str(corr4(idx,1))]);
			saveplot(gcf, [fn_out '_channel_' num2str(idx) '_torque1+_corr_' num2str(corr4(idx,1)) '.eps']);
			plot(torques(nindices1,1), rates(nindices1,idx), '.');
			title(['correlation = ' num2str(corr4(idx,2))]);
			saveplot(gcf, [fn_out '_channel_' num2str(idx) '_torque1-_corr_' num2str(corr4(idx,2)) '.eps']);
			plot(torques(pindices2,2), rates(pindices2,idx), '.');
			title(['correlation = ' num2str(corr4(idx,3))]);
			saveplot(gcf, [fn_out '_channel_' num2str(idx) '_torque2+_corr_' num2str(corr4(idx,3)) '.eps']);
			plot(torques(nindices2,2), rates(nindices2,idx), '.');
			title(['correlation = ' num2str(corr4(idx,4))]);
			saveplot(gcf, [fn_out '_channel_' num2str(idx) '_torque2-_corr_' num2str(corr4(idx,4)) '.eps']);
		end
	end

	%Heat map of correlations and p-values
	if length(fn_out)
		for idx = 1:4
			image(reshape(corr4(:,idx), 8,16), 'CDataMapping', 'scaled');
			zaxis = [-1 1];
   			caxis(zaxis);
   			set(gca,'Zlim',zaxis,'Ztick',zaxis);
   			xlabel('channel');
 	   		title('correlation with torque axis per channel')
    		colorbar
			saveplot(gcf, [fn_out '_axis_' num2str(idx) '_heatmap.eps']);
			image(reshape(p4(:,idx), 8,16), 'CDataMapping', 'scaled');
            zaxis = [0 1];
            caxis(zaxis);
            set(gca,'Zlim',zaxis,'Ztick',zaxis);
            xlabel('channel');
            title('p-value of non-zero correlation with torque axis per channel')
            colorbar
            saveplot(gcf, [fn_out '_axis_' num2str(idx) '_heatmap_p-val.eps']);
		end
	end
end
