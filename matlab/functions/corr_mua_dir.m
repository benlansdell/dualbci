function corr = corr_mua_dir(nevfiles, binsize, fn_out, sigma, offset)
	%corr_mua_dir	Function to compute correlation between spike trains (firing rates) for a set of 
	%		Blackrock recording (nev and ns3) files for each electrode (multi unit) and the cursor direction.
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
	%			corr = [nE x 2] array listing for each electrode the direction (-pi < theta < pi) which has the maximal
	%				correlation. First column lists theta, second column lists correlation measure
	%
	%		Test code:
	%			nevfiles = {'./testdata/20130117SpankyUtah001.nev'};
	%			bs = 0.1;
	%			fn = './worksheets/diagnostics/plots/test_corr_mua_dir_20130117SpankyUtah001';
	%			corr = corr_mua_dir(nevfiles, bs, fn);
	
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
	corr2 = zeros(nE, 2);
	p2 = zeros(nE,2);
	samplerate = 1/binsize;
	rates = [];
	dirs = [];
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

		%Load torque data from NS3 file
		NS3 = openNSx(ns3file, 'read', 'c:138:139');
		nsxtorque = double(NS3.Data);
		nsxsamplerate = double(NS3.MetaTags.SamplingFreq);
		%Normalize torque
		size(nsxtorque)
		for j=1:2
			nsxtorque(j,:) = (nsxtorque(j,:)-mean(nsxtorque(j,:)))/std(nsxtorque(j,:));
		end
		size(nsxtorque)
		%Resample at rate of binsize
		torque = resample(nsxtorque', samplerate, nsxsamplerate);		
		t = size(torque)
		n = size(nevrates)
		if t(1) < n(1)
			torque = [torque; 0 0];
		end
	
		%Compute torque direction
		dir = atan(nsxtorque(2,:)./nsxtorque(1,:));
		%Add pi to angles in second quadrant, subtract pi from angles in 3rd quadrant
		dir = dir + pi*(nsxtorque(1,:)<0).*((nsxtorque(2,:)>0)-(nsxtorque(2,:))<0);

		%Concatenate to previously loaded files
		totaltime = totaltime + duration;
		dirs = [dirs; dir];
		rates = [rates; nevrates];
	end

	%For each electrode compute correlation and make scatter plot
	for idx = 1:nE
		%Fit a cosine function to plots
		%Compute correlation for each electrode
		%[corr2(idx, 1), p2(idx,1)] = corr(torques(:,1), rates(:,idx));
	
		%If filename provided, plot density estimates of distributions
		if (length(fn_out) > 0)
			plot(dirs(:,1), rates(:,idx), '.');
			title(['correlation = ' num2str(corr2(idx,1)) ' p-value ' num2str(p2(idx,1))]);
			saveplot(gcf, [fn_out '_channel_' num2str(idx) '_torque1.eps']);
		end
	end

	%Heat map of correlations and p values
	%if length(fn_out)
	%	for idx = 1:2
	%		image(reshape(corr2(:,idx), 8,16), 'CDataMapping', 'scaled');
	%		zaxis = [-1 1];
   	%		caxis(zaxis);
   	%		%set(gca,'Zlim',zaxis,'Ztick',zaxis, 'NextPlot', 'replacechildren');
   	%		set(gca,'Zlim',zaxis,'Ztick',zaxis);
   	%		xlabel('channel');
 	 %  		title('correlation with torque axis per channel')
    	%		colorbar
	%		saveplot(gcf, [fn_out '_axis_' num2str(idx) '_heatmap.eps']);
%
 %                       image(reshape(p2(:,idx), 8,16), 'CDataMapping', 'scaled');
  %                      zaxis = [0 1];
   %                     caxis(zaxis);
    %                    %set(gca,'Zlim',zaxis,'Ztick',zaxis, 'NextPlot', 'replacechildren');
     %                   set(gca,'Zlim',zaxis,'Ztick',zaxis);
      %                  xlabel('channel');
       %                 title('p-value of non-zero correlation with torque axis per channel')
        %                colorbar
         %               saveplot(gcf, [fn_out '_axis_' num2str(idx) '_heatmap_p-val.eps']);

	%	end
	%end
end
