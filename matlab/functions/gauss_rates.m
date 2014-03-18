function trial_out = gauss_rates(trial_in, sigma, sz)
    %gauss_rates		Compute a smoothed estimate of firing rate
    %
    % Usage:
    %					gauss_rates(trial_in)
    %
    % Input:
    %					trial_in = trial structure from import_trials
    %					sigma = (optional, default = 5) width of Gaussian 
    %					sz = (optional, default = 30) number of points in filter
	%
	% Output:
	%					trial_out = output trial structure with .nevrates field added
    %
    % Examples:
    %					trials = import_trials('Spanky_2013-01-17-1325.mat');
	%					trial = import_spikes(trials(117));
	%					
	if nargin < 1
		throw(MException('Argin:MoreExpected', 'More input arguments expected'));
	elseif nargin < 2
		sigma = 5;
		sz = 30;
	elseif nargin < 3
		sz = 30;
	end

	trial_out = trial_in;
	nE = length(trial_in.spikemuas);
	%Bin spikes into small bins
	binnedspikes = binspikes(trial_in.spikemuas, trial_in.samplerate, [trial_in.starttime, trial_in.endtime]);
	%From this apply gaussian filter to spike train for each electrode
	x = linspace(-sz/2, sz/2, sz);
	gaussFilter = exp(-x.^2/(2*sigma^2));
	gaussFilter = gaussFilter/sum(gaussFilter);

	for idx=1:nE
		gf = conv(binnedspikes(:,idx), gaussFilter, 'same');
		trial_in.nevrates(:,idx) = gf*trial_in.samplerate;
	end

	trial_out = trial_in;
end

