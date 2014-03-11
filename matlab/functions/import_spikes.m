function trial_out = import_spikes(trial_in, timebin)
        %import_spikes       Imports spike information from .nev file corresponding to trial loaded from import_trials
        %
        % Usage:
        %                       import_spikes(trial_in)
        %
        % Input:
        %                       trial_in = trial structure from import_trials
	%			timebin = (optional, default=100ms) size of timebin over which to compute firing rates
	%
	% Output:
	%			trial_out = trial structure from import_trials with nevspikes appended:
	%				nevspikes is nE*nT array listing spike activity for each electrode
	%				in array of nE electrodes, at each time step, sampled at 60Hz, of nT timesteps
        %
        % Examples:
        %                       trials = import_trials('Spanky_2013-01-17-1325.mat');
	%			trial = import_spikes(trials(117));

	if length(trial_in.nevfile) > 0
		trial_in.nevfile
		NEV=openNEV(['./blackrock/' trial_in.nevfile]);
	else
		return;
	end

	nevsamplerate = NEV.MetaTags.TimeRes;
	%Find all spikes that occur within the 1/60s timebin
	labviewsamplerate = 60;
	span = 5;
	nE = length(NEV.MetaTags.ChannelID);
	
	spiketimes = single(NEV.Data.Spikes.TimeStamp)/nevsamplerate + trial_in.offset;
	nT = floor(trial_in.duration*labviewsamplerate)+1;
	trial_in.nevspikes = zeros(nE, nT);
	for i=1:length(spiketimes)
		if (spiketimes(i) > trial_in.starttime) & (spiketimes(i) < trial_in.endtime)
			T = floor((spiketimes(i)-trial_in.starttime)*labviewsamplerate)+1;
			E = NEV.Data.Spikes.Electrode(i);
			trial_in.nevspikes(E,T) = trial_in.nevspikes(E,T) + 1;
		end
	end
	%Compute a smoothed firing rate for each channel
	for i = 1:nE
		trial_in.nevrates(i,:) = smooth(trial_in.nevspikes(i,:), span);
	end

	%Compute firing rate from spikes
	%trial_in.nevrates = trial_	

	trial_out = trial_in;
end
