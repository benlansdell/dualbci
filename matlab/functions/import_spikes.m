function trial = import_spikes(trial)
        %import_spikes       Imports spike information from .nev file corresponding to trial loaded from import_trials
        %
        % Usage:
        %                       import_spikes(trial)
        %
        % Input:
        %                       trial = trial structure from import_trials
	%
	% Output:
	%			trial = trial structure from import_trials with nevspikes appended:
	%				nevspikes is nE*nT array listing spike activity for each electrode
	%				in array of nE electrodes, at each time step, sampled at 60Hz, of nT timesteps
        %
        % Examples:
        %                       trials = import_trials('Spanky_2013-01-17-1325.mat');
	%			trial = import_spikes(trials(117));

	if length(trial.nevfile) > 0
		NEV=openNEV(['./data/' trial.nevfile]);
	else
		return;
	end

	nevsamplerate = NEV.MetaTags.TimeRes;
	%Find all spikes that occur within the 1/60s timebin
	labviewsamplerate = 60;
	span = 5;
	nE = length(NEV.MetaTags.ChannelID);
	
	spiketimes = single(NEV.Data.Spikes.TimeStamp)/nevsamplerate + trial.offset;
	nT = floor(trial.duration*labviewsamplerate)+1;
	trial.nevspikes = zeros(nE, nT);
	for i=1:length(spiketimes)
		if (spiketimes(i) > trial.starttime) & (spiketimes(i) < trial.endtime)
			T = floor((spiketimes(i)-trial.starttime)*labviewsamplerate)+1;
			E = NEV.Data.Spikes.Electrode(i)
			trial.nevspikes(E,T) = trial.nevspikes(E,T) + 1;
		end
	end

	%Compute a smoothed firing rate for each channel
	for i = 1:nE
		trial.nevrates(i,:) = smooth(trial.nevspikes(i,:), span);
	end
end
