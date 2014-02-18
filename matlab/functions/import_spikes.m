function trial = import_spikes(trial)
        %plotsol       Imports all trials from a .mat file output from LabVIEW as a list of structures.
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

	openNEV(['./data/' trial.nevfile]);

	nevsamplerate = NEV.MetaTags.TimeRes;
	%Find all spikes that occur within the 1/60s timebin
	labviewsamplerate = 60;
	nE = length(NEV.MetaTags.ChannelID);
	
	spiketimes = NEV.Data.Spikes.TimeStamp/nevsamplerate + trial.offset;
	nT = floor(trial.duration*labviewsamplerate)+1;
	trial.nevspikes = zeros(nE, nT);
	for i=1:length(spiketimes)
		if (spiketimes(i) > trial.starttime) & (spiketimes(i) < trial.endtime)
			T = (spiketimes-trial.starttime)*labviewsamplerate;
			E = NEV.Data.Spikes.Electrode(i);
			trial.nevspikes(E,T) = trial.nevspikes(E,T) + 1;
		end
	end
end
