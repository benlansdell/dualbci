function trial_out = import_spikes(trial_in)
	%import_spikes		Imports spike information from .nev file corresponding to trial loaded from import_trials
	%
	% Usage:
	%					import_spikes(trial_in)
	%
	% Input:
	%					trial_in = trial structure from import_trials
	%
	% Output:
	%					trial_out = trial structure from import_trials with binnedspikes appended:
	%					binnedspikes is nE*nT array listing spike activity for each electrode
	%					in array of nE electrodes, at each time step, sampled at 60Hz, of nT timesteps
	%
	% Examples:
	%					trials = import_trials('Spanky_2013-01-17-1325.mat');
	%					trial = import_spikes(trials(117));

	if (length(trial_in.nevfile) > 0)
		trial_in.nevfile;
		%NEV=openNEV(['./blackrock/' trial_in.nevfile]);
		NEV=openNEV(['./testdata/' trial_in.nevfile]);
	else
		trial_out = trial_in;
		return;
	end

	nevsamplerate = NEV.MetaTags.TimeRes;
	labviewsamplerate = trial_in.samplerate;
	nE = length(NEV.MetaTags.ChannelID);
	flank = trial_in.flank;
	
	%get spike times in seconds
	spiketimes = double(NEV.Data.Spikes.TimeStamp)/nevsamplerate + double(trial_in.offset);
	nT = floor(trial_in.duration*labviewsamplerate)+1;
	elecs = cell(1,nE);
	trial_in.spikemuas = struct('times', elecs);
	for idx=1:nE
		trial_in.spikemuas(idx).times = [0];
		trial_in.spikemuas_flank(idx).times = [0];
	end
	for i=1:length(spiketimes)
		if (spiketimes(i) > trial_in.starttime) & (spiketimes(i) < trial_in.endtime)
			E = NEV.Data.Spikes.Electrode(i);
			trial_in.spikemuas(E).times = [trial_in.spikemuas(E).times; spiketimes(i)];
		end
		%Import also flanking data for computing tunings with lag
				if (spiketimes(i) > trial_in.starttime-flank) & (spiketimes(i) < trial_in.endtime+flank)
						E = NEV.Data.Spikes.Electrode(i);
						trial_in.spikemuas_flank(E).times = [trial_in.spikemuas_flank(E).times; spiketimes(i)];
				end
	end

	%Bin spikes and estimate the firing rate
	trial_in.binnedspikes = binspikes(trial_in.spikemuas, labviewsamplerate, [trial_in.starttime, trial_in.endtime]);
	trial_in.binnedspikes_flank = binspikes(trial_in.spikemuas_flank, labviewsamplerate, [trial_in.starttime-flank,trial_in.endtime+flank]);
	trial_in = gauss_rates(trial_in);
	trial_out = trial_in;
end
