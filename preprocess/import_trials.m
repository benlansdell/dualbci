function trials = import_trials(fn)
	%import_trials  Imports all trials from a .mat file output from LabVIEW as a list of structures.
	%		By default will not import all raw recording data, or torque data.
	%		Use import_torque, import_spikes, and import_raw to append those data to these structures
	%
	% Usage:
	%       import_trials(fn)
	%
	% Input:
	%       fn = input LabVIEW output filename
	%
	% Output:
	%		trials = list of structures containing the following fields:
	%			starttime, endtime
	%			cursorpos, error
	%			cursorstart, target
	%			success, valid
	%			electrodes, spikes
	%			nev_files, nsx_files
	%
	% Examples:
	%       trials = import_trials('./testdata/Spanky_2013-01-17-1325.mat');

	load(fn);
	n_trials = length(data.trials.time);

	%Debug
	%n_trials = 1;

	trials = [];
	for i=1:n_trials
		trial.duration = data.trials.duration(i);
		trial.duration = data.trials.duration(i);
		trial.endtime = data.trials.time(i);
		trial.starttime = trial.endtime-trial.duration;
		trial.success = data.trials.success(i);
		trial.valid = data.trials.valid(i);
		trial.duration = data.trials.duration(i);
		trial.cursorstart = data.trials.startPos(i,:);
		trial.target = data.trials.targetPos(i,:);
		trial.samplerate = data.sampleRate;
	
		withintrial = (data.stateHist.time < trial.endtime) & (data.stateHist.time > trial.starttime);
		trial.cursor = data.stateHist.cursor(withintrial,:);		
		trial.errors = data.stateHist.error(withintrial,:);		
		trial.velocity = data.stateHist.velocity(withintrial,:);
		trial.times = data.stateHist.time(withintrial);
		trial.spikes = data.stateHist.spikes(withintrial,:); 
		trial.torque = trial.spikes(:,5:7);
		trial.spikes = trial.spikes(:,1:4);
		trial.rates = data.stateHist.rates(withintrial,:);
		trial.lag = data.stateHist.lag(withintrial,:);

		trial.nevfile = '';
		trial.ns3file = '';
		trial.electrodes = [];
		trial.type = '';
		trial.offset = 0;
		trial.flank = 1;
		for j = 1:length(data.nev)
			nevdur = data.nev(j).DurationSec;
			nevoffset = single(data.nev(j).Toffset(1))/60;
			if (trial.starttime > nevoffset) & (trial.starttime < (nevoffset + nevdur))
				%Offset here is measured in seconds
				trial.offset = double(nevoffset);
				trial.nevfile = data.nev(j).nevfile;
				trial.ns3file = [trial.nevfile(1:end-3) 'ns3'];
				trial.electrodes = data.nev(j).chans;
				trial.type = data.nev(j).map;
			end
		end
		trials = [trials trial];
	end
end
