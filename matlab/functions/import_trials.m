function trials = import_trials(fn)
        %plotsol       Imports all trials from a .mat file output from LabVIEW as a list of structures.
	%		By default will not import all raw recording data, or torque data.
	%		Use import_torque, import_spikes, and import_raw to append those data to these structures
        %
        % Usage:
        %                       import_trials(fn)
        %
        % Input:
        %                       fn = input LabVIEW output filename
	%
	% Output:
	%			trials = list of structures containing the following fields:
	%				starttime, endtime
	%				cursorpos, error
	%				cursorstart, target
	%				success, valid
	%				electrodes, spikes
	%				nev_files, nsx_files
        %
        % Examples:
        %                       trials = import_trials('Spanky_2013-01-17-1325.mat');

	load(fn);
	n_trials = length(data.trials.time)

	%Debug
	%n_trials = 1;

	trials = [];
	for i=1:n_trials
		trial.starttime = data.trials.time(i);
		trial.success = data.trials.success(i);
		trial.valid = data.trials.valid(i);
		trial.duration = data.trials.duration(i);
		trial.endtime = trial.starttime+trial.duration;
		trial.cursorstart = data.trials.startPos(i,:);
		trial.target = data.trials.targetPos(i,:);
	
		withintrial = (data.stateHist.time < trial.endtime) & (data.stateHist.time > trial.starttime);
		trial.cursor = data.stateHist.cursor(withintrial,:);		
		trial.errors = data.stateHist.error(withintrial,:);		
		trial.velocity = data.stateHist.velocity(withintrial,:);
		trial.times = data.stateHist.time(withintrial);
		
		trial.spikes = data.stateHist.spikes(withintrial); 
		trial.rates = data.stateHist.rates(withintrial);

		%determine the nev, nsx file to look in. extract electrodes, raw data
		%for each nev file in data.nev, check if the duration and offset mean that it contains the data we want
		trial.nevfile = '';
		trial.electrodes = [];
		for j = 1:length(data.nev)
			nevdur = data.nev(j).DurationSec;
			nevoffset = data.nev(j).Toffset(1)/60;
			if (trial.starttime > nevoffset) & (trial.starttime < (nevoffset + nevdur))
				trial.nevfile = data.nev(j).nevfile;
				trial.ns3file = [trial.nevfile(1:end-3) 'ns3'];
				trial.electrodes = data.nev(j).chans;
			end
		end
		trials = [trials trial];
	end
end
