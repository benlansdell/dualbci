function trials = import_trials(fn)
        %plotsol       Imports all trials from a .mat file output from LabVIEW as a list of structures.
	%		By default will not import all raw recording data, or torque data.
	%		Use import_torque, import_spikes, and import_raw to append those data to trial structure
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



end
