function electrode_check(fn, plotfn)
        %plotsol       Imports all trials from a .mat file output from LabVIEW as a list of structures.
        %               By default will not import all raw recording data, or torque data.
        %               Use import_torque, import_spikes, and import_raw to append those data to these structures
        %
        % Usage:
        %                       import_trials(fn)
        %
        % Input:
        %                       fn = input LabVIEW output filename
        %
        % Examples:
        %                       electrode_check('./data/20130117SpankyUtah001.nev', ./worksheets/diagnostics/elec_check_');

	%Check all electrodes are good...
	%Load nev file and plot number of all 144 channels activity
	NEV = openNEV(fn);

	nelectrodes = 144;
	samplerate = 30000;
	duration = (single(NEV.Data.Spikes.TimeStamp(1))-single(NEV.Data.Spikes.TimeStamp(end)))/samplerate;

	electrode = zeros(nelectrodes,1);
	
	for e=NEV.Data.Spikes.Electrode
		electrode(e) = electrode(e) + 1;
	end

	plot(1:nelectrodes, electrode);
	xlabel('Electrode #');
	ylabel('# spikes');
	title(['number of spikes in ' fn ' over ' num2str(duration) ' of recording'])
	saveplot(gcf, plotfn);
end
