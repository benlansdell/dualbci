function raster_mua(trial, fn)
    %raster_mua       Plot electrodes firing spiking activity as a raster plot for specified trial. Output to eps file
    %
    % Usage:
    %                       raster_mua(trial,fn)
    %
    % Input:
    %                       trial = a trial structure from import_trials
    %                       fn = output file name for plot
    %
    % Examples:
    %                       fn = './worksheets/diagnostics/plots/test_raster.eps';
    %                       trials = import_trials('./testdata/Spanky_2013-01-17-1325.mat');
    %                       raster_mua(trials(117), fn);
    close all;

    if (nargin < 2)
            throw(MException('Argin:MoreExpected', 'More input arguments expected'));
    end

    if ~isfield(trial, 'binnedspikes')
	    trial = import_spikes(trial);
    end

    %convert spike times to correct format
    spikes = cell(length(trial.spikemuas),1);
    for idx = 1:length(trial.spikemuas)
        if length(trial.spikemuas(idx).times)>0
            spikes{idx} = trial.spikemuas(idx).times';
        else
            spikes{idx} = [0];    
        end
    end
    spikes
    plotSpikeRaster(spikes, 'PlotType', 'vertline');
    xlim([trial.starttime, trial.endtime]);
    xlabel('time (s)');
    ylabel('channel');
    saveplot(gcf, fn);

end
