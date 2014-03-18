function plot_mua_rates_time(trial, fn)
        %plot_mua_rates_time       Plot electrode firing rate activity as a heatmap for given trial. Output to gif file
        %              named specified.
        %
        % Usage:
        %                       plot_mua_rates_time(trial,fn)
        %
        % Input:
        %                       trial = a trial structure from import_trials
        %                       fn = output file name for plot
        %
        % Examples:
        %                       fn = './worksheets/diagnostics/plots/test_mua_rates_time.eps';
        %                       trials = import_trials('./testdata/Spanky_2013-01-17-1325.mat');
	%			plot_mua_rates_time(trials(117), fn);

    close all;
	%fig = figure('visible', 'off');

    if (nargin < 2)
        throw(MException('Argin:MoreExpected', 'More input arguments expected'));
    end

	if ~isfield(trial, 'binnedspikes')
		trial = import_spikes(trial);
	end

	zaxis = [0 max(max(trial.nevrates(:,1:128)))];
   	image(trial.nevrates', 'CDataMapping', 'scaled');
   	caxis(zaxis);
   	set(gca,'Zlim',zaxis,'Ztick',zaxis, 'NextPlot', 'replacechildren');
   	xlabel('time (s)');
   	ylabel('channel');
    title('Firing rates')
    colorbar
	saveplot(gcf, fn);
end
