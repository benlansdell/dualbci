function plot_torque(trial, fn)
        %plot_torque       Plot torque measured by ns3 file and labview log file, ensuring they are consistent
        %
        % Usage:
        %                       plot_torque(trial,fn)
        %
        % Input:
        %                       trial = a trial structure from import_trials
        %                       fn = output file name for plot
        %
        % Examples:
        %                       fn = './worksheets/diagnostics/plots/torque_plot.eps';
        %                       trials = import_trials('./testdata/Spanky_2013-01-17-1325.mat');
	%			plot_torque(trials(117), fn);

        close all;
	%fig = figure('visible', 'off');

        if (nargin < 2)
                throw(MException('Argin:MoreExpected', 'More input arguments expected'));
        end

	if ~isfield(trial, 'ns3data')
		trial = import_raw(trial);
	end

	%Plot torque between labview and ns3 file
        %Torque in .nsx file
	plot(trial.times, trial.ns3data(138,:), trial.times, trial.ns3data(139,:));
	xlabel('time (s)');
	saveplot(gcf, [fn '_ns3.eps']);
        %Torque in .mat file
	plot(trial.times, trial.torque(:,1), trial.times, trial.torque(:,2));
	xlabel('time (s)')
	saveplot(gcf, [fn '_labview.eps']);
        %Cursor position
	plot(trial.times, trial.cursor(:,1), trial.times, trial.cursor(:,2));
	xlabel('time (s)')
	saveplot(gcf, [fn '_cursor.eps']);
        %Plot them against one another:
        %nsx torque x vs mat torque x
        plot(trial.ns3data(138,:), trial.torque(:,1));
        xlabel('torque axis 1 ns3 file');
        ylabel('torque axis 1 mat file');
        saveplot(gcf, [fn '_ns3_vs_mat.eps']);
        %nsx torque x vs cursor x
        plot(trial.ns3data(138,:), trial.cursor(:,1));
        xlabel('torque axis 1 ns3 file');
        ylabel('cursor position x mat file');
        saveplot(gcf, [fn '_ns3_vs_cur.eps']);
        %nsx torque x vs mat torque x
        plot(trial.ns3data(138,:), trial.torque(:,1));
        xlabel('torque axis 1 ns3 file');
        ylabel('torque axis 1 mat file');
        saveplot(gcf, [fn '_ns3_vs_mat.eps']);
        %Including lag bw cursor and torque
        lag = round(mean(trial.lag)*(10/16))
        %plot(trial.ns3data(138,1:(end-lag)), trial.torque((lag+1):end,1));
        plot(trial.times, trial.ns3data(138,:)/max(trial.ns3data(138,:)), trial.times, trial.torque(:,1)/max(trial.torque(:,1)));
        xlabel('time (s)');
        ylabel('torque x mat file vs labview');
        saveplot(gcf, [fn '_ns3_vs_cur.eps']);

end
