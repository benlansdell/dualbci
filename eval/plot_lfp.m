function plot_lfp(trial, chans, fn_out)
    %plot_lfp       Plot spectrogram of trial for specified channels. Output to eps file
    %
    % Usage:
    %                       plot_lfp(trial,chans,fn_out)
    %
    % Input:
    %                       trial = trial structure output by import_trials.m
    %                       chans = [chmin chmax] specifies the range of channels to read from
    %                       fn_out = output file name for plot
    %
    % Examples:
    %                       fn_out = './worksheets/diagnostics/plots/test_lfp_trial';
    %                       chans = [55 59];
    %                       trials = import_trials('./testdata/Spanky_2013-01-17-1325.mat');
    %		                plot_lfp(trials(117), chans, fn_out);

    close all;
    fig = figure;

    if (nargin < 3)
            throw(MException('Argin:MoreExpected', 'More input arguments expected'));
    end

    if ~isfield(trial,'ns3data')
        trial = import_raw(trial);
    end

    %Plot the data
    params.Fs = trial.ns3samplerate;
    movingwin = [0.2 0.01];
    params.fpass = [0 400];
    params.tapers = [5 9];
    params.trialave = 0;
    params.err = 0;
    for idx = 1:length(chans)
	    figure
	    %t_in.ns3data = t_in.ns3data(:,withintrial);
        [S,t,f] =  mtspecgramc(trial.ns3data_raw(chans(idx),:)', movingwin, params);
        plot_matrix(S,t,f);
        title(['Frequency content of channel ' num2str(chans(idx))]);
        %Write file
        saveplot(gcf, [fn_out '_ch_' num2str(chans(idx)) '.eps']);
    end
end
