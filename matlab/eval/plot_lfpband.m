function plot_lfpband(trial, bands, chans, fn_out)
    %plot_lfpband       Plot activity in specified bands for a trial. Output to eps file
    %
    % Usage:
    %                       plot_lfpband(trial,bands,chans,fn_out)
    %
    % Input:
    %                       trial = trial structure output by import_trials.m
    %                       bands = {[fmin1 fmax1], [fmin2 fmax2], ..} specifies the freq bands to compute mean activity for
    %                       fn_out = output file name for plot
    %                       chans = channels to compute LFPs for
    %
    % Examples:
    %                       fn_out = './worksheets/diagnostics/plots/test_lfpband_trial';
    %                       bands = {[0 5], [10 40], [45 65], [70 200], [200 400]};
    %                       chans = 1:2;
    %                       trials = import_trials('./testdata/Spanky_2013-01-17-1325.mat');
    %		                plot_lfpband(trials(117), bands, chans, fn_out);
    close all;

    if (nargin < 3)
            throw(MException('Argin:MoreExpected', 'More input arguments expected'));
    end

    if ~isfield(trial,'ns3data')
        trial = import_raw(trial);
    end

    params.Fs = trial.ns3samplerate;
    movingwin = [1 0.01];
    params.fpass = [0 400];
    params.tapers = [5 9];
    params.trialave = 0;
    params.err = 0;
    lfp = [];
    nsxelectrodes = trial.ns3data_raw(chans,:);
    cmap = colormap(hsv(length(bands)));

    %Compute LFP data for all channels
    [S,t,f] =  mtspecgramc(nsxelectrodes', movingwin, params);  
    %Split into frequency bands
    for idx = 1:length(bands)
        figure;
        idxb = (f>bands{idx}(1)) & (f<bands{idx}(2));
        Sb = S(:,idxb,:);
        %Average spectral content over band to get one feature per band per timestep per electrode
        lfp = squeeze(mean(Sb,2));
        %Normalize to z-scores?
        %Plot each bands activity for each channel as a heatmap
        zaxis = [min(min(lfp)), max(max(lfp))];
        image(lfp', 'CDataMapping', 'scaled');
        caxis(zaxis);
        colorbar;
        set(gca,'Zlim',zaxis,'Ztick',zaxis);
        xlabel('time (s)');
        ylabel('channel');
        title(['Activity for band ' num2str(bands{idx}(1)) ' to ' num2str(bands{idx}(2))]);
    end
end
