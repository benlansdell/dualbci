function plot_spikes_muas_nev(nevfile, fn)
        %plot_spikes_muas_nev       Plot spikes for each electrode as a heatmap for given nev file. Output to eps file
        %
        % Usage:
        %                       plot_spikes_muas_nev(nev,fn)
        %
        % Input:
        %                       nev = input nev file
        %                       fn = output file name for plot
        %
        % Examples:
        %                       fn = './worksheets/diagnostics/plots/test_total_spikes_muas_nev.eps';
        %                       nevfile = './testdata/20130117SpankyUtah005.nev';
	%			plot_spikes_muas_nev(nevfile, fn);
        close all;
	fig = figure;
        if (nargin < 2)
                throw(MException('Argin:MoreExpected', 'More input arguments expected'));
        end

        binsize = 1;
        threshold = 5;
        offset = 0;
        pre = preprocess(nevfile, binsize, threshold, offset);

	zaxis = [0 max(max(pre.binnedspikes))];
        image(pre.binnedspikes', 'CDataMapping', 'scaled');
        caxis(zaxis);
        set(gca,'Zlim',zaxis,'Ztick',zaxis, 'NextPlot', 'replacechildren');
        %ylabel(names{1});
        hcb = colorbar;
        drawnow;
	title(['activity for ' nevfile])
        ylabel('Unit #, total spikes')
        xlabel('time (s)')
        colorTitleHandle = get(hcb,'Title');
        titleString = 'spikes/s';
        set(colorTitleHandle ,'String',titleString);
        set(gca, 'YTick', 1:length(pre.unitnames));
        labels = {};
        for idx=1:length(pre.unitnames)
                nspikes = sum(pre.binnedspikes(:,idx));
                labels{idx} = [pre.unitnames{idx} '(' num2str(nspikes) ')'];
        end

        set(gca, 'YTickLabel', labels);
	%Write file
	saveplot(gcf, fn);
end
