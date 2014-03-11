function plot_electrodes_nev(nev, fn)
        %plot_electrodes_nev       Plot total number of spikes for each electrode as a heatmap for given nev file. Output to eps file
        %
        % Usage:
        %                       plot_electrodes_nev(nev,fn)
        %
        % Input:
        %                       nev = input nev file
        %                       fn = output file name for plot
        %
        % Examples:
        %                       fn = './tmp/test_electrodes_plot.eps';
        %                       nev = './testdata/20130117SpankyUtah005.nev';
	%			plot_electrodes_nev(nev, fn);
        close all;
	%fig = figure('visible', 'off');
	fig = figure;

        if (nargin < 2)
                throw(MException('Argin:MoreExpected', 'More input arguments expected'));
        end

	NEV = openNEV(nev);
        nE = length(NEV.MetaTags.ChannelID);
        nS = length(NEV.Data.Spikes.TimeStamp);
	spikes = zeros(1,nE);
        for i=1:nS
                E = NEV.Data.Spikes.Electrode(i);
                spikes(E) = spikes(E) + 1;
        end

	zaxis = [0 max(spikes)];
        clf(fig);
	subplot(2,1,1);
        sol=reshape(spikes(1:96),12,8);
        image(sol, 'CDataMapping', 'scaled');
        caxis(zaxis);
        set(gca,'Zlim',zaxis,'Ztick',zaxis, 'NextPlot', 'replacechildren');
        %ylabel(names{1});
        colorbar, drawnow
	title('Spike counts for electrodes 1-96')

        subplot(2,1,2);
        sol=reshape(spikes(97:128),4,8);
        image(sol, 'CDataMapping', 'scaled');
        caxis(zaxis);
        set(gca,'Zlim',zaxis,'Ztick',zaxis, 'NextPlot', 'replacechildren');
        %ylabel(names{1});
        colorbar, drawnow
        title('Spike counts for electrodes 97-128')
	xlabel(['Recording duration: ' num2str(NEV.MetaTags.DataDurationSec) 's']);

	%Write file to temp png
	saveplot(gcf, fn);
end
