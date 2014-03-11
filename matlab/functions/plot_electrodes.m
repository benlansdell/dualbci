function plot_electrodes(trial, fn)
        %plot_electrodes       Plot electrode activity as a heatmap for all times specified. Output to gif file
        %              named specified.
        %
        % Usage:
        %                       plot_electrodes(trial,fn)
        %
        % Input:
        %                       trial = a trial structure from import_trials
        %                       fn = output file name for plot
        %
        % Examples:
        %                       fn = './tmp/test_electrodes_plot.gif';
        %                       trials = import_trials('./testdata/Spanky_2013-01-17-1325.mat');
	%			plot_electrodes(trials(117), fn);

        close all;
	fig = figure('visible', 'off');

        if (nargin < 2)
                throw(MException('Argin:MoreExpected', 'More input arguments expected'));
        end

	if ~isfield(trial, 'nevspikes')
		trial = import_spikes(trial);
	end

	%Plot animation of trial at 60Hz
        frames = 1:length(trial.cursor);
        %frames = 1:10;
	framen = 1;
	zaxis = [min(min(trial.nevrates)) max(max(trial.nevrates))];
        for i=frames
        	clf(fig);
        	sol=reshape(trial.nevrates(1:96,i),8,12);
		image(sol, 'CDataMapping', 'scaled');
        	caxis(zaxis);
        	set(gca,'Zlim',zaxis,'Ztick',zaxis, 'NextPlot', 'replacechildren');
        	xlabel(['time = ' num2str(trial.times(i))]);
        	%ylabel(names{1});
        	colorbar, drawnow
		title('Electrode firing rates')

		%Write file to temp png
	        saveplot(gcf, [fn '_tmp.png'], 'png');
	        %plotmult(gcf, [fn '_tmp.png'], nsubplots, 'png', [6 4]);

	        [X, m] = imread([fn '_tmp.png'], 'png');
 	        f = im2frame(X, m);
         	if i == 1
            		[im, map] = rgb2ind(f.cdata, 256, 'nodither');
            		im(1,1,1,length(frames)) = 0;
            	else
            		im(:,:,1,framen) = rgb2ind(f.cdata, map, 'nodither');
        	end
        	framen = framen + 1;
        end
	%Save to file
	imwrite(im, map, [fn], 'gif', 'DelayTime', 0, 'LoopCount', inf);
	%Remove tmp file if it exists
	delete([fn '_tmp.png']);
end
