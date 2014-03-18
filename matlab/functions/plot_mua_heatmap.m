function plot_mua_heatmap(trial, fn)
	%plot_mua_heatmap       Plot the spiking activity as a heatmap for each electrode during given trial. If import_spikes hasn't
	%					been run on the trial structure, then function will run import_spikes, adding spiking information 
	%			
	%
	% Usage:
	%					plot_mua_heatmap(trial,fn)
	%
	% Input:
	%                       trial = a trial structure from import_trials
	%                       fn = output file name for plot
	%
	% Examples:
	%                       %plot voltage and A fields
	%                       fn = './worksheets/diagnostics/plots/test_mua_heatmap.gif';
	%                       trials = import_trials('Spanky_2013-01-17-1325.mat');
	%						plot_mua_heatmap(trials(117), fn);
	close all;
	fig = figure('visible', 'off');
	if (nargin < 2)
			throw(MException('Argin:MoreExpected', 'More input arguments expected'));
	end
	%Load spike info, if not already present
	if ~isfield(trial, 'binnedspikes')
		trial = import_spikes(trial);
	end
	nsubplots=2;
	%Plot animation of trial at 60Hz
	 frames = 1:length(trial.cursor);
	%Test run
	%frames = 1:10;
	framen = 1;
	%Reshape to square array for plotting
	shape = size(trial.binnedspikes);
	bcimap = zeros(144,1);
	bcimap(floor(trial.electrodes)) = 1;
	bcimap = reshape(bcimap, 12, 12);
	spks = reshape(trial.binnedspikes, shape(1), 12, 12);
	spikes = permute(spks, [2 3 1]);
	zaxis = [min(min(trial.binnedspikes)), max(max(trial.binnedspikes))];
	for i=frames
		clf(fig);
		subplot(1,nsubplots,1);
		%Highlight BCI electrodes, other electrode info
		image(bcimap, 'CDataMapping', 'scaled');
		title('Electrodes used for BCI');
		%Plot spikes, not smoothed
		subplot(1,nsubplots, 2);
		image(spikes(:,:,i), 'CDataMapping', 'scaled');
		caxis(zaxis);
		set(gca,'Zlim',zaxis,'Ztick',zaxis, 'NextPlot', 'replacechildren');
		colorbar;
		xlabel('Electrode'); zlabel('MUA spike count');
		title(['MUA spike data. time=' num2str(trial.times(i))])
		%Write file to temp png
		plotmult(gcf, [fn '_tmp.png'], nsubplots, 'png', [6 4]);
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
	%Delete tmp file if it exists
	delete([fn '_tmp.png']);
end
