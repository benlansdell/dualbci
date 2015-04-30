function plotEachDimVsTime_colmeantarget(seq, xspec, binWidth, targ, varargin)
%
% plotEachDimVsTime(seq, xspec, binWidth, ...)
%
% Plot each state dimension versus time in a separate panel.
%
% Color traces by octant which trial (target - start) pos lies in 
%
% INPUTS:
%
% seq       - data structure containing extracted trajectories
% xspec     - field name of trajectories in 'seq' to be plotted 
%             (e.g., 'xorth' or 'xsm')
% binWidth  - spike bin width used when fitting model
% octs      - vector indicating which octant trial lies in
%
% OPTIONAL ARGUMENTS:
%
% nPlotMax  - maximum number of trials to plot (default: 20)
% redTrials - vector of trialIds whose trajectories are plotted in red
%             (default: [])
% nCols     - number of subplot columns (default: 4)
%
% @ 2009 Byron Yu -- byronyu@stanford.edu

	nPlotMax  = 20;
	redTrials = [];
	nCols     = 4;
	assignopts(who, varargin);
	alph = 0.5;

	f = figure;
	pos = get(gcf, 'position');
	set(f, 'position', [pos(1) pos(2) 2*pos(3) pos(4)]);

	Xall = [seq.(xspec)];
	xMax = ceil(10 * max(abs(Xall(:)))) / 10; % round max value to next highest 1e-1
	
	Tmax    = max([seq.T]);  
	xtkStep = ceil(Tmax/25)*5;
	xtk     = 1:xtkStep:Tmax;
	xtkl    = 0:(xtkStep*binWidth):(Tmax-1)*binWidth;
	ytk     = [-xMax 0 xMax];

	nC = size(seq(1).xorth,1);
	nT = min(length(seq), nPlotMax);

	nRows   = ceil(size(Xall, 1) / nCols);
	cm = colormap(jet(8));;

	%Find min trial length
	for n = 1:min(length(seq), nPlotMax)
		tlen(n) = seq(n).T;
	end
	Tmin = min(tlen);
	for idx = 1:8
		dat_trunc{idx} = zeros(nT, nC, Tmin);
	end

	%Plot each trial
	for n = 1:nT
		%Projected activity
		dat = seq(n).(xspec);
		target = targ(n);
		dat_trunc{target}(n,:,:) = dat(:,1:Tmin);
	end

	for idx = 1:8
		dat_mean{idx} = squeeze(mean(dat_trunc{idx},1));
		dat_std{idx} = squeeze(std(dat_trunc{idx},1));
	end

	%On each set of axes
	for k = 1:nC
		for j = 1:8
			subplot((nRows+1), nCols, k);
			hold on;
			lw = 0.05;
			col = cm(j,:); 
			%plot(1:Tmin, dat_mean{j}(k,:), 'linewidth', lw, 'color', col);
			%Plot mean plus/minus std w transparency
			plotmeanstd(gcf, 1:Tmin, dat_mean{j}(k,:), dat_std{j}(k,:), col, alph);
		end
	end

	subplot((nRows+1), nCols, (nRows+1)*nCols)
	for k = 1:8
		col = cm(k, :);
		hold on
		plot(zeros(1,8), 'linewidth', lw, 'color', col);
	end
	legend('1', '2', '3', '4', '5', '6', '7', '8');

	%Plot axes
	%for k = 1:size(dat,1)
	%  h = subplot(nRows, nCols, k);
	%  axis([1 Tmax 1.1*min(ytk) 1.1*max(ytk)]);
	%
	%  if isequal(xspec, 'xorth')
	%    str = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$',k);
	%  else
	%    str = sprintf('$${\\mathbf x}_{%d,:}$$',k);
	%  end
	%  %title(str, 'interpreter', 'latex', 'fontsize', 16);
	%      
	%  set(h, 'xtick', xtk, 'xticklabel', xtkl);
	%  set(h, 'ytick', ytk, 'yticklabel', ytk);
	%  xlabel('Time (ms)');
	%end
