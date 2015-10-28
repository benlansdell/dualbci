%Granger analysis on weekly dataset
fn_out = './worksheets/01_17_2015/weeklynevs.mat';
conds = {'2D Manual Position'};
load(fn_out);
modeltype = 'sprc';

%For each condition and for each file in each condition:
for i = 1:length(conds)
	nspikes = [];
	devs = [];
	condition = conds{i};
	nevs = weeklynevs(condition);
	for j = 1:length(nevs)
		recinfo = nevs(j);
		%Get files
		nevfile = ['./blackrock/' recinfo.nevfile];
		matfile = ['./labview/' recinfo.matfile];
		currdate = strrep(recinfo.date, '/', '-');
		curr_fn = ['./worksheets/02_07_2015/glmresults/' modeltype '_' condition '_' currdate '.mat'];
		curr_fn = strrep(curr_fn, ' ', '_')
		if (exist(nevfile, 'file') == 2) & (exist(matfile, 'file') == 2)
			load(curr_fn);
			binsize = processed_mua.binsize;
			T = binsize*size(processed_mua.binnedspikes,1);
			%Extract number of spikes and the deviance of each unit
			nU = size(model.b_hat, 1);
			for u = 1:nU
				dev = model.dev{u};
				nS = sum(processed_mua.binnedspikes(:,u));
				nspikes = [nspikes(:); nS];
				devs = [devs(:); dev];
			end
		end
	end
	%Plot
	plot((nspikes/T), (devs), '.');
	xlabel('firing rate (Hz)');
	xlim([0 80])
	ylabel('deviance of MLE sprc model')
	title([num2str(T) ' seconds of recording'])
	saveplot(gcf, './worksheets/02_14_2015/plots/spikestofit_sprc.eps')
	ylim([0 1e80])
	xlim([0 10])
	saveplot(gcf, './worksheets/02_14_2015/plots/spikestofit_sprc_zoom1.eps')
	ylim([0 1e78])
	xlim([0 2])
	saveplot(gcf, './worksheets/02_14_2015/plots/spikestofit_sprc_zoom2.eps')
	ylim([0 1e76])
	xlim([0 1])
	saveplot(gcf, './worksheets/02_14_2015/plots/spikestofit_sprc_zoom3.eps')
	ylim([0 1e5])
	xlim([0 1])
	saveplot(gcf, './worksheets/02_14_2015/plots/spikestofit_sprc_zoom4.eps')
end