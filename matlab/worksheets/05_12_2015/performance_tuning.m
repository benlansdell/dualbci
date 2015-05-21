files = {'20131031SpankyUtah014', '20131031SpankyUtah016', '20131205SpankyUtah001', '20131205SpankyUtah005' ...
			'20131104SpankyUtah001','20131104SpankyUtah005','20131011SpankyUtah001','20131011SpankyUtah004'};
matfile = {'Spanky_2013-12-05-1400.mat', 'Spanky_2013-11-04-1321.mat', 'Spanky_2013-10-31-1322.mat', 'Spanky_2013-10-11-1440.mat', 'Spanky_2013-09-30-1342.mat'};

files = {'20131031SpankyUtah014'};
matfile = {'Spanky_2013-10-31-1322.mat'};

%Read in spreadsheet
performance = xlsread('./TuningResults1.xlsx')

for idx = 1:length(files);
	%Load test preprocessed data
	nevfile = files{idx};
	nevpath = ['./testdata/', nevfile, '.nev']
	binsize = 1/60;
	offset = 0.0;
	threshold = 5;
	processed = preprocess_spline(nevpath, binsize, threshold, offset);
	%Modify so that it's shorter, then save it
	processed.binnedspikes = processed.binnedspikes(1:10800,:);
	processed.rates = processed.rates(1:10800,:);              
	processed.torque = processed.torque(1:10800,:);
	processed.dtorque = processed.dtorque(1:10800,:);
	processed.ddtorque = processed.ddtorque(1:10800,:);
	%save('./testdata/test_preprocess_brain_spline_60hz_short24.mat', 'processed');
	
	%Run with position filters
	const = 'on';
	pval = 0.001;
	nK_sp = 6;
	nK_pos = 6;
	fn_out = ['./worksheets/05_12_2015/plots/GC_' nevfile '.eps'];
	data = filters_sp_pos_network(processed, nK_sp, nK_pos);
	[GCdevMP, GCpvalMP, GCsigMP] = granger(processed, data, fn_out, pval);
	%Save results
	clear processed;
	save(['./worksheets/05_12_2015/plots/GC_' nevfile '.mat']);
end

for idx = 1:length(files)
	nevfile = files{idx};
	nevpath = ['./testdata/', nevfile, '.nev']
	tuningmatrix = NEVFastAnalysis_NPMK_func(nevpath);
	load(['./worksheets/05_12_2015/plots/GC_' nevfile '.mat']);
	tuningmatrix = sortrows(tuningmatrix);
	nU = size(tuningmatrix, 1);

	fn_out = ['./worksheets/05_12_2015/plots/GC_' nevfile '_wTuning.eps'];
%	unitnames{1} = 'FE'; 
%	unitnames{2} = 'RU';
	unitnames = {};
	for j = 1:nU
		unitnames{j} = [num2str(tuningmatrix(j,1)) '.' num2str(tuningmatrix(j,2))];
		if abs(tuningmatrix(j, 5)) > threshold;
			unitnames{j} = ['*', unitnames{j}];
		end
		if abs(tuningmatrix(j, 6)) > threshold;
			unitnames{j} = ['+', unitnames{j}];
		end
	end
	%GCdevwTuning = [tuningmatrix(:,5:6), GCdevMP];
	%GCdevwTuning = [GCdevwTuning; zeros(2,2), tuningmatrix(:,5:6)'];

	%Plot with Tuning scores
	clf
	colormap(bone);
	imagesc(GCdevMP)
	title('Change in deviance')
	ylabel('Unit')
	xlabel('Unit')
	set(gca,'XTick',1:(nU));
	set(gca,'YTick',1:(nU));
	set(gca,'XTickLabel',unitnames);
	set(gca,'YTickLabel',unitnames);
	rotateXLabels(gca, 90);
	colorbar
	saveplot(gcf, fn_out, 'eps', [6 9])
end