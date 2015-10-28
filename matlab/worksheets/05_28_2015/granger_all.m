files = {'20131031SpankyUtah014', '20131031SpankyUtah016', '20131205SpankyUtah001', '20131205SpankyUtah005' ...
			'20131104SpankyUtah001','20131104SpankyUtah005','20131011SpankyUtah001','20131011SpankyUtah004'};
matfile = {'Spanky_2013-12-05-1400.mat', 'Spanky_2013-11-04-1321.mat', 'Spanky_2013-10-31-1322.mat', 'Spanky_2013-10-11-1440.mat', 'Spanky_2013-09-30-1342.mat'};

%Test run

%Brain
files = {'20130930SpankyUtah015', '20131001SpankyUtah012'};
matfile = {'Spanky_2013-09-30-1342.mat', 'Spanky_2013-10-01-1346.mat'};

%Manual control
files = {'20131001SpankyUtah001'};
matfile = {'Spanky_2013-10-01-1346.mat'};

threshold = 5;
[data, txt] = xlsread('./exptinfo/Results1_excluding2d.xls');

for idx = 1:size(data,1)
	nevfilenms{idx} = [datestr(x2mdate(data(idx,1)),'yyyymmdd'), 'SpankyUtah' num2str(data(idx,4), '%03i') '.nev'];
end

for idx = 1:length(files);
	%Load test preprocessed data
	nevfile = files{idx};
	nevpath = ['./testdata/', nevfile, '.nev']
	binsize = 1/60;
	offset = 0.0;
	threshold = 5;

	nevidx = find(~cellfun('isempty', regexp([nevfile '.nev'], nevfilenms)))	
	BCunit1 = num2str(data(nevidx, 6), '%.1f');
	BCunit2 = num2str(data(nevidx, 7), '%.1f');

	processed = preprocess_spline(nevpath, binsize, threshold, offset);

	%Modify so that it's shorter, then save it
	%Three minutes worth of data
	processed.binnedspikes = processed.binnedspikes(1:10800,:);
	processed.rates = processed.rates(1:10800,:);              
	processed.torque = processed.torque(1:10800,:);
	processed.dtorque = processed.dtorque(1:10800,:);
	processed.ddtorque = processed.ddtorque(1:10800,:);
	%save('./testdata/test_preprocess_brain_spline_60hz_short24.mat', 'processed');
	%Truncate to top 15 units by firing rate
	[totalspikes, indices] = sort(sum(processed.binnedspikes,1), 2, 'descend')

	nU = 15; %size(processed.binnedspikes,2);
	top15 = indices(1:15);
	processed.binnedspikes = processed.binnedspikes(:,top15);
	processed.rates = processed.rates(:,top15);              
	processed.unitnames = processed.unitnames(top15);

	BCunit1idx = find(~cellfun('isempty', regexp(BCunit1, processed.unitnames)));
	BCunit2idx = find(~cellfun('isempty', regexp(BCunit2, processed.unitnames)));

	%Run with position filters
	const = 'on';
	pval = 0.001;
	nK_sp = 6;
	nK_pos = 6;
	fn_out = ['./worksheets/06_04_2015/plots/GC_' nevfile '.eps'];
	gdata = filters_sp_pos_network(processed, nK_sp, nK_pos);
	[GCdev, GCpval, GCsig] = granger(processed, gdata, fn_out, pval);
	%Save results
	results.firingrate = mean(processed.binnedspikes,1);
	results.unitnames = processed.unitnames;
	results.GCdev = GCdev;
	results.GCpval = GCpval;
	results.GCsig = GCsig;
	results.causaldensity = sum(sum(GCsig))/nU/(nU-1);

	results.BC1inSig = sum(GCsig(BCunit1idx,:));
	results.BC1outSig = sum(GCsig(:,BCunit1idx));
	results.BC2inSig = sum(GCsig(BCunit2idx,:));
	results.BC2outSig = sum(GCsig(:,BCunit2idx));

	results.BC1inDev = mean(GCdev(BCunit1idx,:));
	results.BC1outDev = mean(GCdev(:,BCunit1idx));
	results.BC2inDev = mean(GCdev(BCunit2idx,:));
	results.BC2outDev = mean(GCdev(:,BCunit2idx));

	results.BC1toBC2dev = GCdev(BCunit2idx, BCunit1idx);
	results.BC2toBC1dev = GCdev(BCunit1idx, BCunit2idx);

	GCwoutBCrows = GCdev;
	GCwoutBCrows([BCunit1idx, BCunit2idx],:) = [];
	results.BCinAveDev = mean(mean(GCwoutBCrows));
	results.BCinStdDev = std(reshape(GCwoutBCrows, 1, []));

	GCwoutBCcols = GCdev;
	GCwoutBCcols(:,[BCunit1idx, BCunit2idx]) = [];
	results.BCoutAveDev = mean(mean(GCwoutBCcols));

	clear processed;
	save(['./worksheets/06_04_2015/GC_' nevfile '.mat']);
end

