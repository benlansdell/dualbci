%Load data
nevfile = './testdata/20130117SpankyUtah001.nev';
matfile = './testdata/Spanky_2013-01-17-1325.mat';
binsize = 1/60;
nK_pos = 6;
nK_sp = 6;
const = 'on';
offset = 0.0;
threshold = 5;
processed = preprocess_spline_lv(nevfile, matfile, binsize, threshold, offset);

%Chop down to 10 units so this can be easily applied to other recordings that have different # of units
nU = 10;
processed.binnedspikes = processed.binnedspikes(:,1:10);
processed.rates = processed.rates(:,1:10);
processed.unitnames = processed.unitnames(1:10);

%Chop down to 500s' worth
nB = floor(500/processed.binsize);
processed.binnedspikes = processed.binnedspikes(1:nB,:);
processed.rates = processed.rates(1:nB,:);
processed.cursor = processed.cursor(1:nB,:);
processed.dcursor = processed.dcursor(1:nB,:);
processed.ddcursor = processed.ddcursor(1:nB,:);
processed.torque = processed.torque(1:nB,:);
processed.dtorque = processed.dtorque(1:nB,:);
processed.ddtorque = processed.ddtorque(1:nB,:);

data = filters_sp_pos_network(processed, nK_sp, nK_pos);
fn_out = './worksheets/01_27_2015/wholedataset_replace_distinct.mat';

%update number of bins
nB = size(data.y,2);

%Save for later
save(fn_out);

%Make some subsamples vectors of different sizes...
Nfrac = [1 0.5 0.4 0.3 0.2 0.1];
K = 50;

testdevs = zeros(length(Ns), K, nU);
traindevs = zeros(length(Ns), K, nU);

Ns = floor(nB.*Nfrac);
%Generate K samples for each N, record deviance of training and test set
for i = 1:length(Ns)
	N = Ns(i);
	display(['Computing bootstrap estimates of model for N=' num2str(N)])
	%samples = zeros(K, N);
	for j = 1:K
		display(['Repeat number: ' num2str(j)])
		sample = randsample(nB,N,false);
		size(sample);
		%Then, within each subset, compute deviance of a training and test set
		traindata = data;
		traindata.X = traindata.X(sample,:);
		traindata.y = traindata.y(:,sample);
		traindata.torque = traindata.torque(sample,:);
		traindata.dtorque = traindata.dtorque(sample,:);
		traindata.ddtorque = traindata.ddtorque(sample,:);
		model = MLE_glmfit_network(traindata, const);
		%Randomly sample another test set and compute the deviance
	    traindev = deviance_network(model, traindata);
	    allpts = 1:nB;
	    notsampled = allpts(~ismember(allpts, sample));
	    %Take a subsample of these
	    if N < nB
	    	testsample = randsample(notsampled, N, false);
			%Then, within each subset, compute deviance of a training and test set
			testdata = data;
			testdata.X = testdata.X(testsample,:);
			testdata.y = testdata.y(:,testsample);
			testdata.torque = testdata.torque(testsample,:);
			testdata.dtorque = testdata.dtorque(testsample,:);
			testdata.ddtorque = testdata.ddtorque(testsample,:);
	    	testdev = deviance_network(model, testdata);
		    testdevs(i,j,:) = testdev;
		end
	    traindevs(i,j,:) = traindev;
	end
end

%Plot all this stuff
summedtraindevs = sum(traindevs, 3);
summedtestdevs = sum(testdevs, 3);
for i = 1:length(Ns)
	N = Ns(i);
	summedtraindevs(i,:) = summedtraindevs(i,:)/N;
	summedtestdevs(i,:) = summedtestdevs(i,:)/N;
end
summeddiff = summedtraindevs - summedtestdevs;
%Invert so shorter training times are to the left, not the right
summeddiff = flipud(summeddiff);
summedtraindevs = flipud(summedtraindevs);
summedtestdevs = flipud(summedtestdevs);

%month = repmat({'jan' 'feb' 'mar' 'apr' 'may' 'jun' 'jul' 'aug' 'sep' 'oct' 'nov' 'dec'},1,2);
%simobs = [repmat({'sim'},1,12),repmat({'obs'},1,12)];
clf
invNs = fliplr(Ns);
labels = {}; %cellfun(@num2str,num2cell(fliplr(Ns)*binsize),'uniformoutput',0);
for idx = 1:(length(Ns)-1)
	labels{idx} = [num2str(ceil(invNs(idx)*binsize))];% '\newline ' num2str(invNs(idx))];
end
subplot(3,1,1)
boxplot(summedtraindevs(1:end-1,:)', 'labels', labels);
ylabel('deviance train')
subplot(3,1,2)
boxplot(summedtestdevs(1:end-1,:)', 'labels', labels);
ylabel('deviance test')
subplot(3,1,3)
boxplot(summeddiff(1:end-1,:)', 'labels', labels);
xlabel('seconds of training');
ylabel('log-likelihood b/w test and train')
fn_out = './worksheets/01_27_2015/plots/GCbootstrap_replacement_distinct_ll_boxplot.eps';
saveplot(gcf, fn_out, 'eps', [6 15]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute GC for different Ns. Just make the plots and see the difference...%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nfrac = [0.5 0.4 0.3 0.2 0.1];
Ns = floor(nB.*Nfrac);
K = 20;
pval = 0.001;
fn_out = './worksheets/01_27_2015/plots/granger_sp_pos.eps';
GCdevs = {};
for i = 1:length(Ns)
	N = Ns(i);
	display(['Computing bootstrap estimates of model for N=' num2str(N)])
	%samples = zeros(K, N);
	for j = 1:K
		display(['Repeat number: ' num2str(j)])
		sample = randsample(nB,N,false);
		size(sample);
		%Then, within each subset, compute deviance of a training and test set
		traindata = data;
		traindata.X = traindata.X(sample,:);
		traindata.y = traindata.y(:,sample);
		traindata.torque = traindata.torque(sample,:);
		traindata.dtorque = traindata.dtorque(sample,:);
		traindata.ddtorque = traindata.ddtorque(sample,:);
		%model = MLE_glmfit_network(traindata, const);
		GCdev = granger(processed, traindata, fn_out, pval);
		GCdevs{i,j} = GCdev;
	end
end
fn_out = './worksheets/01_27_2015/gc_dev_variance.mat';
save(fn_out);

%Plot a boxplot of each 
GCdata = zeros(length(Ns), K, nU*nU);
for i = 1:length(Ns)
	for j = 1:K
		GCdata(i,j,:) = reshape(GCdevs{i,j},nU*nU,1);
	end
end
for i = 1:length(Ns)
	clf
	boxplot(squeeze(GCdata(i,:,:)));
	title([num2str(ceil(Ns(i)*binsize)) ' seconds of training'])
	fn_out = ['./worksheets/01_27_2015/plots/GCboxplots_' num2str(i) '.eps'];
	saveplot(gcf, fn_out, 'eps', [20 6]);
end


