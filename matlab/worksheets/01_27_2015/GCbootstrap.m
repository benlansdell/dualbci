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
fn_out = './worksheets/01_27_2015/wholedataset.m';

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
		sample = randsample(nB,N,false);
		%Then, within each subset, compute deviance of a training and test set
		testdata = data;
		testdata.X = testdata.X(sample,:);
		testdata.y = testdata.y(:,sample);
		testdata.torque = testdata.torque(sample,:);
		testdata.dtorque = testdata.dtorque(sample,:);
		testdata.ddtorque = testdata.ddtorque(sample,:);
	    testdev = deviance_network(model, testdata);
	    traindevs(i,j,:) = traindev;
	    testdevs(i,j,:) = testdev;
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

%month = repmat({'jan' 'feb' 'mar' 'apr' 'may' 'jun' 'jul' 'aug' 'sep' 'oct' 'nov' 'dec'},1,2);
%simobs = [repmat({'sim'},1,12),repmat({'obs'},1,12)];
clf
subplot(3,1,1)
boxplot(summedtraindevs');
xlabel('N');
ylabel('deviance train')
subplot(3,1,2)
boxplot(summedtestdevs');
xlabel('N');
ylabel('deviance test')
subplot(3,1,3)
boxplot(summeddiff');
xlabel('N');
ylabel('log-likelihood b/w test and train')
fn_out = './worksheets/01_27_2015/plots/GCbootstrap_ll_boxplot.eps';
saveplot(gcf, fn_out, 'eps', [9 6]);

%Compute GC for different Ns. Just make the plots and see the difference...

%[GCdev, GCpval, GCsig] = granger(processed, data, fn_out, pval);
