%k-fold cross validate models

%Fit model to (k-1)/k of data, and compute deviance, encoding and decoding performance on other 1/k data
fold = 10;
fn_out = './worksheets/11_14_2014/data.mat';

%Preprocess data (timebis of 2ms, smooth trajectories)
%Apply no offset
nevfile = './testdata/20130117SpankyUtah001.nev';
binsize = 0.002;
dt_sp = binsize;
dt_pos = 0.2;
offset = 0.0;
threshold = 5;
processed = preprocess_spline(nevfile, binsize, threshold, offset);
nU = length(processed.unitnames);
N = size(processed.rates, 1);

alpha = 0.05;

const = 'on';
nK_sp = 100; 
nK_pos = 4;
data = filters_sp_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);

models_MSP = {};
testdata = {};
devs_MSP = zeros(fold, nU);
encodingccs = zeros(fold, nU);
decodingccs = zeros(fold, nU);

%Do training
display(['Training ' num2str(fold) '-fold cross validation models on datasets'])
for i = 1:fold
	display(['Generating data for fold: ' num2str(i)])
	%Generate training and test data
	[traindata, testd] = makefolds(data, i, fold);
	testdata{i} = testd;
	%Fit to training data
	display('Fitting model')
	models_MSP{i} = MLE_glmfit(traindata, const);
end
display('Done.')

%mid way save
save(fn_out, 'models_MSP', 'testdata', 'processed');

%Do testing
display(['Testing ' num2str(fold) '-fold cross validation models on datasets'])
for i = 1:fold
	display(['Fold: ' num2str(i)])
	order = 1;
	[F, Q, mu] = fit_AR_LS(testdata{i}.torque, order);
	%Deviance of test data
	devs_MSP(i,:) = deviance(models_MSP{i}, testdata{i});
	%Encoding of test data
	display('Estimating encoding performance')
	encodingccs(i,:) = plot_predictions(models_MSP{i}, testdata{i}, processed);
	%Decoding of test data
	display('Estimating decoding performance')
	output = glm_decode(processed, testdata{i}, models_MSP{i}, F, Q, mu, fn_out);
	%Compute GoF via time rescaling
	gof_timerescale(models_MSP{i}, testdata{i}, processed, alpha, fn_out);
end
display('Done.')

%Summarize all runs
%Plot deviance



%save fitted models for later use
save('./worksheets/11_14_2014/crossval_fittedmodels.mat', 'models_MSP', 'devs_MSP', 'encodingccs', 'decodingccs');