%k-fold cross validate models
%Fit model to (k-1)/k of data, and compute deviance, encoding and decoding performance on other 1/k data
fold = 10;
%Change this locally to change the folds that are computed
folds = 1:fold;
lfolds = length(folds);

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
devs_MSP = zeros(lfolds, nU);
encodingccs = zeros(lfolds, nU);
decodingccs = zeros(lfolds, nU);


%Do training
display(['Training ' num2str(fold) '-fold cross validation models on datasets'])
for idx = 1:lfolds
	i = folds(idx);
	display(['Generating data for fold: ' num2str(i)])
	%Generate training and test data
	[traindata, testd] = makefolds(data, i, fold);
	testdata{idx} = testd;
	%Fit to training data
	display('Fitting model')
	models_MSP{idx} = MLE_glmfit(traindata, const);
end
display('Done.')

%mid way save
save(fn_out, 'models_MSP', 'testdata', 'processed', 'fold', 'folds', 'lfolds', '-v7.3');

%Do testing
display(['Testing ' num2str(fold) '-fold cross validation models on datasets'])
for i = 1:lfolds
	display(['Fold: ' num2str(folds(i))])
	%Deviance of test data
	display('Computing deviance of model')
	devs_MSP(i,:) = deviance(models_MSP{i}, testdata{i});
	%Encoding of test data
	display('Estimating encoding performance')
	fn_out = ['./worksheets/11_14_2014/plots/encoding_fold_' num2str(folds(i))];
	encodingccs(i,:) = plot_predictions(models_MSP{i}, testdata{i}, processed, fn_out);
	%Decoding of test data
	%display('Estimating decoding performance')
	%order = 1;
	%[F, Q, mu] = fit_AR_LS(testdata{i}.torque, order);
	%output = glm_decode(processed, testdata{i}, models_MSP{i}, F, Q, mu, fn_out);
end
display('Done.')

display('Do one run where the test code overlaps with the training code to check for overfitting...')
%Deviance of test data
devs_MSP_overlap = deviance(models_MSP{1}, testdata{2});
%Encoding of test data
display('Estimating encoding performance')
fn_out = ['./worksheets/11_14_2014/plots/encoding_overlap_' num2str(folds(i))];
encodingccs_overlap = plot_predictions(models_MSP{1}, testdata{2}, processed, fn_out);
%Decoding of test data
%display('Estimating decoding performance')
%order = 1;
%[F, Q, mu] = fit_AR_LS(testdata{i}.torque, order);
%output = glm_decode(processed, testdata{i}, models_MSP{i}, F, Q, mu, fn_out);
display('Done.')

%%%%%%%%%%%%%%%%%%%%%%%%
%If training was run on different computers then consolidate here before running the plotting code
%%%%%%%%%%%%%%%%%%%%%%%%

load('./worksheets/11_14_2014/crossval_fittedmodels_latte.mat');
load('./worksheets/11_14_2014/data_latte.mat');
devs_MSP_latte = devs_MSP;
devs_MSP_overlap_latte = devs_MSP_overlap;
encodingccs_latte = encodingccs;
encodingccs_overlap_latte = encodingccs_overlap;
models_MSP_latte = models_MSP;

load('./worksheets/11_14_2014/crossval_fittedmodels_mocha.mat');
devs_MSP_mocha = devs_MSP;
devs_MSP_overlap_mocha = devs_MSP_overlap;
encodingccs_mocha = encodingccs;
encodingccs_overlap_mocha = encodingccs_overlap;
models_MSP_mocha = models_MSP;

devs_MSP = [devs_MSP_latte; devs_MSP_mocha];
devs_MSP_overlap = devs_MSP_overlap_latte;
encodingccs = [encodingccs_latte; encodingccs_mocha];
encodingccs_overlap = encodingccs_overlap_latte;
models_MSP = [models_MSP_latte, models_MSP_mocha];

%Summarize all runs
%Plot deviance
clf
boxplot(devs_MSP, 'labels', processed.unitnames, 'labelorientation', 'inline');
xlabel('Unit')
ylabel('Deviance')
set(gca,'XTick',1:length(processed.unitnames));
hold on
plot(devs_MSP_overlap, 'k*');
hold off
saveplot(gcf, './worksheets/11_14_2014/plots/deviance_crossval.eps');
%Plot encoding performance

clf
boxplot(encodingccs, 'labels', processed.unitnames, 'labelorientation', 'inline');
xlabel('Unit')
set(gca,'XTick',1:length(processed.unitnames));
ylabel('Correlation')
hold on
plot(encodingccs_overlap, 'k*');
hold off
saveplot(gcf, './worksheets/11_14_2014/plots/encoding_crossval.eps');
%Plot decoding performance

%save fitted models for later use
save('./worksheets/11_14_2014/crossval_fittedmodels.mat', 'models_MSP', 'devs_MSP', 'encodingccs_overlap', 'devs_MSP_overlap', 'encodingccs', 'decodingccs', '-v7.3');

