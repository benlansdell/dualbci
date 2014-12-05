%Analysis of deviance scripts

%Preprocess data (timebis of 2ms, smooth trajectories)
%Apply no offset
nevfile = './testdata/20130117SpankyUtah001.nev';
labviewfile = './testdata/Spanky_2013-01-17-1325.mat';
fn_out = './worksheets/12_1_2014/plots/test';
binsize = 0.002;
dt_sp = binsize;
dt_pos = 0.2;
dt_tar = 0.2;
dt_vel = 0.1;
offset = 0.0;
threshold = 5;
processed = preprocess_spline_target(nevfile, labviewfile, binsize, threshold, offset);
nU = length(processed.unitnames);
N = size(processed.rates, 1);

%- MPS mean, spike history and position
%Do for a range of filter sizes
const = 'on';
nK_sp = 100; 
nK_poss = [5];
models_MSP = {};
for idx = 1:length(nK_poss)
	idx
	nK_pos = nK_poss(idx);
	data = filters_sp_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
	%Fit each of the above GLMs
	models_MSP{idx} = MLE_glmfit(data, const);
	%Make plots of each filter fitted, predictions of each unit, and record the deviance
	fn_out = ['./worksheets/12_1_2014/plots/AOD_MSP_nK_' num2str(nK_pos)];
	plot_filters(models_MSP{idx}, data, processed, fn_out);
	plot_predictions(models_MSP{idx}, data, processed, fn_out);
end

%- MSPT, spike history, position and target!
const = 'on';
nK_sp = 100; 
nK_poss = [5];
nK_tar = 1;
models_MSPT = {};
for idx = 1:length(nK_poss)
	idx
	nK_pos = nK_poss(idx);
	data = filters_sp_pos_target(processed, nK_sp, nK_pos, nK_tar, dt_sp, dt_pos, dt_tar);
	%Fit each of the above GLMs
	models_MSPT{idx} = MLE_glmfit(data, const);
	%Make plots of each filter fitted, predictions of each unit, and record the deviance
	fn_out = ['./worksheets/12_1_2014/plots/AOD_MSPT_nK_' num2str(nK_pos)];
	plot_filters(models_MSPT{idx}, data, processed, fn_out);
	plot_predictions(models_MSPT{idx}, data, processed, fn_out);
end

%save fitted models for later use
save('./worksheets/12_1_2014/AOD_fittedmodels.mat', 'models_MSP', 'models_MSPT')
load('./worksheets/12_1_2014/AOD_fittedmodels.mat', 'models_MSP', 'models_MSPT')

L = length(nK_poss);
%L = length(nK_vels);
nU = size(model_MS.b_hat,1);
csvMSP = zeros(nU, 5+6*L);
blank = cell(1,L-1);
%Fields to save (nested model MS vs MSP)
%Delta AIC (MSP), Delta BIC (MSP), correlation

MSP_headings = {'Name', 'Dev MSP', 'nMSP', 'Dev MSPT', 'nMSPT', 'Adj Dev MSPT', 'chi2', 'p-val', 'AIC', 'BIC'};
for idx = 1:nU
	%Unit name
	csvMSP(idx, 1) = str2num(processed.unitnames{idx});
	%Dev MSP
	csvMSP(idx, 2) = models_MSP{1}.dev{idx};
	%nMSP
	csvMSP(idx, 3) = size(models_MSP{1}.b_hat,2);
	%Dev MSPT
	csvMSP(idx, 4) = models_MSPT{1}.dev{idx};
	%nMSPT
	csvMSP(idx, 5) = size(models_MSPT{1}.b_hat, 2);
	%Adj Dev MSP
	%Adjusted since the number of data points in the MSP fits is different to the number of data points in the MS fits, which
	%Fewer data points means fewer deviances... which can make MSP fits better than they seem
	NMSP = models_MSP{1}.N;
	NMSPT = models_MSPT{1}.N;
	csvMSP(idx, 6) = models_MSPT{1}.dev{idx}*NMSP/NMSPT;
	%Compute chi^2 value
	csvMSP(idx, 7) = csvMSP(idx,2)-csvMSP(idx,6);
	%p-value
	csvMSP(idx, 8) = 1-chi2cdf(csvMSP(idx, 7), csvMSP(idx, 5) - csvMSP(idx, 3));
	%Report change in AIC and change in BIC for each unit when different covariates are added
	%Change in AIC is the twice the change in number of parameters, minus twice the change in the maximized likelihood (deviance)
	%csvMSP(idx, 9) = 2*csvMSP(idx, 5+L+j) - 2*csvMSP(idx, 5) + csvMSP(idx, 5+j) - csvMSP(idx, 4);
	%Change in BIC
	%csvMSP(idx, 10) = (csvMSP(idx, 5+L+j) - csvMSP(idx, 5))*log(NMSV) + csvMSP(idx, 5+j) - csvMSP(idx, 4);
end
%Save all data as a csv for analysis in excel or similar
csvwrite_heading('./worksheets/12_1_2014/AOD_MSPT.csv', csvMSP, MSP_headings);
