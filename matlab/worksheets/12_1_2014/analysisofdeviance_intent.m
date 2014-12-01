%Analysis of deviance scripts

%Preprocess data (timebis of 2ms, smooth trajectories)
%Apply no offset
nevfile = './testdata/20130117SpankyUtah001.nev';
binsize = 0.002;
dt_sp = binsize;
dt_pos = 0.2;
dt_tar = 0.2;
dt_vel = 0.1;
offset = 0.0;
threshold = 5;
processed = preprocess_spline(nevfile, binsize, threshold, offset);
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
	models_MST{idx} = MLE_glmfit(data, const);
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

MSP_headings = {'Name', 'Dev M', 'nM', 'Dev MS', 'nMS', 'Dev MSP', blank{:}, 'nMSP', blank{:}, 'chi2', blank{:}, 'p-val', blank{:}, 'AIC', blank{:}, 'BIC'};
for idx = 1:nU
	%Unit name
	csvMSP(idx, 1) = str2num(processed.unitnames{idx});
	%Dev M
	csvMSP(idx, 2) = model_M.dev{idx};
	%nM
	csvMSP(idx, 3) = size(model_M.b_hat,2);
	%Dev MS
	csvMSP(idx, 4) = model_MS.dev{idx};
	%nMS
	csvMSP(idx, 5) = size(model_MS.b_hat, 2);
	%Number of data points in MS model
	NMS = model_MS.N;
	for j = 1:L
		%Adj Dev MSP 6
		%Adjusted since the number of data points in the MSP fits is different to the number of data points in the MS fits, which
		%Fewer data points means fewer deviances... which can make MSP fits better than they seem
		NMSP = models_MSP{j}.N;
		csvMSP(idx, 5+j) = models_MSP{j}.dev{idx}*NMS/NMSP;
		%nMSP 7
		csvMSP(idx, 5+L+j) = size(models_MSP{j}.b_hat, 2);
		%Compute chi^2 value 8
		csvMSP(idx, 5+2*L+j) = csvMSP(idx,4)-csvMSP(idx,5+j);
		%p-value 9
		csvMSP(idx, 5+3*L+j) = 1-chi2cdf(csvMSP(idx, 5+2*L+j), csvMSP(idx, 5+L+j) - csvMSP(idx, 5));
		%Report change in AIC and change in BIC for each unit when different covariates are added
		%Change in AIC is the twice the change in number of parameters, minus twice the change in the maximized likelihood (deviance)
		%10
		csvMSP(idx, 5+4*L+j) = 2*csvMSP(idx, 5+L+j) - 2*csvMSP(idx, 5) + csvMSP(idx, 5+j) - csvMSP(idx, 4);
		%Change in BIC
		%11
		csvMSP(idx, 5+5*L+j) = (csvMSP(idx, 5+L+j) - csvMSP(idx, 5))*log(NMSV) + csvMSP(idx, 5+j) - csvMSP(idx, 4);
	end
end
%Save all data as a csv for analysis in excel or similar
csvwrite_heading('./worksheets/12_1_2014/AOD_MSPT.csv', csvMSP, MSP_headings);
