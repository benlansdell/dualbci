%Analysis of deviance scripts

%Preprocess data (timebis of 2ms, smooth trajectories)
%Apply no offset
nevfile = './testdata/20130117SpankyUtah001.nev';
binsize = 0.002;
dt_sp = binsize;
dt_pos = 0.2;
dt_vel = 0.1;
offset = 0.0;
threshold = 5;
processed = preprocess_spline(nevfile, binsize, threshold, offset);
nU = length(processed.unitnames);
N = size(processed.rates, 1);

%- M mean only (no filters at all)
const = 'on';
%Prepare data for fitting GLM
%Prepare:
nK_sp = 0; 
nK_pos = 0;
data = filters_sp_pos(processed, nK_sp, nK_pos);
%Fit each of the above GLMs
model_M = MLE_glmfit(data, const);
%Make plots of each filter fitted, predictions of each unit, and record the deviance
fn_out = './worksheets/11_11_2014/plots/AOD_M';
plot_predictions(model_M, data, processed, fn_out);

%- MS mean plus spike history
const = 'on';
nK_sp = 100; 
nK_pos = 0;
data = filters_sp_pos(processed, nK_sp, nK_pos);
%Fit each of the above GLMs
model_MS = MLE_glmfit(data, const);
%Make plots of each filter fitted, predictions of each unit, and record the deviance
fn_out = './worksheets/11_11_2014/plots/AOD_MS';
plot_filters(model_MS, data, processed, fn_out);
plot_predictions(model_MS, data, processed, fn_out);


%- MSP mean, spike history and position
%Do for a range of filter sizes
const = 'on';
nK_sp = 100; 
nK_poss = [1 2 3 4 5 6 7 8 9 10];
models_MSP = {};
for idx = 1:length(nK_poss)
	idx
	nK_pos = nK_poss(idx);
	data = filters_sp_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
	%Fit each of the above GLMs
	models_MSP{idx} = MLE_glmfit(data, const);
	%Make plots of each filter fitted, predictions of each unit, and record the deviance
	fn_out = ['./worksheets/11_11_2014/plots/AOD_MSP_nK_' num2str(nK_pos)];
	plot_filters(models_MSP{idx}, data, processed, fn_out);
	plot_predictions(models_MSP{idx}, data, processed, fn_out);
end

%- MSV mean, spike history and velocity
const = 'on';
nK_sp = 100; 
nK_vels = [1 2 3 4 5 6 7 8 9 10];
models_MSV = {};
for idx = 1:length(nK_vels)
	nK_vel = nK_vels(idx);
	data = filters_sp_vel(processed, nK_sp, nK_vel, dt_sp, dt_vel);
	%Fit each of the above GLMs
	models_MSV{idx} = MLE_glmfit(data, const);
	%Make plots of each filter fitted, predictions of each unit, and record the deviance
	fn_out = ['./worksheets/11_11_2014/plots/AOD_MSV_nK_' num2str(nK_vel)];
	plot_filters(models_MSV{idx}, data, processed, fn_out);
	plot_predictions(models_MSV{idx}, data, processed, fn_out);
end

%save fitted models for later use
%save('./worksheets/11_11_2014/AOD_fittedmodels.mat', 'model_M', 'model_MS', 'models_MSP')
save('./worksheets/11_11_2014/AOD_fittedmodels_vel.mat', 'model_M', 'model_MS', 'models_MSV')

L = length(nK_poss);
L = length(nK_vels);
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
		csvMSP(idx, 5+5*L+j) = (csvMSP(idx, 5+L+j) - csvMSP(idx, 5))*log(N) + csvMSP(idx, 5+j) - csvMSP(idx, 4);
	end
end
%Save all data as a csv for analysis in excel or similar
csvwrite_heading('./worksheets/11_11_2014/AOD_MSP.csv', csvMSP, MSP_headings);

csvMSV = zeros(nU, 11);
MSV_headings = {'Name', 'Dev M', 'nM', 'Dev MS', 'nMS', 'Dev MSV', blank{:}, 'nMSV', blank{:}, 'chi2', blank{:}, 'p-val', blank{:}, 'AIC', blank{:}, 'BIC'};
for idx = 1:nU
	%Unit name
	csvMSV(idx, 1) = str2num(processed.unitnames{idx});
	%Dev M
	csvMSV(idx, 2) = model_M.dev{idx};
	%nM
	csvMSV(idx, 3) = size(model_M.b_hat,2);
	%Dev MS
	csvMSV(idx, 4) = model_MS.dev{idx};
	%nMS
	csvMSV(idx, 5) = size(model_MS.b_hat, 2);
	%Number of data points in MS model
	NMS = model_MS.N;
	for j = 1:L
		%Adj Dev MSP 6
		%Adjusted since the number of data points in the MSV fits is different to the number of data points in the MS fits, which
		%Fewer data points means fewer deviances... which can make MSV fits better than they seem
		NMSV = models_MSV{j}.N;
		csvMSV(idx, 5+j) = models_MSV{j}.dev{idx}*NMS/NMSV;
		%nMSP 7
		csvMSV(idx, 5+L+j) = size(models_MSV{j}.b_hat, 2);
		%Compute chi^2 value 8
		csvMSV(idx, 5+2*L+j) = csvMSV(idx,4)-csvMSV(idx,5+j);
		%p-value 9
		csvMSV(idx, 5+3*L+j) = 1-chi2cdf(csvMSV(idx, 5+2*L+j), csvMSV(idx, 5+L+j) - csvMSV(idx, 5));
		%Report change in AIC and change in BIC for each unit when different covariates are added
		%Change in AIC is the twice the change in number of parameters, minus twice the change in the maximized likelihood (deviance)
		%10
		csvMSV(idx, 5+4*L+j) = 2*csvMSV(idx, 5+L+j) - 2*csvMSV(idx, 5) + csvMSV(idx, 5+j) - csvMSV(idx, 4);
		%Change in BIC
		%11
		csvMSV(idx, 5+5*L+j) = (csvMSV(idx, 5+L+j) - csvMSV(idx, 5))*log(N) + csvMSV(idx, 5+j) - csvMSV(idx, 4);
	end
end
%Save all data as a csv for analysis in excel or similar
csvwrite_heading('./worksheets/11_11_2014/AOD_MSV.csv', csvMSV, MSV_headings);


%Populate model structures with number of data points used for model, N
for j = 1:L
	j
	nK_vel = nK_vels(j);
	data = filters_sp_vel(processed, nK_sp, nK_vel, dt_sp, dt_vel);
	models_MSV{j}.N = size(data.y, 2);
end

%Clear current structures of model stats like residuals that take up memory
%for idx = 1:nU
%	for j = 1:L
%		models_MSP{j}.stats{idx} = rmfield(models_MSP{j}.stats{idx}, {'resid', 'residp', 'residd', 'resida', 'wts'});
%	end
%end

%for idx = 1:nU
%	model_M.stats{idx} = rmfield(model_M.stats{idx}, {'resid', 'residp', 'residd', 'resida', 'wts'});
%end