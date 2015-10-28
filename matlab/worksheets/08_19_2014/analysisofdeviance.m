%Analysis of deviance scripts

%Preprocess data (timebis of 2ms, smooth trajectories)
%Apply no offset
nevfile = './testdata/20130117SpankyUtah001.nev';
binsize = 0.002;
offset = 0.0;
threshold = 5;
processed = preprocess_spline(nevfile, binsize, threshold, offset);
nU = length(processed.unitnames);
N = size(processed.rates, 1);

%Prepare data for fitting GLM
%Prepare:
%- M mean only (no filters at all)
const = 'on';
nK_sp = 0; 
nK_pos = 0;
data = filters_sp_pos(processed, nK_sp, nK_pos);
%Fit each of the above GLMs
model_M = MLE_glmfit(data, const);
%Make plots of each filter fitted, predictions of each unit, and record the deviance
fn_out = './worksheets/08_19_2014/plots/AOD_M';
plot_predictions(model_M, data, processed, fn_out);

%- MS mean plus spike history
const = 'on';
nK_sp = 100; 
nK_pos = 0;
data = filters_sp_pos(processed, nK_sp, nK_pos);
%Fit each of the above GLMs
model_MS = MLE_glmfit(data, const);
%Make plots of each filter fitted, predictions of each unit, and record the deviance
fn_out = './worksheets/08_19_2014/plots/AOD_MS';
plot_filters(model_MS, data, processed, fn_out);
plot_predictions(model_MS, data, processed, fn_out);

%- MSP mean, spike history and position
const = 'on';
nK_sp = 100; 
nK_pos = 1;
data = filters_sp_pos(processed, nK_sp, nK_pos);
%Fit each of the above GLMs
model_MSP = MLE_glmfit(data, const);
%Make plots of each filter fitted, predictions of each unit, and record the deviance
fn_out = './worksheets/08_19_2014/plots/AOD_MSP';
plot_filters(model_MSP, data, processed, fn_out);
plot_predictions(model_MSP, data, processed, fn_out);

%- MSV mean, spike history and velocity
const = 'on';
nK_sp = 100; 
nK_vel = 1;
data = filters_sp_vel(processed, nK_sp, nK_vel);
%Fit each of the above GLMs
model_MSV = MLE_glmfit(data, const);
%Make plots of each filter fitted, predictions of each unit, and record the deviance
fn_out = './worksheets/08_19_2014/plots/plot_filters';
plot_filters(model_MSV, data, processed, fn_out);
plot_predictions(model_MSV, data, processed, fn_out);

%- MSA mean, spike history and acceleration
const = 'on';
nK_sp = 100; 
nK_acc = 1;
data = filters_sp_accel(processed, nK_sp, nK_acc);
%Fit each of the above GLMs
model_MSA = MLE_glmfit(data, const);
%Make plots of each filter fitted, predictions of each unit, and record the deviance
fn_out = './worksheets/08_19_2014/plots/AOD_MSA';
plot_filters(model_MSA, data, processed, fn_out);
plot_predictions(model_MSA, data, processed, fn_out);

%- MSD mean, spike history and direction and speed
const = 'on';
nK_sp = 100; 
nK_acc = 1;
data = filters_sp_dir(processed, nK_sp, nK_acc);
%Fit each of the above GLMs
model_MSD = MLE_glmfit(data, const);
%Make plots of each filter fitted, predictions of each unit, and record the deviance
fn_out = './worksheets/08_19_2014/plots/AOD_MSD';
plot_filters(model_MSD, data, processed, fn_out);
plot_predictions(model_MSD, data, processed, fn_out);

%- MSP mean, spike history and position with larger filter
const = 'on';
nK_sp = 100; 
nK_pos = 5;
dt_sp = 0.002; dt_pos = 0.05;
data = filters_sp_pos(processed, nK_sp, nK_pos, dt_pos, dt_pos);
%Fit each of the above GLMs
model_MSP_nK_pos_5 = MLE_glmfit(data, const);
%Make plots of each filter fitted, predictions of each unit, and record the deviance
fn_out = './worksheets/08_19_2014/plots/AOD_MSP_nK_pos_5';
plot_filters(model_MSP, data, processed, fn_out);
plot_predictions(model_MSP, data, processed, fn_out);


%save fitted models
save('./worksheets/08_19_2014/AOD_fittedmodels.mat', 'model_M', 'model_MS', 'model_MSP', 'model_MSV', 'model_MSA')

csvMSP = zeros(nU, 11);
%Fields to save (nested model MS vs MSP)
%Delta AIC (MSP), Delta BIC (MSP), correlation
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
	%Dev MSP
	csvMSP(idx, 6) = model_MSP.dev{idx};
	%nMSP
	csvMSP(idx, 7) = size(model_MSP.b_hat, 2);
	%Compute chi^2 value
	csvMSP(idx, 8) = csvMSP(idx,4)-csvMSP(idx,6);
	%p-value
	csvMSP(idx, 9) = 1-chi2cdf(csvMSP(idx, 8), csvMSP(idx, 7) - csvMSP(idx, 5));
	%Report change in AIC and change in BIC for each unit when different covariates are added
	%Change in AIC is the twice the change in number of parameters, minus twice the change in the maximized likelihood (deviance)
	csvMSP(idx, 10) = 2*csvMSP(idx, 7) - 2*csvMSP(idx, 5) + csvMSP(idx, 6) - csvMSP(idx, 4);
	%Change in BIC
	csvMSP(idx, 11) = (csvMSP(idx, 7) - csvMSP(idx, 5))*log(N) + csvMSP(idx, 6) - csvMSP(idx, 4);
end

%Save all data as a csv for analysis in excel or similar
csvwrite('./worksheets/08_19_2014/AOD_MSP.csv', csvMSP);

csvMSV = zeros(nU, 11);
%Fields to save (nested model MS vs MSV)
%Delta AIC (MSV), Delta BIC (MSV), correlation
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
	%Dev MSV
	csvMSV(idx, 6) = model_MSV.dev{idx};
	%nMSV
	csvMSV(idx, 7) = size(model_MSV.b_hat, 2);
	%Compute chi^2 value
	csvMSV(idx, 8) = csvMSV(idx,4)-csvMSV(idx,6);
	%p-value
	csvMSV(idx, 9) = 1-chi2cdf(csvMSV(idx, 8), csvMSV(idx, 7) - csvMSV(idx, 5));
	%Report change in AIC and change in BIC for each unit when different covariates are added
	%Change in AIC is the twice the change in number of parameters, minus twice the change in the maximized likelihood (deviance)
	csvMSV(idx, 10) = 2*csvMSV(idx, 7) - 2*csvMSV(idx, 5) + csvMSV(idx, 6) - csvMSV(idx, 4);
	%Change in BIC
	csvMSV(idx, 11) = (csvMSV(idx, 7) - csvMSV(idx, 5))*log(N) + csvMSV(idx, 6) - csvMSV(idx, 4);
end

%Save all data as a csv for analysis in excel or similar
csvwrite('./worksheets/08_19_2014/AOD_MSV.csv', csvMSV);

csvMSA = zeros(nU, 11);
%Fields to save (nested model MS vs MSA)
%Delta AIC (MSA), Delta BIC (MSA), correlation
for idx = 1:nU
	%Unit name
	csvMSA(idx, 1) = str2num(processed.unitnames{idx});
	%Dev M
	csvMSA(idx, 2) = model_M.dev{idx};
	%nM
	csvMSA(idx, 3) = size(model_M.b_hat,2);
	%Dev MS
	csvMSA(idx, 4) = model_MS.dev{idx};
	%nMS
	csvMSA(idx, 5) = size(model_MS.b_hat, 2);
	%Dev MSA
	csvMSA(idx, 6) = model_MSA.dev{idx};
	%nMSA
	csvMSA(idx, 7) = size(model_MSA.b_hat, 2);
	%Compute chi^2 value
	csvMSA(idx, 8) = csvMSA(idx,4)-csvMSA(idx,6);
	%p-value
	csvMSA(idx, 9) = 1-chi2cdf(csvMSA(idx, 8), csvMSA(idx, 7) - csvMSA(idx, 5));
	%Report change in AIC and change in BIC for each unit when different covariates are added
	%Change in AIC is the twice the change in number of parameters, minus twice the change in the maximized likelihood (deviance)
	csvMSA(idx, 10) = 2*csvMSA(idx, 7) - 2*csvMSA(idx, 5) + csvMSA(idx, 6) - csvMSA(idx, 4);
	%Change in BIC
	csvMSA(idx, 11) = (csvMSA(idx, 7) - csvMSA(idx, 5))*log(N) + csvMSA(idx, 6) - csvMSA(idx, 4);
end

%Save all data as a csv for analysis in excel or similar
csvwrite('./worksheets/08_19_2014/AOD_MSA.csv', csvMSA);