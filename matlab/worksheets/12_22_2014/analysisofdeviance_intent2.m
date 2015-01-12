%Analysis of deviance scripts for intentional signals

%Preprocess data (timebis of 2ms, smooth trajectories)
%Apply no offset
nevfile = './testdata/20130117SpankyUtah001.nev';
labviewfile = './testdata/Spanky_2013-01-17-1325.mat';
fn_out = './worksheets/12_22_2014/plots/test';
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
	fn_out = ['./worksheets/12_22_2014/plots/AOD_MSP_nK_' num2str(nK_pos)];
	plot_filters(models_MSP{idx}, data, processed, fn_out);
	plot_predictions(models_MSP{idx}, data, processed, fn_out);
end

%- MSR, spike history, relative position
const = 'on';
nK_sp = 100; 
nK_poss = [5];
models_MSR = {};
for idx = 1:length(nK_poss)
	idx
	nK_pos = nK_poss(idx);
	data = filters_sp_relpos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
	%Fit each of the above GLMs
	models_MSR{idx} = MLE_glmfit(data, const);
	%Make plots of each filter fitted, predictions of each unit, and record the deviance
	fn_out = ['./worksheets/12_22_2014/plots/AOD_MSR_nK_' num2str(nK_pos)];
	plot_filters(models_MSR{idx}, data, processed, fn_out);
	plot_predictions(models_MSR{idx}, data, processed, fn_out);
end

%Plot log likelihood of one vs another
clf
fn_out = './worksheets/12_22_2014/plots/logliRelVsAbsPos.eps';
xx = linspace(min([models_MSR{1}.logli, models_MSP{1}.logli]), max([models_MSR{1}.logli, models_MSP{1}.logli]), 100);
plot(models_MSR{1}.logli, models_MSP{1}.logli, '.', xx, xx);
xlabel('Rel pos');
ylabel('Abs pos');
title('Log likelihood');
saveplot(gcf, fn_out);

%save fitted models for later use
save('./worksheets/12_22_2014/AOD_fittedmodels.mat', 'models_MSP', 'models_MSR')
load('./worksheets/12_22_2014/AOD_fittedmodels.mat', 'models_MSP', 'models_MSR')

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Do the same with velocity%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%MSV mean, spike history and velocity
%Do for a range of filter sizes
const = 'on';
nK_sp = 100; 
nK_poss = [5];
models_MSV = {};
for idx = 1:length(nK_poss)
	idx
	nK_pos = nK_poss(idx);
	data = filters_sp_vel(processed, nK_sp, nK_pos, dt_sp, dt_vel);
	%Fit each of the above GLMs
	models_MSV{idx} = MLE_glmfit(data, const);
	%Make plots of each filter fitted, predictions of each unit, and record the deviance
	fn_out = ['./worksheets/12_22_2014/plots/AOD_MSV_nK_' num2str(nK_pos)];
	plot_filters(models_MSV{idx}, data, processed, fn_out);
	plot_predictions(models_MSV{idx}, data, processed, fn_out);
end

%MSVP
%MSV mean, spike history and velocity
%Do for a range of filter sizes
const = 'on';
nK_sp = 100; 
nK_poss = [5];
nK_vel = 5;
models_MSVP = {};
for idx = 1:length(nK_poss)
	idx
	nK_pos = nK_poss(idx);
	data = filters_sp_pos_vel(processed, nK_sp, nK_pos, nK_vel, dt_sp, dt_pos, dt_vel);
	%Fit each of the above GLMs
	models_MSVP{idx} = MLE_glmfit(data, const);
	%Make plots of each filter fitted, predictions of each unit, and record the deviance
	fn_out = ['./worksheets/12_22_2014/plots/AOD_MSVP_nK_' num2str(nK_pos)];
	plot_filters(models_MSVP{idx}, data, processed, fn_out);
	plot_predictions(models_MSVP{idx}, data, processed, fn_out);
end

%MSVPr


%MSPrT


%Plot log likelihood of one vs another
clf
fn_out = './worksheets/12_22_2014/plots/logliRelVsAbsVel.eps';
xx = linspace(min([models_MSVr{1}.logli, models_MSV{1}.logli]), max([models_MSVr{1}.logli, models_MSV{1}.logli]), 100);
plot(models_MSVr{1}.logli, models_MSV{1}.logli, '.', xx, xx);
xlabel('Rel vel');
ylabel('Abs vel');
title('Log likelihood');
saveplot(gcf, fn_out);

%Plot sum of likelihoods across all units, and compare for each of the above models (use the AIC and the BIC, here)
