%Comparing linear and GLM models

%%%%%%%%%%%%%%%%%%%%%%%%%%
%At 60Hz timebins network%
%%%%%%%%%%%%%%%%%%%%%%%%%%

const = 'on';
nK_sp = 6; 
nK_pos = 6;
%Preprocess spline at 60Hz
pre = load('./testdata/test_preprocess_spline_60hz_short24.mat');
processed = pre.processed;
%Load data for models
data = filters_sp_pos_network(processed, nK_sp, nK_pos);
%Fit a GLM to the network
modelGLM = MLE_glmfit_network(data, const);
%Fit a LM
modelLM = MLE_lmfit_network(data, const);
nU = length(processed.unitnames);
%Plot the deviance of each unit as LM vs GLM
llGLM = modelGLM.logli;
llLM = modelLM.logli;
%Save plot
%for i = 1:nU
%	llGLM(i) = modelGLM.logli{i};
%	llLM(i) = modelLM.logli{i};
%end
fn_out = './worksheets/12_9_2014/plots/llLMvsGLM_bs_16.eps';
mx = 10000;
plot(llGLM, llLM, '.', -mx:0, -mx:0)
xlabel('Log-likelihood GLM')
ylabel('Log-likelihood LM')
title('dt = 0.016')
xlim([-mx 0])
ylim([-mx 0])
saveplot(gcf, fn_out);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%At 500Hz timebins, individual units%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

const = 'on';
nK_sp = 100; 
nK_pos = 6;
dt_sp = 0.002;
dt_pos = 0.05;
%Preprocess spline at 500Hz
pre = load('./testdata/test_preprocess_spline.mat');
processed = pre.processed;
%Load data for models
data = filters_sp_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
%Fit a GLM to the network
modelGLM = MLE_glmfit(data, const);
%Fit a LM
modelLM = MLE_lmfit(data, const);
nU = length(processed.unitnames);
%Plot the deviance of each unit as LM vs GLM
llGLM = modelGLM.logli;
llLM = modelLM.logli;
fn_out = './worksheets/12_9_2014/plots/llLMvsGLM_bs_2.eps';
mx = 20000;
plot(llGLM, llLM, '.', -mx:0, -mx:0)
xlabel('Log-likelihood GLM')
ylabel('Log-likelihood LM')
title('dt = 0.002')
xlim([-mx 0])
ylim([-mx 0])
saveplot(gcf, fn_out);