%%%%%%%%%%%%%%%%%%
%test_optimizer.m%
%%%%%%%%%%%%%%%%%%

%Code to see what the difference between the glmfit IRLS, and the fminuc MLE methods for different types of filters


%%%%%%%%%%%%
%First part%
%%%%%%%%%%%%

%Compare very simple filters
fn_out = './worksheets/10_11_2014/plots/compare_simple.eps';
const = 'on';
nK_sp = 0; 
nK_pos = 1;
%Load test preprocessed data
pre = load('./testdata/test_preprocess_spline_short.mat');
data = filters_sp_pos(pre.processed, nK_sp, nK_pos);
[model_SD, intermediates] = MLE_SD(data, const);
model_IRLS = MLE_glmfit(data, const);

%Compare the values
subplot(2,1,1)
plot(model_SD.b_hat)
title('fminunc filters')
subplot(2,1,2)
plot(model_IRLS.b_hat)
title('IRLS filters')
saveplot(gcf, fn_out);

%Compare the deviances
dev_SD = deviance(model_SD, data)
dev_IRLS = deviance(model_IRLS, data)

%%%%%%%%%%%%%
%Second part%
%%%%%%%%%%%%%

%Compare intermediate filters
fn_out = './worksheets/10_11_2014/plots/compare_intermediate.eps';
const = 'on';
nK_sp = 10; 
nK_pos = 5;
dt_sp = 0.002;
dt_pos = 0.05;
%Load test preprocessed data
data = filters_sp_pos(pre.processed, nK_sp, nK_pos, dt_sp, dt_pos);
[model_SD, intermediates] = MLE_SD(data, const);
model_IRLS = MLE_glmfit(data, const);

%Compare the values
subplot(2,1,1)
plot(model_SD.b_hat)
title('fminunc filters. Smoothed torque data')
subplot(2,1,2)
plot(model_IRLS.b_hat)
title('IRLS filters. Smoothed torque data')
saveplot(gcf, fn_out);

%Compare the deviances
dev_SD = deviance(model_SD, data)
dev_IRLS = deviance(model_IRLS, data)

%%%%%%%%%%%%
%Third part%
%%%%%%%%%%%%

%Compare intermediate filters with noisier input data
fn_out = './worksheets/10_11_2014/plots/compare_intermediate_unsmoothed.eps';
const = 'on';
nK_sp = 10; 
nK_pos = 5;
dt_sp = 0.002;
dt_pos = 0.05;
%Load test preprocessed data
pre = load('./testdata/test_preprocess_short.mat');
data = filters_sp_pos(pre.processed, nK_sp, nK_pos, dt_sp, dt_pos);
[model_SD, intermediates] = MLE_SD(data, const);
model_IRLS = MLE_glmfit(data, const);

%Compare the values
subplot(2,1,1)
plot(model_SD.b_hat)
title('fminunc filters. Unsmoothed torque data')
subplot(2,1,2)
plot(model_IRLS.b_hat)
title('IRLS filters. Unsmoothed torque data')
saveplot(gcf, fn_out);

%Compare the deviances
dev_SD = deviance(model_SD, data)
dev_IRLS = deviance(model_IRLS, data)

%The plots here look a lot more similar... the smoothing doesn't help the IRLS method, it would seem!

%%%%%%%%%%%%%
%Fourth part%
%%%%%%%%%%%%%

%Compare larger filters
fn_out = './worksheets/10_11_2014/plots/compare_larger_unsmoothed.eps';
const = 'on';
nK_sp = 100; 
nK_pos = 100;
%Load test preprocessed data
pre = load('./testdata/test_preprocess_short.mat');
data = filters_sp_pos(pre.processed, nK_sp, nK_pos);
[model_SD, intermediates] = MLE_SD(data, const);
model_IRLS = MLE_glmfit(data, const);

%Compare the values
subplot(2,1,1)
plot(model_SD.b_hat)
title('fminunc filters. Unsmoothed torque data')
subplot(2,1,2)
plot(model_IRLS.b_hat)
title('IRLS filters. Unsmoothed torque data')
saveplot(gcf, fn_out);

%Compare the deviances
dev_SD = deviance(model_SD, data)
dev_IRLS = deviance(model_IRLS, data)

%%%%%%%%%%%%
%Fifth part%
%%%%%%%%%%%%

%Compare larger filters, smoothed
fn_out = './worksheets/10_11_2014/plots/compare_larger_smoothed.eps';
const = 'on';
nK_sp = 100; 
nK_pos = 100;
%Load test preprocessed data
pre = load('./testdata/test_preprocess_spline_short.mat');
data = filters_sp_pos(pre.processed, nK_sp, nK_pos);
[model_SD, intermediates] = MLE_SD(data, const);
[model_IRLS, intermediates] = MLE_glmfit(data, const);

%Compare the values
subplot(2,1,1)
plot(model_SD.b_hat)
title('fminunc filters. Smoothed torque data')
subplot(2,1,2)
plot(model_IRLS.b_hat)
title('IRLS filters. Smoothed torque data')
saveplot(gcf, fn_out);

%Compare the deviances
dev_SD = deviance(model_SD, data)
dev_IRLS = deviance(model_IRLS, data)



%%%%%%%%%%%%
%Sixth part%
%%%%%%%%%%%%

%Play with different basis functions...

%Save the outcome of all these trials...
save('./worksheets/10_11_2014/data.mat');