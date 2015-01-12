%Load test preprocessed data
fn_out = './worksheets/01_12_2015/plots/testdecode.eps';
pre = load('./testdata/test_preprocess_spline_short24.mat');
%pre = load('./testdata/test_preprocess_spline.mat');
const = 'on';
nK_sp = 50; 
nK_pos = 1;
dt_sp = 0.002;
dt_pos = 0.2;
R = 1;
data = filters_sp_pos(pre.processed, nK_sp, nK_pos, dt_sp, dt_pos);
model = MLE_glmfit(data, const);
[F, Q, mu] = fit_AR_LS(data.torque, R);
decoded_torque = glm_decode(pre.processed, data, model, F, Q, mu, fn_out);

%Switch RU and FE axes to see if correlation (performance) changes accordingly. This will tell us
%that the large difference in RU and FE performance is due to the data, and not how it's processed,
%which would be indicative of a bug somewhere.
%This difference, if real, is interesting, as it may be because of the fact that many more trials in 
%which Spanky was performing manual or brain RU control were performed, compared with FE control 
%(horizontal, compared with vertical). The neurons may represent more horizontal position than vertical, 
%hence the improved ability to decode horizontal position

%Switch and repeat
processed_switch = pre.processed;
processed_switch.torque = fliplr(pre.processed.torque);
processed_switch.dtorque = fliplr(pre.processed.dtorque);
processed_switch.ddtorque = fliplr(pre.processed.ddtorque);
data = filters_sp_pos(processed_switch, nK_sp, nK_pos, dt_sp, dt_pos);
model = MLE_glmfit(data, const);
[F, Q, mu] = fit_AR_LS(data.torque, R);
%load('./testdata/testglm.mat')
decoded_torque = glm_decode(processed_switch, data, model, F, Q, mu, fn_out);
