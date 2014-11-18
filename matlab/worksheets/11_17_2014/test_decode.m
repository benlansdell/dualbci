%Load test preprocessed data
fn_out = './worksheets/11_17_2014/plots/testdecode.eps';
pre = load('./testdata/test_preprocess_spline_short24.mat');
%pre = load('./testdata/test_preprocess_spline.mat');
const = 'on';
nK_sp = 50; 
nK_pos = 4;
dt_sp = 0.002;
dt_pos = 0.2;
R = 5;
P = nK_pos;
d = dt_pos/dt_sp;
data = filters_sp_pos(pre.processed, nK_sp, nK_pos, dt_sp, dt_pos);
model = MLE_glmfit(data, const);
[F, Q, mu] = fit_AR_LS_multi(data.torque, R, P, d);
%load('./testdata/testglm.mat')
decoded_torque = glm_decode_multi(pre.processed, data, model, F, Q, mu, R, P, d, fn_out);
