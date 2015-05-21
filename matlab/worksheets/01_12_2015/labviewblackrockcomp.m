nevfile = './testdata/20130117SpankyUtah001.nev';
matfile = './testdata/Spanky_2013-01-17-1325.mat';
binsize = 0.002;
offset = 0.0;
threshold = 5;
fn_out = './worksheets/01_12_2015/test_spline_pre_lv.eps';
processed = preprocess_spline_lv(nevfile, matfile, binsize, threshold, offset, fn_out);

%Check that deviance and filters fit for both torque and labview datasets are similar...
nK_sp = 100;
nK_pos = 5;
nK_vel = 5;
dt_sp = 0.002;
dt_pos = 0.2;
dt_vel = 0.1;
const = 'on';

%Position model
%Torque data:
data = filters_sp_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
model = MLE_glmfit(data, const);
fn_out = './worksheets/01_12_2015/plots/torque';
plot_filters(model, data, processed, fn_out);
clear data; %Save memory

%Lab view data:
datalv = filters_sp_pos_lv(processed, nK_sp, nK_pos, dt_sp, dt_pos);
modellv = MLE_glmfit(datalv, const);
fn_out = './worksheets/01_12_2015/plots/cursor';
plot_filters(modellv, datalv, processed, fn_out);
clear datalv;


%%%%Vel

%Torque data:
data = filters_sp_vel(processed, nK_sp, nK_vel, dt_sp, dt_vel);
model = MLE_glmfit(data, const);
fn_out = './worksheets/01_12_2015/plots/torque_vel';
plot_filters(model, data, processed, fn_out);
clear data; %Save memory

%Lab view data:
datalv = filters_sp_vel_lv(processed, nK_sp, nK_vel, dt_sp, dt_vel);
modellv = MLE_glmfit(datalv, const);
fn_out = './worksheets/01_12_2015/plots/cursor_vel';
plot_filters(modellv, datalv, processed, fn_out);
clear datalv;

%Compare deviances:
modelcomp = ...
    [3.2719e+04-3.2723e+04;
    2.6631e+04-2.6621e+04;
    3.7305e+04-3.7322e+04;
    5.2619e+04-5.2619e+04;
    8.6347e+04-8.6348e+04;
    6.1947e+04-6.2009e+04;
    6.8468e+04-6.8506e+04;
    5.2079e+04-5.2077e+04;
    5.9829e+04-5.9841e+04;
    6.3084e+04-6.3086e+04;
    7.3217e+04-7.3258e+04;
    5.9806e+04-5.9893e+04;
    4.3210e+04-4.3221e+04;
    2.2259e+04-2.2282e+04;
    6.8030e+04-6.8061e+04;
    2.9697e+04-2.9709e+04;
    2.8945e+04-2.8956e+04;
    6.0222e+04-6.0224e+04;
    9.9663e+04-9.9663e+04;
    6.0407e+04-6.0410e+04;
    7.1175e+04-7.1161e+04;
    5.7007e+04-5.7005e+04;
    4.3475e+04-4.3548e+04;
    3.6168e+04-3.6165e+04;]
%Pretty similar, in general... but not identical

%Also a file that uses a BCI:

%1D horizontal position (I believe)
nevfile = './testdata/20130117SpankyUtah005.nev';
matfile = './testdata/Spanky_2013-01-17-1325.mat';
binsize = 0.002;
offset = 0.0;
threshold = 5;
fn_out = './worksheets/01_12_2015/test_spline_pre_BCI_lv.eps';
processed = preprocess_spline_lv(nevfile, matfile, binsize, threshold, offset, fn_out);


