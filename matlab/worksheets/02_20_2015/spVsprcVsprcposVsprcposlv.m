const = 'on';
pre = load('./testdata/test_preprocess_spline_lv_55_1.mat');
dt_sp = 0.002;
dt_pos = 0.2;

%sp
nK_sp = 100; 
nK_pos = 0;
data = filters_sp_pos(pre.processed, nK_sp, nK_pos, dt_sp, dt_pos);
model_sp = MLE_glmfit(data, const);

%sprc
nK_sp = 100; 
nK_pos = 0;
%Load test preprocessed data
data = filters_sprc_pos(pre.processed, nK_sp, nK_pos, dt_sp, dt_pos);
model_sprc = MLE_glmfit(data, const);

%sprc_pos
nK_sp = 100; 
nK_pos = 6;
%Load test preprocessed data
data = filters_sprc_pos(pre.processed, nK_sp, nK_pos, dt_sp, dt_pos);
model_sprcpos = MLE_glmfit(data, const);

%sprc_pos_lv
nK_sp = 100; 
nK_pos = 6;
%Load test preprocessed data
data = filters_sp_pos_lv(pre.processed, nK_sp, nK_pos, dt_sp, dt_pos);
model_sprcposlv = MLE_glmfit(data, const);


%Compare sp and sprc deviances:
