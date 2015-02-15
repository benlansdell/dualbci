const = 'on';
fn_out = './worksheets/02_05_2015/20130117SpankyUtah001';
nK_sp = 100; 
nK_pos = 5;
nK_vel = 5;
nK_tar = 5;
nK_acc = 5;
dt_pos = 0.2;
dt_vel = 0.2;
dt_acc = 0.2;
dt_tar = 0.2;
dt_sp = 0.002;
%Load test preprocessed data
pre = load('./testdata/test_preprocess_spline_55_1.mat');
%Create raised cosine spike history filters
data = filters_sprc_pos(pre.processed, nK_sp, nK_pos);
data = filters_sprc_vel(pre.processed, nK_sp, nK_vel);
data = filters_sprc_revvel(pre.processed, nK_sp, nK_vel);
data = filters_sprc_accel(pre.processed, nK_sp, nK_acc);
data = filters_sprc_dir(pre.processed, nK_sp, nK_pos);
data = filters_sprc_revpos(pre.processed, nK_sp, nK_pos);
data = filters_sprc_pos_vel(pre.processed, nK_sp, nK_pos, nK_vel);
data = filters_sprc_pos_cosine(pre.processed, nK_sp, nK_pos);

%Create raised cosine spike history filters with labview cursor info
pre = load('./testdata/test_preprocess_spline_target_55_1.mat');
data = filters_sprc_relpos(pre.processed, nK_sp, nK_pos);
data = filters_sprc_pos_target(pre.processed, nK_sp, nK_pos, nK_tar);
data = filters_sprc_accel_lv(pre.processed, nK_sp, nK_pos);
data = filters_sprc_dir_lv(pre.processed, nK_sp, nK_pos);
data = filters_sprc_vel_lv(pre.processed, nK_sp, nK_pos);
data = filters_sprc_pos_target_lv(pre.processed, nK_sp, nK_pos, nK_tar);
data = filters_sprc_pos_lv(pre.processed, nK_sp, nK_pos);
data = filters_sprc_relpos_lv(pre.processed, nK_sp, nK_pos);
data = filters_sprc_revvel_lv(pre.processed, nK_sp, nK_pos);
data = filters_sprc_revpos_lv(pre.processed, nK_sp, nK_pos);
data = filters_sprc_pos_cosine_lv(pre.processed, nK_sp, nK_pos);
data = filters_sprc_pos_network_lv(pre.processed, nK_sp, nK_pos);

%Create raised cosine spike history filters with network models
data = filters_sprc_pos_network(pre.processed, nK_sp, nK_pos);