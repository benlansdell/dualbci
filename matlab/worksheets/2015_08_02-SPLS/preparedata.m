%Prepare data for analysis in R
nevfile = './testdata/20130117SpankyUtah001.nev';
matfile = './testdata/Spanky_2013-01-17-1325.mat';
binsize = 1/60;
offset = -0.075;
threshold = 3;
verbose = 0;
processed = preprocess_spline_lv(nevfile, matfile, binsize, threshold, offset, fn_out, verbose);

fn_out = './worksheets/2015_08_02-SPLS/SPLS.eps';
nK_sp = 0; 
nK_pos = 6;
dt_sp = 1/60;
dt_pos = 1/60;
data = filters_sp_pos_network_lv(processed, nK_sp, nK_pos, dt_sp, dt_pos);
data.y = data.y';
%Save data
%           X: [33366x36 double]
%           k: {26x3 cell}
%           y: [24x33366 double]
%      torque: [33366x2 double]
%     dtorque: [33366x2 double]
%    ddtorque: [33366x2 double]
%      cursor: [33366x2 double]
%     dcursor: [33366x2 double]
%    ddcursor: [33366x2 double]

save('./worksheets/2015_08_02-SPLS/20130117SpankyUtah001.mat', 'data')

%Run script in R


%Reimport data

%Plot results

