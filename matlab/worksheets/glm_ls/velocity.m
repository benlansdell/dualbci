nevfile = './testdata/20130117SpankyUtah001.nev';
matfile = './testdata/Spanky_2013-01-17-1325.mat';
units = {'55.1'};
threshold = 5;

binsize = 0.002;
nkt = 100;
fn_out = './worksheets/glm_ls/plots/nk_100_bs_2ms_vel';
b = fitGLM_LS_vel(nevfile, matfile, fn_out, binsize, nkt, threshold, units);