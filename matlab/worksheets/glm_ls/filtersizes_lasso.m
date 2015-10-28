nevfile = './testdata/20130117SpankyUtah001.nev';
matfile = './testdata/Spanky_2013-01-17-1325.mat';
units = {'55.1'};
threshold = 5;

binsize = 0.002;
nkt = 100;
fn_out = './worksheets/glm_ls/plots/nk_100_bs_2ms_lasso';
b = fitGLM_lasso(nevfile, matfile, fn_out, binsize, nkt, threshold, units);

binsize = 0.005;
nkt = 50;
fn_out = './worksheets/glm_ls/plots/nk_50_bs_5ms_lasso';
b = fitGLM_lasso(nevfile, matfile, fn_out, binsize, nkt, threshold, units);

binsize = 0.01;
nkt = 30;
fn_out = './worksheets/glm_ls/plots/nk_30_bs_10ms_lasso';
b = fitGLM_lasso(nevfile, matfile, fn_out, binsize, nkt, threshold, units);

binsize = 0.02;
nkt = 15;
fn_out = './worksheets/glm_ls/plots/nk_15_bs_20ms_lasso';
b = fitGLM_lasso(nevfile, matfile, fn_out, binsize, nkt, threshold, units);

binsize = 0.05;
nkt = 6;
fn_out = './worksheets/glm_ls/plots/nk_6_bs_50ms_lasso';
b = fitGLM_lasso(nevfile, matfile, fn_out, binsize, nkt, threshold, units);