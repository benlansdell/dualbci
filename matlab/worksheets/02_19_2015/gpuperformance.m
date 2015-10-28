%Test GPU performance...
const = 'on';
nK_sp = 100; 
nK_pos = 100;
pre = load('./testdata/test_preprocess_spline_short.mat');
data = filters_sp_pos(pre.processed, nK_sp, nK_pos);

tic;
model = MLE_glmfit(data, const);
tCPU = toc

tic;
modelGPU = MLE_glmfit_GPU(data, const);
tGPU = toc
