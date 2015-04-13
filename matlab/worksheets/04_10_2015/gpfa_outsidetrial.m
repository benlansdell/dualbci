%Make data structure
runIdx = 8;
fn_out = './worksheets/04_10_2015/gpfa_outsidetrial';

nevfile = './testdata/20130117SpankyUtah001.nev';
labviewfile = './testdata/Spanky_2013-01-17-1325.mat';
binsize = 0.001;
offset = 0.0;
threshold = 5;
processed = preprocess_spline_target(nevfile, labviewfile, binsize, threshold, offset);

spikes = processed.binnedspikes;
%Time bins which occur within a trial
trials = sum(abs(processed.target),2)==0;
dat = makedata(trials, spikes, processed.cursor, processed.target);

%Split into quadrants
datquads = {};
nelem = zeros(4,1);
for idx = 1:length(dat)
	q = dat(idx).quadrant;
	el = nelem(q)+1;
	datquads{q}(el).trialId = el;
	datquads{q}(el).quadrant = q;
	datquads{q}(el).spikes = dat(idx).spikes;
	nelem(q) = el;
end

method = 'gpfa';
% Select number of latent dimensions
xDim = 8;
kernSD = 30;
% Extract neural trajectories
result = neuralTraj(runIdx, dat, 'method', method, 'xDim', xDim, 'kernSDList', kernSD);
% Orthonormalize neural trajectories
[estParams, seqTrain, seqTest, DD] = postprocess(result, 'kernSD', kernSD);
% Plot each dimension of neural trajectories versus time
plotEachDimVsTime(seqTrain, 'xorth', result.binWidth);
saveplot(gcf, [fn_out '_EachDimVsTime.eps'], 'eps', [10 6])
plot3D(seqTrain, 'xorth', 'dimsToPlot', 1:3, 'nPlotMax', 100);
saveplot(gcf, [fn_out '_plot3D.eps'])
plot2D(seqTrain, 'xorth', 'dimsToPlot', 1:3, 'nPlotMax', 100);
saveplot(gcf, [fn_out '_plot2D.eps'])


% ========================================================
% 2) Full cross-validation to find:
%  - optimal state dimensionality for all methods
%  - optimal smoothing kernel width for two-stage methods
% ========================================================

% Select number of cross-validation folds
numFolds = 4;

% Perform cross-validation for different state dimensionalities.
% Results are saved in mat_results/runXXX/, where XXX is runIdx.
for xDim = [2 5 8]
  neuralTraj(runIdx, dat, 'method',  'pca', 'xDim', xDim, 'numFolds', numFolds);
  neuralTraj(runIdx, dat, 'method', 'ppca', 'xDim', xDim, 'numFolds', numFolds);
  neuralTraj(runIdx, dat, 'method',   'fa', 'xDim', xDim, 'numFolds', numFolds);
  neuralTraj(runIdx, dat, 'method', 'gpfa', 'xDim', xDim, 'numFolds', numFolds);
end
fprintf('\n');

% Plot prediction error versus state dimensionality.
% Results files are loaded from mat_results/runXXX/, where XXX is runIdx.
kernSD = 30; % select kernSD for two-stage methods
plotPredErrorVsDim(runIdx, kernSD);

% Plot prediction error versus kernelSD.
% Results files are loaded from mat_results/runXXX/, where XXX is runIdx.
xDim = 5; % select state dimensionality
plotPredErrorVsKernSD(runIdx, xDim);
