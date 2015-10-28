%Make data structure
runIdx = 19;
nevfile = './testdata/20140116SpankyUtah001.nev';
fn_out = './worksheets/04_21_2015/gpfa_centerout_2014-01-16';
labviewfile = './testdata/20140116_ExperimentalParamCorrected.mat';

binsize = 0.001;
offset = 0.0;
threshold = 5;
processed_bci = preprocess_centerout(nevfile, labviewfile, binsize, threshold, offset);

%Time bins which occur within a trial
spikes_bci = processed_bci.binnedspikes;
trials_bci = processed_bci.target>0;

%Run all trials together
[dat_bci, targ] = makedata_centerout(trials_bci, spikes_bci, processed_bci.target);

method = 'gpfa';
% Select number of latent dimensions
xDim = 8;
kernSD = 30;
% Extract neural trajectories
result_bci = neuralTraj(runIdx, dat_bci, 'method', method, 'xDim', xDim, 'kernSDList', kernSD);
% Orthonormalize neural trajectories
[estParams_bci, seqTrain_bci, seqTest_bci, DD] = postprocess(result_bci, 'kernSD', kernSD);
%Power
pow = 100*diag(DD.^2)/sum(diag(DD.^2));

% Plot each dimension of neural trajectories versus time
plotEachDimVsTime_coloct(seqTrain_bci, 'xorth', result_bci.binWidth, targ, 'nPlotMax', 100);
saveplot(gcf, [fn_out '_EachDimVsTime_coltarg.eps'], 'eps', [15 15])
plotEachDimVsTime_colmeantarget(seqTrain_bci, 'xorth', result_bci.binWidth, targ, 'nPlotMax', 10000);
plot2svg([fn_out '_EachDimVsTime_colmeantarget.svg'])
plot3D_coloct(seqTrain_bci, 'xorth', targ, 'dimsToPlot', 1:3, 'nPlotMax', 100);
plot2svg([fn_out '_plot3D_coloct.svg'])
%saveplot(gcf, [fn_out '_plot3D_coloct.eps'])
plot2D_coloct(seqTrain_bci, 'xorth', targ, 'dimsToPlot', 1:3, 'nPlotMax', 100);
plot2svg([fn_out '_plot2D_coloct.svg'])
%saveplot(gcf, [fn_out '_plot2D_coloct.eps'])

%saveplot(gcf, [fn_out '_plot3D_coloct.eps'])
plot2D_coloct_mean(seqTrain_bci, 'xorth', targ, 'dimsToPlot', 1:3, 'nPlotMax', 10000);
plot2svg([fn_out '_plot2D_coloct_mean.svg'])

figure
hold on
cm = lines(3);
for idx = 1:3
	%plot((estParams_bci.Corth(:,1:3)), '-o')
	x = idx*0.05+(1:size(estParams_bci.Corth, 1));
	stem(x,(estParams_bci.Corth(:,idx)), 'Color', cm(idx,:))
end
%label according to unit
nU = length(processed_bci.unitnames);
set(gca,'XTick',1:nU);
set(gca,'XTickLabel',processed_bci.unitnames);
legend(['PC1. ' num2str(pow(1)) '% energy'], ['PC2. ' num2str(pow(2)) '% energy'], ['PC3. ' num2str(pow(3)) '% energy'])
saveplot(gcf, [fn_out '_basisvecs.eps'], 'eps', [10 7])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

runIdx = 20;
nevfile = './testdata/20140203SpankyUtah001.nev';
fn_out = './worksheets/04_21_2015/gpfa_centerout_2014-02-03';
labviewfile = './testdata/20140203_ExperimentalParamCorrected.mat';

binsize = 0.001;
offset = 0.0;
threshold = 5;
processed_bci = preprocess_centerout(nevfile, labviewfile, binsize, threshold, offset);

%Time bins which occur within a trial
spikes_bci = processed_bci.binnedspikes;
trials_bci = processed_bci.target>0;

%Run all trials together
[dat_bci, targ] = makedata_centerout(trials_bci, spikes_bci, processed_bci.target);

method = 'gpfa';
% Select number of latent dimensions
xDim = 8;
kernSD = 30;
% Extract neural trajectories
result_bci = neuralTraj(runIdx, dat_bci, 'method', method, 'xDim', xDim, 'kernSDList', kernSD);
% Orthonormalize neural trajectories
[estParams_bci, seqTrain_bci, seqTest_bci, DD] = postprocess(result_bci, 'kernSD', kernSD);

%Power
pow = 100*diag(DD.^2)/sum(diag(DD.^2));

% Plot each dimension of neural trajectories versus time
plotEachDimVsTime_coloct(seqTrain_bci, 'xorth', result_bci.binWidth, targ, 'nPlotMax', 100);
saveplot(gcf, [fn_out '_EachDimVsTime_coltarg.eps'], 'eps', [15 15])
plotEachDimVsTime_colmeantarget(seqTrain_bci, 'xorth', result_bci.binWidth, targ, 'nPlotMax', 10000);
plot2svg([fn_out '_EachDimVsTime_colmeantarget.svg'])
plot3D_coloct(seqTrain_bci, 'xorth', targ, 'dimsToPlot', 1:3, 'nPlotMax', 100);
plot2svg([fn_out '_plot3D_coloct.svg'])
%saveplot(gcf, [fn_out '_plot3D_coloct.eps'])
plot2D_coloct(seqTrain_bci, 'xorth', targ, 'dimsToPlot', 1:3, 'nPlotMax', 100);
plot2svg([fn_out '_plot2D_coloct.svg'])
%saveplot(gcf, [fn_out '_plot2D_coloct.eps'])

%saveplot(gcf, [fn_out '_plot3D_coloct.eps'])
plot2D_coloct_mean(seqTrain_bci, 'xorth', targ, 'dimsToPlot', 1:3, 'nPlotMax', 100);
plot2svg([fn_out '_plot2D_coloct_mean.svg'])


figure
hold on
cm = lines(3);
for idx = 1:3
	%plot((estParams_bci.Corth(:,1:3)), '-o')
	x = idx*0.05+(1:size(estParams_bci.Corth, 1));
	stem(x,(estParams_bci.Corth(:,idx)), 'Color', cm(idx,:))
end
%label according to unit
nU = length(processed_bci.unitnames);
set(gca,'XTick',1:nU);
set(gca,'XTickLabel',processed_bci.unitnames);
legend(['PC1. ' num2str(pow(1)) '% energy'], ['PC2. ' num2str(pow(2)) '% energy'], ['PC3. ' num2str(pow(3)) '% energy'])
saveplot(gcf, [fn_out '_basisvecs.eps'], 'eps', [10 7])