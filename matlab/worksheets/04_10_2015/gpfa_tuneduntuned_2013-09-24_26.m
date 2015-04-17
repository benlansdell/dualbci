%Make data structure

%%%%%%%%%%%%%%%%%
%Untuned control%
%%%%%%%%%%%%%%%%%
runIdx = 14;
nevfile = './testdata/20130924SpankyUtah004.nev';
fn_out = './worksheets/04_10_2015/gpfa_bci_2013-09-24';
labviewfile = './testdata/Spanky_2013-09-24-1318.mat';
binsize = 0.001;
offset = 0.0;
threshold = 5;
processed_bci = preprocess_spline_target(nevfile, labviewfile, binsize, threshold, offset);

%Time bins which occur within a trial
spikes_bci = processed_bci.binnedspikes;
trials_bci = sum(abs(processed_bci.target),2)>0;

%Run all trials together
[dat_bci, octs_bci, quads_bci] = makedata(trials_bci, spikes_bci, processed_bci.cursor, processed_bci.target);
method = 'gpfa';
% Select number of latent dimensions
xDim = 8;
kernSD = 30;
% Extract neural trajectories
result_bci = neuralTraj(runIdx, dat_bci, 'method', method, 'xDim', xDim, 'kernSDList', kernSD);
% Orthonormalize neural trajectories
[estParams_bci, seqTrain_bci, seqTest_bci, DD] = postprocess(result_bci, 'kernSD', kernSD);
% Plot each dimension of neural trajectories versus time
plotEachDimVsTime_colhoriz(seqTrain_bci, 'xorth', result_bci.binWidth, quads_bci, 'nPlotMax', 100);
saveplot(gcf, [fn_out '_EachDimVsTime_colhoriz.eps'], 'eps', [10 6])
plotEachDimVsTime_specificquad(seqTrain_bci, 'xorth', result_bci.binWidth, quads, 'nPlotMax', 100);
saveplot(gcf, [fn_out '_EachDimVsTime_specificquad.eps'], 'eps', [20 6])
plot3D_coloct(seqTrain_bci, 'xorth', octs, 'dimsToPlot', 1:3, 'nPlotMax', 100);
saveplot(gcf, [fn_out '_plot3D_coloct.eps'])
plot2D_coloct(seqTrain_bci, 'xorth', octs, 'dimsToPlot', 1:3, 'nPlotMax', 100);
saveplot(gcf, [fn_out '_plot2D_coloct.eps'])
figure
plot(abs(estParams_bci.Corth(:,1:3)), '-o')
%label according to unit
nU = length(processed_bci.unitnames);
set(gca,'XTick',1:nU);
set(gca,'XTickLabel',processed_bci.unitnames);
legend('|PC1|', '|PC2|', '|PC3|')
saveplot(gcf, [fn_out '_basisvecs.eps'], 'eps', [10 7])


%%%%%%%%%%%%%%%
%Tuned control%
%%%%%%%%%%%%%%%

runIdx = 15;
fn_out = './worksheets/04_10_2015/gpfa_manual_2013-09-26';
nevfile = './testdata/20130926SpankyUtah014.nev';
labviewfile = './testdata/Spanky_2013-09-26-1334.mat';
binsize = 0.001;
offset = 0.0;
threshold = 5;
processed_man = preprocess_spline_target(nevfile, labviewfile, binsize, threshold, offset);

%Time bins which occur within a trial
spikes_man = processed_man.binnedspikes;
trials_man = sum(abs(processed_man.target),2)>0;

%Run all trials together
[dat_man, octs_man, quads_man] = makedata(trials_man, spikes_man, processed_man.cursor, processed_man.target);
method = 'gpfa';
% Select number of latent dimensions
xDim = 8;
kernSD = 30;
% Extract neural trajectories
result_man = neuralTraj(runIdx, dat_man, 'method', method, 'xDim', xDim, 'kernSDList', kernSD);
% Orthonormalize neural trajectories
[estParams_man, seqTrain_man, seqTest, DD] = postprocess(result_man, 'kernSD', kernSD);
% Plot each dimension of neural trajectories versus time
plotEachDimVsTime_colhoriz(seqTrain_man, 'xorth', result_man.binWidth, quads_man, 'nPlotMax', 100);
saveplot(gcf, [fn_out '_EachDimVsTime_colhoriz.eps'], 'eps', [10 6])
plotEachDimVsTime_specificquad(seqTrain_man, 'xorth', result_man.binWidth, quads_man, 'nPlotMax', 100);
saveplot(gcf, [fn_out '_EachDimVsTime_specificquad.eps'], 'eps', [20 6])
plot3D_coloct(seqTrain_man, 'xorth', octs_man, 'dimsToPlot', 1:3, 'nPlotMax', 100);
saveplot(gcf, [fn_out '_plot3D_coloct.eps'])
plot2D_coloct(seqTrain_man, 'xorth', octs_man, 'dimsToPlot', 1:3, 'nPlotMax', 100);
saveplot(gcf, [fn_out '_plot2D_coloct.eps'])

figure
plot(abs(estParams_man.Corth(:,1:3)), '-o')
%label according to unit
nU = length(processed_man.unitnames);
set(gca,'XTick',1:nU);
set(gca,'XTickLabel',processed_man.unitnames);
legend('|PC1|', '|PC2|', '|PC3|')
saveplot(gcf, [fn_out '_basisvecs.eps'], 'eps', [10 7])
