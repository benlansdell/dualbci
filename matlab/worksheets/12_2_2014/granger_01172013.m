%%%%%%%%%%%%%%%%
%Manual control%
%%%%%%%%%%%%%%%%

%Load test preprocessed data
pre = load('./testdata/test_preprocess_spline_60hz_short24.mat');
processed = pre.processed;

%Run with position filters
const = 'on';
pval = 0.001;
nK_sp = 6; 
nK_pos = 6;
fn_out = './worksheets/12_2_2014/plots/granger_manual_pos.eps';
data = filters_sp_pos_network(processed, nK_sp, nK_pos);
[GCdevMP, GCpvalMP, GCsigMP] = granger(processed, data, fn_out, pval);

%Run without position filters
const = 'on';
pval = 0.001;
nK_sp = 6; 
nK_pos = 0;
fn_out = './worksheets/12_2_2014/plots/granger_manual_wout_pos.eps';
data = filters_sp_pos_network(processed, nK_sp, nK_pos);
[GCdevM, GCpvalM, GCsigM] = granger(processed, data, fn_out, pval);

%%%%%%%%%%%%%%%
%Brain control%
%%%%%%%%%%%%%%%

%nevfile = './testdata/20130117SpankyUtah005.nev';
%binsize = 1/60;
%offset = 0.0;
%threshold = 5;
%processed = preprocess_spline(nevfile, binsize, threshold, offset);
%Modify so that it's shorter, then save it
%save('./testdata/test_preprocess_brain_spline_60hz_short24.mat', 'processed');
pre = load('./testdata/test_preprocess_brain_spline_60hz_short24.mat');
processed = pre.processed;
%Run on brain control data with position filters
const = 'on';
pval = 0.001;
nK_sp = 6; 
nK_pos = 6;
fn_out = './worksheets/12_2_2014/plots/granger_brain_pos.eps';
data = filters_sp_pos_network(processed, nK_sp, nK_pos);
[GCdevBP, GCpvalBP, GCsigBP] = granger(processed, data, fn_out, pval);

%Run on brain control data without position filters
const = 'on';
pval = 0.001;
nK_sp = 6; 
nK_pos = 0;
fn_out = './worksheets/12_2_2014/plots/granger_brain_wout_pos.eps';
data = filters_sp_pos_network(processed, nK_sp, nK_pos);
[GCdevB, GCpvalB, GCsigB] = granger(processed, data, fn_out, pval);

%Save results
save('./worksheets/12_2_2014/GLMGranger.mat', 'GCdevB', 'GCdevBP', 'GCdevM', 'GCdevMP', 'GCpvalB', 'GCpvalBP', 'GCpvalM', 'GCpvalMP', 'GCsigB', 'GCsigBP', 'GCsigM', 'GCsigMP');

%Compare brain control and non-brain control datasets

%Compare with Ivana's code: (Manual control, linear model on spike trains, testing based on F-stat)
clear
pre = load('./testdata/test_preprocess_spline_60hz_short24.mat');
processed = pre.processed;
load('./worksheets/12_2_2014/ParamforBen.mat');
%Add identity to sig matrix (remove NaN's from main diagonal)
nU = 24;
for i = 1:nU
	sig(i,i) = 1;
end
%Plot
clf
fn_out = './worksheets/12_2_2014/plots/IvanaGC_LM_wout_pos.eps';
subplot(3,1,1)
colormap(bone);
imagesc(F')
title('Change in deviance')
ylabel('Unit')
xlabel('Unit')
set(gca,'XTick',1:nU);
set(gca,'YTick',1:nU);
set(gca,'XTickLabel',processed.unitnames);
set(gca,'YTickLabel',processed.unitnames);
rotateXLabels(gca, 90);
colorbar
subplot(3,1,2)
imagesc(pval')
title(['p-value'])
ylabel('Unit')
xlabel('Unit')
set(gca,'XTick',1:nU);
set(gca,'YTick',1:nU);
set(gca,'XTickLabel',processed.unitnames);
set(gca,'YTickLabel',processed.unitnames);
%rotate x labels
rotateXLabels(gca, 90);
colorbar
subplot(3,1,3)
imagesc(sig')
title(['Significance test @ p<' num2str(pval)])
ylabel('Unit')
xlabel('Unit')
set(gca,'XTick',1:nU);
set(gca,'YTick',1:nU);
set(gca,'XTickLabel',processed.unitnames);
set(gca,'YTickLabel',processed.unitnames);
rotateXLabels(gca, 90);
colorbar
saveplot(gcf, fn_out, 'eps', [6 18])

%Compare LM's F-statistc to GLM's chi2-statistic to see if at least related...
clf
Fvec = reshape(F, nU*nU,1);
chi2vec = reshape(GCdev', nU*nU,1);
plot(Fvec', chi2vec', '.');
xlabel('F-stat')
ylabel('chi2 statistic')
saveplot(gcf, './worksheets/12_2_2014/plots/fstatvschi2stat.eps')

pre = load('./testdata/test_preprocess_brain_spline_60hz_short24.mat');
processed = pre.processed;
unitnames = processed.unitnames;
unitnamesB = unitnames
[GCdevpermB, GCpvalpermB, GCsigpermB, clustersB, namespermB] = granger_cluster(GCdevB, GCpvalB, GCsigB, unitnames, fn_out);
[GCdevpermBP, GCpvalpermBP, GCsigpermBP, clustersBP, namespermBP] = granger_cluster(GCdevBP, GCpvalBP, GCsigBP, unitnames, fn_out);
save('./worksheets/12_2_2014/GLMGrangerB.mat', 'GCdevB', 'GCpvalB', 'GCsigB', 'GCdevpermB', 'GCpvalpermB', 'GCsigpermB', 'clustersB', 'namespermB', 'unitnames');
save('./worksheets/12_2_2014/GLMGrangerBP.mat', 'GCdevBP', 'GCpvalBP', 'GCsigBP', 'GCdevpermBP', 'GCpvalpermBP', 'GCsigpermBP', 'clustersBP', 'namespermBP', 'unitnames');

pre = load('./testdata/test_preprocess_spline_60hz_short24.mat');
processed = pre.processed;
unitnames = processed.unitnames;
unitnamesM = unitnames
[GCdevpermM, GCpvalpermM, GCsigpermM, clustersM, namespermM] = granger_cluster(GCdevM, GCpvalM, GCsigM, unitnames, fn_out);
[GCdevpermMP, GCpvalpermMP, GCsigpermMP, clustersMP, namespermMP] = granger_cluster(GCdevMP, GCpvalMP, GCsigMP, unitnames, fn_out);
save('./worksheets/12_2_2014/GLMGrangerM.mat', 'GCdevM', 'GCpvalM', 'GCsigM', 'GCdevpermM', 'GCpvalpermM', 'GCsigpermM', 'clustersM', 'namespermM', 'unitnames');
save('./worksheets/12_2_2014/GLMGrangerMP.mat', 'GCdevMP', 'GCpvalMP', 'GCsigMP', 'GCdevpermMP', 'GCpvalpermMP', 'GCsigpermMP', 'clustersMP', 'namespermMP', 'unitnames');


%Find which units are in common between brain and manual control
nUB = length(unitnamesB);
nUM = length(unitnamesM);
MinB = [];
BinM = [];
for i = 1:nUM
	if any(ismember(unitnamesB, unitnamesM{i}))
		MinB = [MinB, i];
	end
end

for i = 1:nUB
	if any(ismember(unitnamesM, unitnamesB{i}))
		BinM = [BinM, i];
	end
end

unitnames = unitnamesB(BinM);
nU = length(unitnames);
GCdevMcom = GCdevM(MinB,MinB);
GCpvalMcom = GCpvalM(MinB,MinB);
GCsigMcom = GCsigM(MinB,MinB);
GCdevMPcom = GCdevMP(MinB,MinB);
GCpvalMPcom = GCpvalMP(MinB,MinB);
GCsigMPcom = GCsigMP(MinB,MinB);

GCdevBcom = GCdevB(BinM,BinM);
GCpvalBcom = GCpvalB(BinM,BinM);
GCsigBcom = GCsigB(BinM,BinM);
GCdevBPcom = GCdevBP(BinM,BinM);
GCpvalBPcom = GCpvalBP(BinM,BinM);
GCsigBPcom = GCsigBP(BinM,BinM);

%Plot deviance for brain vs manual, brain position vs manual position
%Plot deviance for MP vs M, and BP vs B
clf
X = reshape(GCdevBcom, nU*nU,1);
Y = reshape(GCdevBPcom, nU*nU,1);
plot(X, Y, '.');
xlabel('Brain, no pos filter')
ylabel('Brian, pos filter')
saveplot(gcf, './worksheets/12_2_2014/plots/BvsBP.eps')

clf
X = reshape(GCdevMcom, nU*nU,1);
Y = reshape(GCdevMPcom, nU*nU,1);
plot(X, Y, '.');
xlabel('Manual, no pos filter')
ylabel('Manual, pos filter')
saveplot(gcf, './worksheets/12_2_2014/plots/MvsMP.eps')

%%%%%%%%%%%%%%%

clf
X = reshape(GCdevMcom, nU*nU,1);
Y = reshape(GCdevBcom, nU*nU,1);
plot(X, Y, '.');
xlabel('Manual')
ylabel('Brain')
saveplot(gcf, './worksheets/12_2_2014/plots/MvsB.eps')

clf
X = reshape(GCdevMPcom, nU*nU,1);
Y = reshape(GCdevBPcom, nU*nU,1);
plot(X, Y, '.');
xlabel('Manual, pos filter')
ylabel('Brain, pos filter')
saveplot(gcf, './worksheets/12_2_2014/plots/MPvsBP.eps')
