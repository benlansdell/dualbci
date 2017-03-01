control = load('bcitTEchanges_bootstrap_rotated.mat')
noncontrol = load('bcitTEchanges_bootstrap_rotated_noncontrol.mat')

%Combine
dcontrolMCBC = [control.dcotunedMCBC; control.dotherMCBC];
dnoncontrolMCBC = [noncontrol.dcotunedMCBC; noncontrol.dotherMCBC];
dcontrolMCDC = [control.dcotunedMCDC; control.dotherMCDC];         
dnoncontrolMCDC = [noncontrol.dcotunedMCDC; noncontrol.dotherMCDC];

[pMCBCrs, hMCBCrs] = ranksum((dcontrolMCBC), (dnoncontrolMCBC))
[pMCDCrs, hMCDCrs] = ranksum((dcontrolMCDC), (dnoncontrolMCDC))

[pMCBCsrC, hMCBCsrC] = signrank(dcontrolMCBC)
[pMCBCsrN, hMCBCsrN] = signrank(dnoncontrolMCBC)
[pMCDCsrC, hMCDCsrC] = signrank(dcontrolMCDC)
[pMCDCsrN, hMCDCsrN] = signrank(dnoncontrolMCDC)

figure
bar([mean((dcontrolMCBC)), mean((dnoncontrolMCBC)), mean((dcontrolMCDC)), mean((dnoncontrolMCDC))]);
hold on 
title(['MCBC: (dc)vs (do) p-value: ' num2str(pMCBCrs) ', MCDC: (dc)vs (do) p-value: ' num2str(pMCDCrs)])

errorbar([mean((dcontrolMCBC)), mean((dnoncontrolMCBC)), mean((dcontrolMCDC)), mean((dnoncontrolMCDC))],...
	[std(dcontrolMCBC)/sqrt(length(dcontrolMCBC)), std(dnoncontrolMCBC)/sqrt(length(dnoncontrolMCBC)), std(dcontrolMCDC)/sqrt(length(dcontrolMCDC)), std(dnoncontrolMCDC)/sqrt(length(dnoncontrolMCDC))]);
ylim([-0.00018, 0])

saveplot(gcf, './worksheets/2016_06_10-resultsforpaper/TE-decoupling-bargraph_bootstrap_rotated_signed_sem_control+noncontrol.eps')

figure
groups = [ones(size(dcontrolMCBC)); 2*ones(size(dnoncontrolMCBC));...
			3*ones(size(dcontrolMCDC));4*ones(size(dnoncontrolMCDC))];

h = boxplot([dcontrolMCBC; dnoncontrolMCBC; dcontrolMCDC; dnoncontrolMCDC], groups);
set(h(7,:), 'Visible', 'off')
ylim([-0.0005, 0.0005])

title(['(rank sum) MCBC: (dc)vs (do) p-value: ' num2str(pMCBCrs) ', MCDC: (dc)vs (do) p-value: ' num2str(pMCDCrs)])
saveplot(gcf, './worksheets/2016_06_10-resultsforpaper/TE-decoupling-boxplot_bootstrap_rotated_signed_sem_control+noncontrol.eps')
