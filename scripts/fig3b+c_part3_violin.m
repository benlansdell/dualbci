load('./scripts/fig3b+c_part2.mat')

noncontrol_dcotunedMCBC = dcotunedMCBC;
noncontrol_dotherMCBC = dotherMCBC;
noncontrol_dcotunedMCDC = dcotuned;
noncontrol_dotherMCDC = dother;

load('./scripts/fig3b+c_part1.mat')

control_dcotunedMCBC = dcotunedMCBC;
control_dotherMCBC = dotherMCBC;
control_dcotunedMCDC = dcotuned;
control_dotherMCDC = dother;

%control = load('bcitTEchanges_bootstrap_rotated.mat')
%noncontrol = load('bcitTEchanges_bootstrap_rotated_noncontrol.mat')

%Combine
dcontrolMCBC = [control_dcotunedMCBC; control_dotherMCBC];
dnoncontrolMCBC = [noncontrol_dcotunedMCBC; noncontrol_dotherMCBC];
dcontrolMCDC = [control_dcotunedMCDC; control_dotherMCDC];         
dnoncontrolMCDC = [noncontrol_dcotunedMCDC; noncontrol_dotherMCDC];

[pMCBCrs, hMCBCrs] = ranksum((dcontrolMCBC), (dnoncontrolMCBC))
[pMCDCrs, hMCDCrs] = ranksum((dcontrolMCDC), (dnoncontrolMCDC))

[pMCBCsrC, hMCBCsrC] = signrank(dcontrolMCBC)
[pMCBCsrN, hMCBCsrN] = signrank(dnoncontrolMCBC)
[pMCDCsrC, hMCDCsrC] = signrank(dcontrolMCDC)
[pMCDCsrN, hMCDCsrN] = signrank(dnoncontrolMCDC)

figure
subplot(2,2,1)
violinplot(dcontrolMCBC, 'c1', 'ShowData', false);
title('MC-BC control')
ylim([-0.001, 0.0002])
subplot(2,2,2)
violinplot(dnoncontrolMCBC, 'c1', 'ShowData', false);
title('MC-BC non control')
ylim([-0.001, 0.0002])
subplot(2,2,3)
violinplot(dcontrolMCDC, 'c1', 'ShowData', false);
title('MC-BC control')
ylim([-0.001, 0.0002])
subplot(2,2,4)
violinplot(dnoncontrolMCDC, 'c1', 'ShowData', false);
title('MC-BC non control')
ylim([-0.001, 0.0002])

saveplot(gcf, './figures/TE-decoupling-violin_bootstrap_rotated_signed_sem_control+noncontrol.eps')

figure
groups = [ones(size(dcontrolMCBC)); 2*ones(size(dnoncontrolMCBC));...
			3*ones(size(dcontrolMCDC));4*ones(size(dnoncontrolMCDC))];

h = boxplot([dcontrolMCBC; dnoncontrolMCBC; dcontrolMCDC; dnoncontrolMCDC], groups);
set(h(7,:), 'Visible', 'off')
ylim([-0.0005, 0.0005])

title(['(rank sum) MCBC: (dc)vs (do) p-value: ' num2str(pMCBCrs) ', MCDC: (dc)vs (do) p-value: ' num2str(pMCDCrs)])
saveplot(gcf, './figures/TE-decoupling-boxplot_bootstrap_rotated_signed_sem_control+noncontrol.eps')
