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

figure
BMIN = -0.001;
BMAX = 0.001;
x = 0;

subplot(4,2,1)
histogram(control_dcotunedMCBC,'BinLimits',[BMIN,BMAX], 'Normalization', 'countdensity')
title('control cotuned MCBC')
hold on;
line([x, x], ylim, 'LineWidth', 2, 'Color', 'r');
x1 = mean(control_dcotunedMCBC);
line([x1, x1], ylim, 'LineWidth', 2, 'Color', 'k');

subplot(4,2,2)
histogram(control_dcotunedMCDC,'BinLimits',[BMIN,BMAX], 'Normalization', 'countdensity')
title('control cotuned MCDC')
hold on;
line([x, x], ylim, 'LineWidth', 2, 'Color', 'r');
x1 = mean(control_dcotunedMCDC);
line([x1, x1], ylim, 'LineWidth', 2, 'Color', 'k');

subplot(4,2,3)
histogram(noncontrol_dcotunedMCBC,'BinLimits',[BMIN,BMAX], 'Normalization', 'countdensity')
title('non control cotuned MCBC')
hold on;
line([x, x], ylim, 'LineWidth', 2, 'Color', 'r');
x1 = mean(noncontrol_dcotunedMCBC);
line([x1, x1], ylim, 'LineWidth', 2, 'Color', 'k');

subplot(4,2,4)
histogram(noncontrol_dcotunedMCDC,'BinLimits',[BMIN,BMAX], 'Normalization', 'countdensity')
title('non control cotuned MCDC')
hold on;
line([x, x], ylim, 'LineWidth', 2, 'Color', 'r');
x1 = mean(noncontrol_dcotunedMCDC);
line([x1, x1], ylim, 'LineWidth', 2, 'Color', 'k');

subplot(4,2,5)
histogram(control_dotherMCBC,'BinLimits',[BMIN,BMAX], 'Normalization', 'countdensity')
title('control random MCBC')
hold on;
line([x, x], ylim, 'LineWidth', 2, 'Color', 'r');
x1 = mean(control_dotherMCBC);
line([x1, x1], ylim, 'LineWidth', 2, 'Color', 'k');


subplot(4,2,6)
histogram(control_dotherMCDC,'BinLimits',[BMIN,BMAX], 'Normalization', 'countdensity')
title('control random MCDC')
hold on;
line([x, x], ylim, 'LineWidth', 2, 'Color', 'r');
x1 = mean(noncontrol_dotherMCDC);
line([x1, x1], ylim, 'LineWidth', 2, 'Color', 'k');

subplot(4,2,7)
histogram(noncontrol_dotherMCBC,'BinLimits',[BMIN,BMAX], 'Normalization', 'countdensity')
title('non control random MCBC')
hold on;
line([x, x], ylim, 'LineWidth', 2, 'Color', 'r');
x1 = mean(noncontrol_dotherMCBC);
line([x1, x1], ylim, 'LineWidth', 2, 'Color', 'k');

subplot(4,2,8)
histogram(noncontrol_dotherMCDC,'BinLimits',[BMIN,BMAX], 'Normalization', 'countdensity')
title('non control random MCDC')
hold on;
line([x, x], ylim, 'LineWidth', 2, 'Color', 'r');
x1 = mean(noncontrol_dotherMCDC);
line([x1, x1], ylim, 'LineWidth', 2, 'Color', 'k');

%bar([mean((dcotunedMCBC)), mean((dotherMCBC)), mean((dcotunedMCDC)), mean((dotherMCDC))]);

saveplot(gcf, './figures/fig3b_TE_changes.eps', 'eps', [8, 6])

%Plot the bar graphs