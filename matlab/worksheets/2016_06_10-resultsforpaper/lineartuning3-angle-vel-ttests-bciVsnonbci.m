conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%Tuning angles, BCI units (velocity)
bci_data = fetch(exec(conn, ['SELECT (fl1.dir+fl3.dir)/2, fl2.dir, fl5.dir, et1.`tuning_type`, rec1.`successrate`, rec2.`successrate`, '...
' IF(EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et1.`1DBCrecording` AND bci.unit = flin1.unit), '...
' 1, 0) FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`manualrecording`'...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'INNER JOIN `fits` flin2 '...
'ON flin2.`nev file` = et1.`1DBCrecording`'...
'INNER JOIN `fits_linear` fl2 '...
'ON flin2.id = fl2.id '...
'INNER JOIN `fits` flin3 '...
'ON flin3.`nev file` = et1.`manualrecordingafter`'...
'INNER JOIN `fits_linear` fl3 '...
'ON flin3.id = fl3.id '...
'INNER JOIN `fits` flin5 '...
'ON flin5.`nev file` = et1.`dualrecording`'...
'INNER JOIN `fits_linear` fl5 '...
'ON flin5.id = fl5.id '...
'INNER JOIN `recordings` rec1 '...
'ON rec1.`nev file` = et1.`1DBCrecording` '...
'INNER JOIN `recordings` rec2 '...
'ON rec2.`nev file` = et1.`dualrecording` '...
'WHERE flin1.modelID = 30 AND flin2.modelID = 30 AND flin3.modelID = 30 AND flin5.modelID = 30 ' ...
'AND flin1.unit = flin2.unit AND flin2.unit = flin3.unit AND flin2.unit = flin5.unit ' ...
'AND fl1.r2 > .01 AND fl3.r2 > .01']));

all_r2 = cell2mat(bci_data.Data(:,1:3));
bcituningtype = cell2mat(bci_data.Data(:,4));
bciperformance = 50*cell2mat(bci_data.Data(:,5:6))+10;
bciunit = cell2mat(bci_data.Data(:,7));

rot = (bcituningtype == 1 | bcituningtype == 3| bcituningtype == 4);
bci = (bciunit == 1) & rot;
nonbci = (bciunit == 0) & rot;

dof_bci = sum(bci)-1;
dof_nonbci = sum(nonbci)-1;

%1 = mc1 
%2 = bc
%3 = dc 

mc1bc_difftheta = abs(center_angles(all_r2(rot,1), all_r2(rot,2)));
mcdc_difftheta = abs(center_angles(all_r2(rot,1), all_r2(rot,3)));
bcdc_difftheta = abs(center_angles(all_r2(rot,2), all_r2(rot,3)));

bci_mc1bc_difftheta = abs(center_angles(all_r2(bci,1), all_r2(bci,2)));
bci_mcdc_difftheta = abs(center_angles(all_r2(bci,1), all_r2(bci,3)));
bci_bcdc_difftheta = abs(center_angles(all_r2(bci,2), all_r2(bci,3)));

nonbci_mc1bc_difftheta = abs(center_angles(all_r2(nonbci,1), all_r2(nonbci,2)));
nonbci_mcdc_difftheta = abs(center_angles(all_r2(nonbci,1), all_r2(nonbci,3)));
nonbci_bcdc_difftheta = abs(center_angles(all_r2(nonbci,2), all_r2(nonbci,3)));

bcidifftheta = [mean(bci_mc1bc_difftheta), mean(nonbci_mc1bc_difftheta), mean(bci_mcdc_difftheta), mean(nonbci_mcdc_difftheta)];
stdbcidifftheta = [std(bci_mc1bc_difftheta), std(nonbci_mc1bc_difftheta), std(bci_mcdc_difftheta), std(nonbci_mcdc_difftheta)];

figure 
bar(180*bcidifftheta/pi)
hold on 
errorbar(180*bcidifftheta/pi, 180*stdbcidifftheta/pi)
saveplot(gcf, './worksheets/2016_06_10-resultsforpaper/tuningangleBCI-nonBCI-bargraph.eps')

figure 
subplot(1,2,1)
grp = [zeros(size(bci_mc1bc_difftheta)); ones(size(nonbci_mc1bc_difftheta))];
boxplot([180*bci_mc1bc_difftheta/pi; 180*nonbci_mc1bc_difftheta/pi], grp)
ylim([0 180])

subplot(1,2,2)
grp = [zeros(size(bci_mcdc_difftheta)); ones(size(nonbci_mcdc_difftheta))];
boxplot([180*bci_mcdc_difftheta/pi; 180*nonbci_mcdc_difftheta/pi], grp)
ylim([0 180])

mu_bci_mc1bc_difftheta = mean(bci_mc1bc_difftheta);
mu_bci_mcdc_difftheta = mean(bci_mcdc_difftheta);
mu_bci_bcdc_difftheta = mean(bci_bcdc_difftheta);
mu_nonbci_mc1bc_difftheta  = mean(nonbci_mc1bc_difftheta);
mu_nonbci_mcdc_difftheta = mean(nonbci_mcdc_difftheta);
mu_nonbci_bcdc_difftheta = mean(nonbci_bcdc_difftheta);

[h_mc1bc, p_mc1bc] = ttest2(bci_mc1bc_difftheta, nonbci_mc1bc_difftheta, 0.05)
[h_mcdc, p_mcdc] = ttest2(bci_mcdc_difftheta, nonbci_mcdc_difftheta, 0.05)
[h_bcdc, p_bcdc] = ttest2(bci_bcdc_difftheta, nonbci_bcdc_difftheta, 0.05)

%[h_mc1bc_bci, p_mc1bc_bci] = ttest2(bci_mc1mc2_difftheta, bci_mc1bc_difftheta, 0.05)
%[h_mcdc_bci, p_mcdc_bci] = ttest2(bci_mc1mc2_difftheta, bci_mcdc_difftheta, 0.05)
%[h_bcdc_bci, p_bcdc_bci] = ttest2(bci_mc1mc2_difftheta, bci_bcdc_difftheta, 0.05)

%[h_mc1bc_nonbci, p_mc1bc_nonbci] = ttest2(nonbci_mc1mc2_difftheta, nonbci_mc1bc_difftheta, 0.05);
%[h_mcdc_nonbci, p_mcdc_nonbci] = ttest2(nonbci_mc1mc2_difftheta, nonbci_mcdc_difftheta, 0.05);
%[h_bcdc_nonbci, p_bcdc_nonbci] = ttest2(nonbci_mc1mc2_difftheta, nonbci_bcdc_difftheta, 0.05);

saveplot(gcf, './worksheets/2016_06_10-resultsforpaper/tuningangleBCI-nonBCI-boxplot.eps')


corrsbci(1) = corr(all_r2(bci,1), translate_angles(all_r2(bci,1), all_r2(bci,2)));
corrsbci(2) = corr(all_r2(bci,1), translate_angles(all_r2(bci,1), all_r2(bci,3)));
corrsbci(3) = corr(all_r2(bci,2), translate_angles(all_r2(bci,2), all_r2(bci,3)));

corrsnonbci(1) = corr(all_r2(nonbci,1), translate_angles(all_r2(nonbci,1), all_r2(nonbci,2)));
corrsnonbci(2) = corr(all_r2(nonbci,1), translate_angles(all_r2(nonbci,1), all_r2(nonbci,3)));
corrsnonbci(3) = corr(all_r2(nonbci,2), translate_angles(all_r2(nonbci,2), all_r2(nonbci,3)));

clf
cc = ones(size(bci));
cc(bci) = 1;
cc(nonbci) = 2;
colors = [1 0 0; 0 0 1];
c = [];
for idx = 1:size(cc, 1);
	c(idx,1:3) = colors(cc(idx),:);
end
colormap(colors)
subplot(2,2,1)
subplot(2,3,1)
scatter(all_r2(rot,1), translate_angles(all_r2(rot,1), all_r2(rot,2)), [], cc(rot,:))
xlabel('\theta MC1')
ylabel('\theta BC1')
title(['corr bci: ' num2str(corrsbci(1)) ' corr nonbci: ' num2str(corrsnonbci(1))])
subplot(2,2,2)
scatter(all_r2(rot,1), translate_angles(all_r2(rot,1), all_r2(rot,3)), [], cc(rot,:))
xlabel('\theta MC1')
ylabel('\theta DC')
title(['corr bci: ' num2str(corrsbci(2)) ' corr nonbci: ' num2str(corrsnonbci(2))])
subplot(2,2,3)
scatter(all_r2(rot,2), translate_angles(all_r2(rot,2), all_r2(rot,3)), [], cc(rot,:))
xlabel('\theta BC1')
ylabel('\theta DC')
title(['corr bci: ' num2str(corrsbci(3)) ' corr nonbci: ' num2str(corrsnonbci(3))])
saveplot(gcf, './worksheets/2016_06_10-resultsforpaper/tuningangleVR2R2out.eps', 'eps', [10 6])
