conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%Tuning angles, BCI units (velocity)
bci_data = fetch(exec(conn, ['SELECT fl1.size, fl2.size, fl3.size, fl5.size, et1.`tuning_type`, rec1.`successrate`, rec2.`successrate`, '...
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

all_r2 = cell2mat(bci_data.Data(:,1:4));
bcituningtype = cell2mat(bci_data.Data(:,5));
bciperformance = 50*cell2mat(bci_data.Data(:,6:7))+10;
bciunit = cell2mat(bci_data.Data(:,8));

unrot = (bcituningtype == 5);
bci = (bciunit == 1) & unrot;
nonbci = (bciunit == 0) & unrot;

dof_bci = sum(bci)-1;
dof_nonbci = sum(nonbci)-1;

%1 = mc1 
%2 = bc
%3 = dc 

mc1bc_diffsize = all_r2(unrot,1) - all_r2(unrot,2);
mc1mc2_diffsize = all_r2(unrot,1)- all_r2(unrot,3);
mc2bc_diffsize = all_r2(unrot,3)- all_r2(unrot,2);
mcdc_diffsize = all_r2(unrot,1)- all_r2(unrot,4);
bcdc_diffsize = all_r2(unrot,2)- all_r2(unrot,4);

[h_mc1bc_bci, p_mc1bc_bci] = ttest2(mc1mc2_diffsize, mc1bc_diffsize, 0.05)
[h_mcdc_bci, p_mcdc_bci] = ttest2(mc1mc2_diffsize, mcdc_diffsize, 0.05)
[h_bcdc_bci, p_bcdc_bci] = ttest2(mc1mc2_diffsize, bcdc_diffsize, 0.05)

figure 
bar([mean(abs(mc1mc2_diffsize)), mean(abs(mc1bc_diffsize)), mean(abs(mcdc_diffsize))])
hold on 
errorbar([mean(abs(mc1mc2_diffsize)), mean(abs(mc1bc_diffsize)), mean(abs(mcdc_diffsize))], ...
	[std(abs(mc1mc2_diffsize)), std(abs(mc1bc_diffsize)), std(abs(mcdc_diffsize))])
ylabel('\Delta size')
%ylim([0 120])
saveplot(gcf, './worksheets/2016_06_10-resultsforpaper/tuningstrength-conditions-bargraph-unrotated.eps')


%[h_mc1bc_nonbci, p_mc1bc_nonbci] = ttest2(nonbci_mc1mc2_difftheta, nonbci_mc1bc_difftheta, 0.05);
%[h_mcdc_nonbci, p_mcdc_nonbci] = ttest2(nonbci_mc1mc2_difftheta, nonbci_mcdc_difftheta, 0.05);
%[h_bcdc_nonbci, p_bcdc_nonbci] = ttest2(nonbci_mc1mc2_difftheta, nonbci_bcdc_difftheta, 0.05);

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
scatter(all_r2(unrot,1), all_r2(unrot,2), [], cc(unrot))
hold on 
plot([0 3], [0 3], 'k')
xlim([0 3])
ylim([0 3])
xlabel('strength MC1')
ylabel('strength BC1')
%title(['corr bci: ' num2str(corrsbci(1)) ' corr nonbci: ' num2str(corrsnonbci(1))])
subplot(2,2,2)
scatter(all_r2(unrot,1), all_r2(unrot,3), [], cc(unrot))
hold on 
plot([0 3], [0 3], 'k')
xlim([0 3])
ylim([0 3])
xlabel('strength MC1')
ylabel('strength MC2')
%title(['corr bci: ' num2str(corrsbci(2)) ' corr nonbci: ' num2str(corrsnonbci(2))])
subplot(2,2,3)
scatter(all_r2(unrot,2), all_r2(unrot,4), [], cc(unrot))
hold on 
plot([0 3], [0 3], 'k')
xlim([0 3])
ylim([0 3])
xlabel('strength BC1')
ylabel('strength DC')
subplot(2,2,4)
scatter(all_r2(unrot,1), all_r2(unrot,4), [], cc(unrot))
hold on 
plot([0 3], [0 3], 'k')
xlim([0 3])
ylim([0 3])
xlabel('strength MC1')
ylabel('strength DC')
%title(['corr bci: ' num2str(corrsbci(3)) ' corr nonbci: ' num2str(corrsnonbci(3))])
saveplot(gcf, './worksheets/2016_06_10-resultsforpaper/tuningstrength-conditions-scatter-unrotated.eps')
