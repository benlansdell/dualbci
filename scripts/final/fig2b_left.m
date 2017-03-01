conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%Tuning angles, BCI units (velocity)
bci_data = fetch(exec(conn, ['SELECT fl1.dir, fl2.dir, fl3.dir, fl5.dir, et1.`tuning_type`, rec1.`successrate`, rec2.`successrate`, '...
' IF(EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et1.`1DBCrecording` AND bci.unit = flin1.unit), '...
' 1, 0), flin1.`nev file` FROM '...
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
'AND flin1.unit = flin2.unit AND flin2.unit = flin3.unit AND flin2.unit = flin5.unit '...
'AND fl1.r2 > .01 AND fl3.r2 > .01']));
all_r2 = cell2mat(bci_data.Data(:,1:4));
bcituningtype = cell2mat(bci_data.Data(:,5));
bciperformance = 50*cell2mat(bci_data.Data(:,6:7))+10;
bciunit = cell2mat(bci_data.Data(:,8));
nevfiles = bci_data.Data(:,9);

rot = (bcituningtype == 1 | bcituningtype == 3| bcituningtype == 4);
bci = (bciunit == 1) & rot;
nonbci = (bciunit == 0) & rot;

dof_bci = sum(bci)-1;
dof_nonbci = sum(nonbci)-1;

%1 = mc1 
%2 = bc
%3 = mc2 
%4 = dc 

mc1bc_difftheta = abs(center_angles(all_r2(rot,1), all_r2(rot,2)));
mc1mc2_difftheta = abs(center_angles(all_r2(rot,1), all_r2(rot,3)));
mc2bc_difftheta = abs(center_angles(all_r2(rot,3), all_r2(rot,2)));
mcdc_difftheta = abs(center_angles(all_r2(rot,1), all_r2(rot,4)));
bcdc_difftheta = abs(center_angles(all_r2(rot,2), all_r2(rot,4)));

bci_mc1bc_difftheta = abs(center_angles(all_r2(bci,1), all_r2(bci,2)));
bci_mc1mc2_difftheta = abs(center_angles(all_r2(bci,1), all_r2(bci,3)));
bci_mc2bc_difftheta = abs(center_angles(all_r2(bci,3), all_r2(bci,2)));
bci_mcdc_difftheta = abs(center_angles(all_r2(bci,1), all_r2(bci,4)));
bci_bcdc_difftheta = abs(center_angles(all_r2(bci,2), all_r2(bci,4)));

nonbci_mc1bc_difftheta = abs(center_angles(all_r2(nonbci,1), all_r2(nonbci,2)));
nonbci_mc1mc2_difftheta = abs(center_angles(all_r2(nonbci,1), all_r2(nonbci,3)));
nonbci_mc2bc_difftheta = abs(center_angles(all_r2(nonbci,3), all_r2(nonbci,2)));
nonbci_mcdc_difftheta = abs(center_angles(all_r2(nonbci,1), all_r2(nonbci,4)));
nonbci_bcdc_difftheta = abs(center_angles(all_r2(nonbci,2), all_r2(nonbci,4)));

mu_mc1bc_difftheta = mean(mc1bc_difftheta);
mu_mc1mc2_difftheta = mean(mc1mc2_difftheta);
mu_mc2bc_difftheta = mean(mc2bc_difftheta);
mu_mcdc_difftheta = mean(mcdc_difftheta);
mu_bcdc_difftheta = mean(bcdc_difftheta);

std_mc1bc_difftheta = std(mc1bc_difftheta);
std_mc1mc2_difftheta = std(mc1mc2_difftheta);
std_mc2bc_difftheta = std(mc2bc_difftheta);
std_mcdc_difftheta = std(mcdc_difftheta);
std_bcdc_difftheta = std(bcdc_difftheta);

mu_bci_mc1bc_difftheta = mean(bci_mc1bc_difftheta);
mu_bci_mc1mc2_difftheta = mean(bci_mc1mc2_difftheta);
mu_bci_mc2bc_difftheta = mean(bci_mc2bc_difftheta);
mu_bci_mcdc_difftheta = mean(bci_mcdc_difftheta);
mu_bci_bcdc_difftheta = mean(bci_bcdc_difftheta);

std_bci_mc1bc_difftheta = std(bci_mc1bc_difftheta);
std_bci_mc1mc2_difftheta = std(bci_mc1mc2_difftheta);
std_bci_mc2bc_difftheta = std(bci_mc2bc_difftheta);
std_bci_mcdc_difftheta = std(bci_mcdc_difftheta);
std_bci_bcdc_difftheta = std(bci_bcdc_difftheta);

mu_nonbci_mc1bc_difftheta  = mean(nonbci_mc1bc_difftheta);
mu_nonbci_mc1mc2_difftheta = mean(nonbci_mc1mc2_difftheta);
mu_nonbci_mc2bc_difftheta  = mean(nonbci_mc2bc_difftheta);
mu_nonbci_mcdc_difftheta = mean(nonbci_mcdc_difftheta);
mu_nonbci_bcdc_difftheta = mean(nonbci_bcdc_difftheta);

std_nonbci_mc1bc_difftheta = std(nonbci_mc1bc_difftheta);
std_nonbci_mc1mc2_difftheta = std(nonbci_mc1mc2_difftheta);
std_nonbci_mc2bc_difftheta = std(nonbci_mc2bc_difftheta);
std_nonbci_mcdc_difftheta = std(nonbci_mcdc_difftheta);
std_nonbci_bcdc_difftheta = std(nonbci_bcdc_difftheta);

[h_mc1bc, p_mc1bc] = ttest2(bci_mc1bc_difftheta, nonbci_mc1bc_difftheta, 0.05)
[h_mc1mc2, p_mc1mc2] = ttest2(bci_mc1mc2_difftheta, nonbci_mc1mc2_difftheta, 0.05)
[h_mc2bc, p_mc2bc] = ttest2(bci_mc2bc_difftheta, nonbci_mc2bc_difftheta, 0.05)
[h_mcdc, p_mcdc] = ttest2(bci_mcdc_difftheta, nonbci_mcdc_difftheta, 0.05)
[h_bcdc, p_bcdc] = ttest2(bci_bcdc_difftheta, nonbci_bcdc_difftheta, 0.05)

[h_mc1bc, p_mc1bc] = ttest2(mc1mc2_difftheta, mc1bc_difftheta, 0.05)
[h_mc2bc, p_mc2bc] = ttest2(mc1mc2_difftheta, mc2bc_difftheta, 0.05)
[h_mcdc, p_mcdc] = ttest2(mc1mc2_difftheta, mcdc_difftheta, 0.05)
[h_bcdc, p_bcdc] = ttest2(mc1mc2_difftheta, bcdc_difftheta, 0.05)

[h_bcdc2, p_bcdc2] = ttest2(mc1bc_difftheta, mcdc_difftheta, 0.05);

[h_mc1bc_bci, p_mc1bc_bci] = ttest2(bci_mc1mc2_difftheta, bci_mc1bc_difftheta, 0.05)
[h_mc2bc_bci, p_mc2bc_bci] = ttest2(bci_mc1mc2_difftheta, bci_mc2bc_difftheta, 0.05)
[h_mcdc_bci, p_mcdc_bci] = ttest2(bci_mc1mc2_difftheta, bci_mcdc_difftheta, 0.05)
[h_bcdc_bci, p_bcdc_bci] = ttest2(bci_mc1mc2_difftheta, bci_bcdc_difftheta, 0.05)

[h_mc1bc_nonbci, p_mc1bc_nonbci] = ttest2(nonbci_mc1mc2_difftheta, nonbci_mc1bc_difftheta, 0.05);
[h_mc2bc_nonbci, p_mc2bc_nonbci] = ttest2(nonbci_mc1mc2_difftheta, nonbci_mc2bc_difftheta, 0.05);
[h_mcdc_nonbci, p_mcdc_nonbci] = ttest2(nonbci_mc1mc2_difftheta, nonbci_mcdc_difftheta, 0.05);
[h_bcdc_nonbci, p_bcdc_nonbci] = ttest2(nonbci_mc1mc2_difftheta, nonbci_bcdc_difftheta, 0.05);

figure 
bar(180*[mean(mc1mc2_difftheta), mean(mc1bc_difftheta), mean(mcdc_difftheta)]/pi)
hold on 
errorbar(180*[mean(mc1mc2_difftheta), mean(mc1bc_difftheta), mean(mcdc_difftheta)]/pi, ...
	180*[std(mc1mc2_difftheta), std(mc1bc_difftheta), std(mcdc_difftheta)]/pi)
ylabel('|\Delta \theta|')
ylim([0 120])
saveplot(gcf, './figures/tuningangle-conditions-bargraph-rotated.eps')
