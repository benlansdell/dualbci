%conn = database('',databaseuser,databasepwd,'com.mysql.jdbc.Driver', ...
%	databaseurl);
%
%%Tuning angles, BCI units (velocity)
%bci_data = fetch(exec(conn, ['SELECT fl1.size, fl2.size, fl5.size, et1.`tuning_type`, rec1.`successrate`, rec2.`successrate`, '...
%' IF(EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et1.`1DBCrecording` AND bci.unit = flin1.unit), '...
%' 1, 0), ep1.value, ep2.value, ep5.value, fl3.size FROM '...
%'`experiment_tuning` et1 '...
%'INNER JOIN `fits` flin1 '...
%'ON flin1.`nev file` = et1.`manualrecording`'...
%'INNER JOIN `fits_linear` fl1 '...
%'ON flin1.id = fl1.id '...
%'INNER JOIN `estimates_parameters` ep1 '...
%'ON flin1.id = ep1.id AND ep1.label = "const" '...
%'INNER JOIN `fits` flin2 '...
%'ON flin2.`nev file` = et1.`1DBCrecording`'...
%'INNER JOIN `estimates_parameters` ep2 '...
%'ON flin2.id = ep2.id AND ep2.label = "const" '...
%'INNER JOIN `fits_linear` fl2 '...
%'ON flin2.id = fl2.id '...
%'INNER JOIN `fits` flin3 '...
%'ON flin3.`nev file` = et1.`manualrecordingafter`'...
%'INNER JOIN `estimates_parameters` ep3 '...
%'ON flin3.id = ep3.id AND ep3.label = "const" '...
%'INNER JOIN `fits_linear` fl3 '...
%'ON flin3.id = fl3.id '...
%'INNER JOIN `fits` flin5 '...
%'ON flin5.`nev file` = et1.`dualrecording`'...
%'INNER JOIN `estimates_parameters` ep5 '...
%'ON flin5.id = ep5.id AND ep5.label = "const" '...
%'INNER JOIN `fits_linear` fl5 '...
%'ON flin5.id = fl5.id '...
%'INNER JOIN `recordings` rec1 '...
%'ON rec1.`nev file` = et1.`1DBCrecording` '...
%'INNER JOIN `recordings` rec2 '...
%'ON rec2.`nev file` = et1.`dualrecording` '...
%'WHERE flin1.modelID = 30 AND flin2.modelID = 30 AND flin3.modelID = 30 AND flin5.modelID = 30 ' ...
%'AND flin1.unit = flin2.unit AND flin2.unit = flin3.unit AND flin2.unit = flin5.unit ' ...
%'AND fl1.r2 > .01 AND fl3.r2 > .01']));
%
%save('./scripts/figS3b.mat')
load('./scripts/figS3b.mat')

all_r2 = cell2mat(bci_data.Data(:,1:3));
bcituningtype = cell2mat(bci_data.Data(:,4));
bciperformance = 50*cell2mat(bci_data.Data(:,5:6))+10;
bciunit = cell2mat(bci_data.Data(:,7));
mc2size = cell2mat(bci_data.Data(:,end));

unrot = (bcituningtype == 5);
bci = (bciunit == 1) & unrot;
nonbci = (bciunit == 0) & unrot;

dof_bci = sum(bci)-1;
dof_nonbci = sum(nonbci)-1;

%1 = mc1 
%2 = bc
%3 = dc 

mc1bc_diffsize = abs(all_r2(unrot,1) - all_r2(unrot,2));
mcdc_diffsize = abs(all_r2(unrot,1)- all_r2(unrot,3));
bcdc_diffsize = abs(all_r2(unrot,2)- all_r2(unrot,3));
mc1mc2_diffsize = abs(all_r2(unrot,1)- mc2size(unrot,1));

bci_mc1bc_diffsize = abs(all_r2(bci,1)- all_r2(bci,2));
bci_mcdc_diffsize = abs(all_r2(bci,1)- all_r2(bci,3));
bci_bcdc_diffsize = abs(all_r2(bci,2)- all_r2(bci,3));

nonbci_mc1bc_diffsize = abs(all_r2(nonbci,1)- all_r2(nonbci,2));
nonbci_mcdc_diffsize = abs(all_r2(nonbci,1)- all_r2(nonbci,3));
nonbci_bcdc_diffsize = abs(all_r2(nonbci,2)- all_r2(nonbci,3));

mu_bci_mc1bc_diffsize = mean(bci_mc1bc_diffsize);
mu_bci_mcdc_diffsize = mean(bci_mcdc_diffsize);
mu_bci_bcdc_diffsize = mean(bci_bcdc_diffsize);
mu_nonbci_mc1bc_diffsize  = mean(nonbci_mc1bc_diffsize);
mu_nonbci_mcdc_diffsize = mean(nonbci_mcdc_diffsize);
mu_nonbci_bcdc_diffsize = mean(nonbci_bcdc_diffsize);

std_bci_mc1bc_diffsize = std(bci_mc1bc_diffsize);
std_bci_mcdc_diffsize = std(bci_mcdc_diffsize);
std_bci_bcdc_diffsize = std(bci_bcdc_diffsize);
std_nonbci_mc1bc_diffsize  = std(nonbci_mc1bc_diffsize);
std_nonbci_mcdc_diffsize = std(nonbci_mcdc_diffsize);
std_nonbci_bcdc_diffsize = std(nonbci_bcdc_diffsize);

[h_mc1bc, p_mc1bc] = ttest2(bci_mc1bc_diffsize, nonbci_mc1bc_diffsize, 0.05)
[h_mcdc, p_mcdc] = ttest2(bci_mcdc_diffsize, nonbci_mcdc_diffsize, 0.05)
[h_bcdc, p_bcdc] = ttest2(bci_bcdc_diffsize, nonbci_bcdc_diffsize, 0.05)

[h_mc1bc, p_mc1bc] = ttest2(mc1mc2_diffsize, mc1bc_diffsize, 0.05)
[h_mcdc, p_mcdc] = ttest2(mc1mc2_diffsize, mcdc_diffsize, 0.05)
%[h_bcdc, p_bcdc] = ttest2(mc1mc2_diffsize, bcdc_diffsize, 0.05)

[h_bcdc2, p_bcdc2] = ttest2(mc1bc_diffsize, mcdc_diffsize, 0.05)

figure 
bar([mu_bci_mc1bc_diffsize, mu_nonbci_mc1bc_diffsize, mu_bci_mcdc_diffsize, mu_nonbci_mcdc_diffsize])
hold on 
errorbar([mu_bci_mc1bc_diffsize, mu_nonbci_mc1bc_diffsize, mu_bci_mcdc_diffsize, mu_nonbci_mcdc_diffsize], [std_bci_mc1bc_diffsize, std_nonbci_mc1bc_diffsize, std_bci_mcdc_diffsize, std_nonbci_mcdc_diffsize])
saveplot(gcf, './figures/tuningstrengthBCI-nonBCI-bargraph-unrotated.eps')

%correlations
corrsbci(1) = corr(all_r2(bci,1), all_r2(bci,2));
corrsbci(2) = corr(all_r2(bci,1), all_r2(bci,3));
corrsbci(3) = corr(all_r2(bci,2), all_r2(bci,3));

corrsnonbci(1) = corr(all_r2(nonbci,1),all_r2(nonbci,2));
corrsnonbci(2) = corr(all_r2(nonbci,1),all_r2(nonbci,3));
corrsnonbci(3) = corr(all_r2(nonbci,2),all_r2(nonbci,3));

%[h_mc1bc_bci, p_mc1bc_bci] = ttest2(bci_mc1mc2_difftheta, bci_mc1bc_difftheta, 0.05)
%[h_mcdc_bci, p_mcdc_bci] = ttest2(bci_mc1mc2_difftheta, bci_mcdc_difftheta, 0.05)
%[h_bcdc_bci, p_bcdc_bci] = ttest2(bci_mc1mc2_difftheta, bci_bcdc_difftheta, 0.05)

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
scatter(all_r2(unrot,1), all_r2(unrot,2), [], cc)
xlabel('strength MC1')
ylabel('strength BC1')
title(['corr bci: ' num2str(corrsbci(1)) ' corr nonbci: ' num2str(corrsnonbci(1))])
subplot(2,2,2)
scatter(all_r2(unrot,1), all_r2(unrot,3), [], cc)
xlabel('strength MC1')
ylabel('strength DC')
title(['corr bci: ' num2str(corrsbci(2)) ' corr nonbci: ' num2str(corrsnonbci(2))])
subplot(2,2,3)
scatter(all_r2(unrot,2), all_r2(unrot,3), [], cc)
xlabel('strength BC1')
ylabel('strength DC')
title(['corr bci: ' num2str(corrsbci(3)) ' corr nonbci: ' num2str(corrsnonbci(3))])
saveplot(gcf, './figures/tuningstrengthBCI-nonBCI-scatter-unrotated.eps')