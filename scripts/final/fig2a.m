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


corrsdir(1) = corr(all_r2(rot,1), translate_angles(all_r2(rot,1), all_r2(rot,2)));
corrsdir(2) = corr(all_r2(rot,1), translate_angles(all_r2(rot,1), all_r2(rot,3)));
corrsdir(3) = corr(all_r2(rot,3), translate_angles(all_r2(rot,3), all_r2(rot,2)));
corrsdir(4) = corr(all_r2(rot,1), translate_angles(all_r2(rot,1), all_r2(rot,4)));
corrsdir(5) = corr(all_r2(rot,2), translate_angles(all_r2(rot,2), all_r2(rot,4)));

corrsbci(1) = corr(all_r2(bci,1), translate_angles(all_r2(bci,1), all_r2(bci,2)));
corrsbci(2) = corr(all_r2(bci,1), translate_angles(all_r2(bci,1), all_r2(bci,3)));
corrsbci(3) = corr(all_r2(bci,3), translate_angles(all_r2(bci,3), all_r2(bci,2)));
corrsbci(4) = corr(all_r2(bci,1), translate_angles(all_r2(bci,1), all_r2(bci,4)));
corrsbci(5) = corr(all_r2(bci,2), translate_angles(all_r2(bci,2), all_r2(bci,4)));

corrsnonbci(1) = corr(all_r2(nonbci,1), translate_angles(all_r2(nonbci,1), all_r2(nonbci,2)));
corrsnonbci(2) = corr(all_r2(nonbci,1), translate_angles(all_r2(nonbci,1), all_r2(nonbci,3)));
corrsnonbci(3) = corr(all_r2(nonbci,3), translate_angles(all_r2(nonbci,3), all_r2(nonbci,2)));
corrsnonbci(4) = corr(all_r2(nonbci,1), translate_angles(all_r2(nonbci,1), all_r2(nonbci,4)));
corrsnonbci(5) = corr(all_r2(nonbci,2), translate_angles(all_r2(nonbci,2), all_r2(nonbci,4)));

figure 
clf
cc = ones(size(bciunit));
cc(bciunit == 1) = 2;
colors = [1 0 0; 0 0 1];
c = [];
for idx = 1:size(cc, 1);
	c(idx,1:3) = colors(cc(idx),:);
end
colormap(colors)
subplot(2,3,1)
scatter(180/pi*all_r2(rot,1), 180/pi*translate_angles(all_r2(rot,1), all_r2(rot,2)), [], cc(rot), '.')
xlabel('\theta MC1')
ylabel('\theta BC1')
title(['corr bci: ' num2str(corrsbci(1)) ' corr nonbci: ' num2str(corrsnonbci(1))])
subplot(2,3,2)
scatter(180/pi*all_r2(rot,1), 180/pi*translate_angles(all_r2(rot,1), all_r2(rot,3)), [], cc(rot), '.')
xlabel('\theta MC1')
ylabel('\theta MC2')
ylim([-200 600])
title(['corr bci: ' num2str(corrsbci(2)) ' corr nonbci: ' num2str(corrsnonbci(2))])
subplot(2,3,3)
scatter(180/pi*all_r2(rot,3), 180/pi*translate_angles(all_r2(rot,3), all_r2(rot,2)), [], cc(rot), '.')
ylabel('\theta BC1')
xlabel('\theta MC2')
title(['corr bci: ' num2str(corrsbci(3)) ' corr nonbci: ' num2str(corrsnonbci(3))])
subplot(2,3,5)
scatter(180/pi*all_r2(rot,1), 180/pi*translate_angles(all_r2(rot,1), all_r2(rot,4)), [], cc(rot), '.')
xlabel('\theta MC1')
ylabel('\theta DC')
title(['corr bci: ' num2str(corrsbci(4)) ' corr nonbci: ' num2str(corrsnonbci(4))])
subplot(2,3,6)
scatter(180/pi*all_r2(rot,2), 180/pi*translate_angles(all_r2(rot,2), all_r2(rot,4)), [], cc(rot), '.')
xlabel('\theta BC1')
ylabel('\theta DC')
title(['corr bci: ' num2str(corrsbci(5)) ' corr nonbci: ' num2str(corrsnonbci(5))])
saveplot(gcf, './figures/tuningangle-bciVnonbci-rotated.eps', 'eps', [10 6])