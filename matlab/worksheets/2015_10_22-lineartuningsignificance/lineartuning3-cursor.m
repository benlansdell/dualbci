conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gather data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_data = fetch(exec(conn, ['SELECT fl1.size, fl2.size, fl3.size, fl4.size, fl1.dir, fl2.dir, fl3.dir, fl4.dir FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`manualrecording` '...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'INNER JOIN `fits` flin2 '...
'ON flin2.`nev file` = et1.`1DBCrecording` '...
'INNER JOIN `fits_linear` fl2 '...
'ON flin2.id = fl2.id '...
'INNER JOIN `fits` flin3 '...
'ON flin3.`nev file` = et1.`manualrecordingafter` '...
'INNER JOIN `fits_linear` fl3 '...
'ON flin3.id = fl3.id '...
'INNER JOIN `fits` flin4 '...
'ON flin4.`nev file` = et1.`1DBCrecordingafter` '...
'INNER JOIN `fits_linear` fl4 '...
'ON flin4.id = fl4.id '...
'WHERE flin1.modelID = 7 AND flin2.modelID = 7 AND flin3.modelID = 7 AND flin4.modelID = 7 ' ...
'AND flin1.unit = flin2.unit AND flin2.unit = flin3.unit AND flin2.unit = flin4.unit']));
all_d = cell2mat(all_data.Data(:,1:4));

MI(1) = mutualinfo(all_d(:,1), all_d(:,2));
MI(2) = mutualinfo(all_d(:,1), all_d(:,3));
MI(3) = mutualinfo(all_d(:,3), all_d(:,2));
MI(4) = mutualinfo(all_d(:,2), all_d(:,4));

clf
subplot(2,3,1)
plot(all_d(:,1), all_d(:,2), '.')
xlabel('Tuning strength MC1')
ylabel('Tuning strength BC1')
title(['Mutual information: ' num2str(MI(1))])
subplot(2,3,2)
plot(all_d(:,1), all_d(:,3), '.')
xlabel('Tuning strength MC1')
ylabel('Tuning strength MC2')
title(['Mutual information: ' num2str(MI(2))])
subplot(2,3,3)
plot(all_d(:,3), all_d(:,2), '.')
ylabel('Tuning strength BC1')
xlabel('Tuning strength MC2')
title(['Mutual information: ' num2str(MI(3))])
subplot(2,3,4)
plot(all_d(:,2), all_d(:,4), '.')
xlabel('Tuning strength BC1')
ylabel('Tuning strength BC2')
xlim([0 3])
ylim([0 8])
title(['Mutual information: ' num2str(MI(4))])
saveplot(gcf, './worksheets/2015_10_22-lineartuningsignificance/tuning-cursor-strengthMCBCMC2.eps')