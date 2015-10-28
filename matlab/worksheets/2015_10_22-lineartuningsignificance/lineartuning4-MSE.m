conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gather data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_data = fetch(exec(conn, ['SELECT flin1.`mse out`, flin2.`mse out`, flin3.`mse out`, flin4.`mse out`, flin5.`mse out` FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`manualrecording` '...
'INNER JOIN `fits` flin2 '...
'ON flin2.`nev file` = et1.`1DBCrecording` '...
'INNER JOIN `fits` flin3 '...
'ON flin3.`nev file` = et1.`manualrecordingafter` '...
'INNER JOIN `fits` flin4 '...
'ON flin4.`nev file` = et1.`1DBCrecordingafter` '...
'INNER JOIN `fits` flin5 '...
'ON flin5.`nev file` = et1.`dualrecording` '...
'WHERE flin1.modelID = 1 AND flin2.modelID = 1 AND flin3.modelID = 1 AND flin4.modelID = 1 AND flin5.modelID = 1 ' ...
'AND flin1.unit = flin2.unit AND flin2.unit = flin3.unit AND flin2.unit = flin4.unit AND flin2.unit = flin5.unit']));
all_d = cell2mat(all_data.Data(:,1:5));

MI(1) = mutualinfo(all_d(:,1), all_d(:,2));
MI(2) = mutualinfo(all_d(:,1), all_d(:,3));
MI(3) = mutualinfo(all_d(:,3), all_d(:,2));
MI(4) = mutualinfo(all_d(:,2), all_d(:,4));
MI(5) = mutualinfo(all_d(:,1), all_d(:,5));
MI(6) = mutualinfo(all_d(:,2), all_d(:,5));
MI(7) = mutualinfo(all_d(:,3), all_d(:,5));
MI(8) = mutualinfo(all_d(:,4), all_d(:,5));

clf
subplot(3,3,1)
plot(all_d(:,1), all_d(:,2), '.')
xlabel('Tuning MSE MC1')
ylabel('Tuning MSE BC1')
title(['Mutual information: ' num2str(MI(1))])
subplot(3,3,2)
plot(all_d(:,1), all_d(:,3), '.')
xlabel('Tuning MSE MC1')
ylabel('Tuning MSE MC2')
title(['Mutual information: ' num2str(MI(2))])
subplot(3,3,3)
plot(all_d(:,3), all_d(:,2), '.')
ylabel('Tuning MSE BC1')
xlabel('Tuning MSE MC2')
title(['Mutual information: ' num2str(MI(3))])
subplot(3,3,4)
plot(all_d(:,2), all_d(:,4), '.')
xlabel('Tuning MSE BC1')
ylabel('Tuning MSE BC2')
xlim([0 3])
ylim([0 8])
title(['Mutual information: ' num2str(MI(4))])
subplot(3,3,5)
plot(all_d(:,1), all_d(:,5), '.')
xlabel('Tuning MSE MC1')
ylabel('Tuning MSE DC')
title(['Mutual information: ' num2str(MI(5))])
subplot(3,3,6)
plot(all_d(:,2), all_d(:,5), '.')
xlabel('Tuning MSE BC1')
ylabel('Tuning MSE DC')
title(['Mutual information: ' num2str(MI(6))])
subplot(3,3,7)
plot(all_d(:,3), all_d(:,5), '.')
xlabel('Tuning MSE MC2')
ylabel('Tuning MSE DC')
title(['Mutual information: ' num2str(MI(7))])
subplot(3,3,8)
plot(all_d(:,4), all_d(:,5), '.')
xlabel('Tuning MSE BC2')
ylabel('Tuning MSE DC')
title(['Mutual information: ' num2str(MI(8))])
saveplot(gcf, './worksheets/2015_10_22-lineartuningsignificance/tuningMSEMCBCMC2.eps')


