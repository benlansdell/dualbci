conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gather data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_data = fetch(exec(conn, ['SELECT flin1.dev, flin2.dev, flin3.dev, flin4.dev, flin5.dev FROM '...
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
'WHERE flin1.modelID = 31 AND flin2.modelID = 31 AND flin3.modelID = 31 AND flin4.modelID = 31 AND flin5.modelID = 31 ' ...
'AND flin1.unit = flin2.unit AND flin2.unit = flin3.unit AND flin2.unit = flin4.unit AND flin2.unit = flin5.unit']));
all_r2 = cell2mat(all_data.Data(:,1:5));

corrs(1) = corr(all_r2(:,1), all_r2(:,2));
corrs(2) = corr(all_r2(:,1), all_r2(:,3));
corrs(3) = corr(all_r2(:,3), all_r2(:,2));
corrs(4) = corr(all_r2(:,2), all_r2(:,4));
corrs(5) = corr(all_r2(:,1), all_r2(:,5));
corrs(6) = corr(all_r2(:,2), all_r2(:,5));

clf
subplot(2,3,1)
plot(all_r2(:,1), all_r2(:,2), '.')
xlabel('Dev MC1')
ylabel('Dev BC1')
title(['corr: ' num2str(corrs(1))])
subplot(2,3,2)
plot(all_r2(:,1), all_r2(:,3), '.')
xlabel('Dev MC1')
ylabel('Dev MC2')
title(['corr: ' num2str(corrs(2))])
subplot(2,3,3)
plot(all_r2(:,3), all_r2(:,2), '.')
ylabel('Dev BC1')
xlabel('Dev MC2')
title(['corr: ' num2str(corrs(3))])
subplot(2,3,4)
plot(all_r2(:,2), all_r2(:,4), '.')
xlabel('Dev BC1')
ylabel('Dev BC2')
%xlim([0 3])
%ylim([0 8])
title(['corr: ' num2str(corrs(4))])
subplot(2,3,5)
plot(all_r2(:,1), all_r2(:,5), '.')
xlabel('Dev MC1')
ylabel('Dev DC')
title(['corr: ' num2str(corrs(5))])
subplot(2,3,6)
plot(all_r2(:,2), all_r2(:,5), '.')
xlabel('Dev BC1')
ylabel('Dev DC')
title(['corr: ' num2str(corrs(6))])
saveplot(gcf, './worksheets/2015_11_03-bcireportplots/glmconsttuningdev.eps', 'eps', [10 6])