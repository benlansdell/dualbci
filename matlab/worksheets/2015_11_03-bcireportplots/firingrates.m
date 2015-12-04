conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gather data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_data = fetch(exec(conn, ['SELECT flin1.firingrate, flin2.firingrate, flin3.firingrate, flin4.firingrate, flin5.firingrate FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `units` flin1 '...
'ON flin1.`nev file` = et1.`manualrecording`'...
'INNER JOIN `units` flin2 '...
'ON flin2.`nev file` = et1.`1DBCrecording`'...
'INNER JOIN `units` flin3 '...
'ON flin3.`nev file` = et1.`manualrecordingafter`'...
'INNER JOIN `units` flin4 '...
'ON flin4.`nev file` = et1.`1DBCrecordingafter`'...
'INNER JOIN `units` flin5 '...
'ON flin5.`nev file` = et1.`dualrecording`'...
'WHERE flin1.unit = flin2.unit AND flin2.unit = flin3.unit AND flin2.unit = flin4.unit AND flin2.unit = flin5.unit']));
all_fr = cell2mat(all_data.Data(:,1:5));

corrs(1) = corr(all_fr(:,1), all_fr(:,2));
corrs(2) = corr(all_fr(:,1), all_fr(:,3));
corrs(3) = corr(all_fr(:,3), all_fr(:,2));
corrs(4) = corr(all_fr(:,2), all_fr(:,4));
corrs(5) = corr(all_fr(:,1), all_fr(:,5));
corrs(6) = corr(all_fr(:,2), all_fr(:,5));

clf
subplot(2,3,1)
plot(all_fr(:,1), all_fr(:,2), '.')
xlabel('R^2 MC1')
ylabel('R^2 BC1')
title(['corr: ' num2str(corrs(1))])
subplot(2,3,2)
plot(all_fr(:,1), all_fr(:,3), '.')
xlabel('R^2 MC1')
ylabel('R^2 MC2')
title(['corr: ' num2str(corrs(2))])
subplot(2,3,3)
plot(all_fr(:,3), all_fr(:,2), '.')
ylabel('R^2 BC1')
xlabel('R^2 MC2')
title(['corr: ' num2str(corrs(3))])
subplot(2,3,4)
plot(all_fr(:,2), all_fr(:,4), '.')
xlabel('R^2 BC1')
ylabel('R^2 BC2')
%xlim([0 3])
%ylim([0 8])
title(['corr: ' num2str(corrs(4))])
subplot(2,3,5)
plot(all_fr(:,1), all_fr(:,5), '.')
xlabel('R^2 MC1')
ylabel('R^2 DC')
title(['corr: ' num2str(corrs(5))])
subplot(2,3,6)
plot(all_fr(:,2), all_fr(:,5), '.')
xlabel('R^2 BC1')
ylabel('R^2 DC')
title(['corr: ' num2str(corrs(6))])
saveplot(gcf, './worksheets/2015_11_03-bcireportplots/firingratecorr.eps', 'eps', [10 6])