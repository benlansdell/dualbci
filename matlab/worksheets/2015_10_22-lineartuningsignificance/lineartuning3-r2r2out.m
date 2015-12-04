conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gather data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_data = fetch(exec(conn, ['SELECT fl1.r2, fl2.r2, fl3.r2, fl4.r2, fl5.r2, fl1.r2out, fl2.r2out, fl3.r2out, fl4.r2out, fl5.r2out FROM '...
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
'INNER JOIN `fits` flin4 '...
'ON flin4.`nev file` = et1.`1DBCrecordingafter`'...
'INNER JOIN `fits_linear` fl4 '...
'ON flin4.id = fl4.id '...
'INNER JOIN `fits` flin5 '...
'ON flin5.`nev file` = et1.`dualrecording`'...
'INNER JOIN `fits_linear` fl5 '...
'ON flin5.id = fl5.id '...
'WHERE flin1.modelID = 1 AND flin2.modelID = 1 AND flin3.modelID = 1 AND flin4.modelID = 1 AND flin5.modelID = 1 ' ...
'AND flin1.unit = flin2.unit AND flin2.unit = flin3.unit AND flin2.unit = flin4.unit AND flin2.unit = flin5.unit']));
all_r2 = cell2mat(all_data.Data(:,1:5));
all_r2out = cell2mat(all_data.Data(:,6:10));

corrs(1) = corr(all_r2(:,1), all_r2(:,2));
corrs(2) = corr(all_r2(:,1), all_r2(:,3));
corrs(3) = corr(all_r2(:,3), all_r2(:,2));
%corrs(4) = corr(all_r2(:,2), all_r2(:,4));
corrs(5) = corr(all_r2(:,1), all_r2(:,5));
corrs(6) = corr(all_r2(:,2), all_r2(:,5));

MI(1) = mutualinfo(all_r2(:,1)*100, all_r2(:,2)*100);
MI(2) = mutualinfo(all_r2(:,1)*100, all_r2(:,3)*100);
MI(3) = mutualinfo(all_r2(:,3)*100, all_r2(:,2)*100);
%MI(4) = mutualinfo(all_r2(:,2)*100, all_r2(:,4)*100);
MI(5) = mutualinfo(all_r2(:,1)*100, all_r2(:,5)*100);
MI(6) = mutualinfo(all_r2(:,2)*100, all_r2(:,5)*100);

clf
subplot(2,3,1)
plot(all_r2(:,1), all_r2(:,2), '.')
xlabel('R^2 MC1')
ylabel('R^2 BC1')
title(['Mutual information: ' num2str(MI(1)) '. corr: ' num2str(corrs(1))])
subplot(2,3,2)
plot(all_r2(:,1), all_r2(:,3), '.')
xlabel('R^2 MC1')
ylabel('R^2 MC2')
title(['Mutual information: ' num2str(MI(2)) '. corr: ' num2str(corrs(2))])
subplot(2,3,3)
plot(all_r2(:,3), all_r2(:,2), '.')
ylabel('R^2 BC1')
xlabel('R^2 MC2')
title(['Mutual information: ' num2str(MI(3)) '. corr: ' num2str(corrs(3))])
subplot(2,3,4)
plot(all_r2(:,2), all_r2(:,4), '.')
xlabel('R^2 BC1')
ylabel('R^2 BC2')
%xlim([0 3])
%ylim([0 8])
title(['Mutual information: ' num2str(MI(4)) '. corr: ' num2str(corrs(4))])
subplot(2,3,5)
plot(all_r2(:,1), all_r2(:,5), '.')
xlabel('R^2 MC1')
ylabel('R^2 DC')
title(['Mutual information: ' num2str(MI(5)) '. corr: ' num2str(corrs(5))])
subplot(2,3,6)
plot(all_r2(:,2), all_r2(:,5), '.')
xlabel('R^2 BC1')
ylabel('R^2 DC')
title(['Mutual information: ' num2str(MI(6)) '. corr: ' num2str(corrs(6))])
saveplot(gcf, './worksheets/2015_10_22-lineartuningsignificance/tuningstrengthR2R2out.eps', 'eps', [10 6])