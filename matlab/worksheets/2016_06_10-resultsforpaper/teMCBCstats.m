conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
teMCBC = fetch(exec(conn, ['SELECT fBC.analyses_id, fMC.`nev file`, fMC.unit tounit, '...
' egMC.fromunit, egMC.score egMCscore, egBC.score egBCscore, EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et.`dualrecording` AND bci.unit = fMC.unit), et.`tuning_type` FROM `experiment_tuning` et '...
'INNER JOIN `fits` fMC '...
'ON et.`manualrecording` = fMC.`nev file` '...
'INNER JOIN `estimates_te` egMC '...
'ON egMC.`id` = fMC.`id` '...
'INNER JOIN `fits` fBC '...
'ON et.`1DBCrecording` = fBC.`nev file` '...
'INNER JOIN `estimates_te` egBC '...
'ON egBC.`id` = fBC.`id` '...
'WHERE fMC.unit = fBC.unit AND egMC.fromunit = egBC.fromunit AND fMC.modelID = 37 AND fBC.modelID = 37 AND fMC.analyses_id = fBC.analyses_id']));
teMC = cell2mat(teMCBC.Data(:,5));
teBC = cell2mat(teMCBC.Data(:,6));
bcunitMCBC = cell2mat(teMCBC.Data(:,7));
tuningtype = cell2mat(teMCBC.Data(:,8));
rotMCBC = (tuningtype == 1 | tuningtype == 3| tuningtype == 4);

teMCMC2 = fetch(exec(conn, ['SELECT fBC.analyses_id, fMC.`nev file`, fMC.unit tounit, '...
' egMC.fromunit, egMC.score egMCscore, egBC.score egBCscore, EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et.`dualrecording` AND bci.unit = fMC.unit), et.`tuning_type` FROM `experiment_tuning` et '...
'INNER JOIN `fits` fMC '...
'ON et.`manualrecording` = fMC.`nev file` '...
'INNER JOIN `estimates_te` egMC '...
'ON egMC.`id` = fMC.`id` '...
'INNER JOIN `fits` fBC '...
'ON et.`manualrecordingafter` = fBC.`nev file` '...
'INNER JOIN `estimates_te` egBC '...
'ON egBC.`id` = fBC.`id` '...
'WHERE fMC.unit = fBC.unit AND egMC.fromunit = egBC.fromunit AND fMC.modelID = 37 AND fBC.modelID = 37 AND fMC.analyses_id = fBC.analyses_id']));
teMCb = cell2mat(teMCMC2.Data(:,5));
teMC2 = cell2mat(teMCMC2.Data(:,6));
bcunitMCMC2 = cell2mat(teMCMC2.Data(:,7));
tuningtype = cell2mat(teMCMC2.Data(:,8));
rotMCMC2 = (tuningtype == 1 | tuningtype == 3| tuningtype == 4);

teMCDC = fetch(exec(conn, ['SELECT fBC.analyses_id, fMC.`nev file`, fMC.unit tounit, '...
' egMC.fromunit, egMC.score egMCscore, egBC.score egBCscore,EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et.`dualrecording` AND bci.unit = fMC.unit), et.`tuning_type` FROM `experiment_tuning` et '...
'INNER JOIN `fits` fMC '...
'ON et.`manualrecording` = fMC.`nev file` '...
'INNER JOIN `estimates_te` egMC '...
'ON egMC.`id` = fMC.`id` '...
'INNER JOIN `fits` fBC '...
'ON et.`dualrecording` = fBC.`nev file` '...
'INNER JOIN `estimates_te` egBC '...
'ON egBC.`id` = fBC.`id` '...
'WHERE fMC.unit = fBC.unit AND egMC.fromunit = egBC.fromunit AND fMC.modelID = 37 AND fBC.modelID = 37 AND fMC.analyses_id = fBC.analyses_id']));
teMCc = cell2mat(teMCDC.Data(:,5));
teDC = cell2mat(teMCDC.Data(:,6));
bcunitMCDC = cell2mat(teMCDC.Data(:,7));
tuningtype = cell2mat(teMCDC.Data(:,8));
rotMCDC = (tuningtype == 1 | tuningtype == 3| tuningtype == 4);

teBCDC = fetch(exec(conn, ['SELECT fBC.analyses_id, fMC.`nev file`, fMC.unit tounit, '...
' egMC.fromunit, egMC.score egMCscore, egBC.score egBCscore,EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et.`dualrecording` AND bci.unit = fMC.unit), et.`tuning_type` FROM `experiment_tuning` et '...
'INNER JOIN `fits` fMC '...
'ON et.`1DBCrecording` = fMC.`nev file` '...
'INNER JOIN `estimates_te` egMC '...
'ON egMC.`id` = fMC.`id` '...
'INNER JOIN `fits` fBC '...
'ON et.`dualrecording` = fBC.`nev file` '...
'INNER JOIN `estimates_te` egBC '...
'ON egBC.`id` = fBC.`id` '...
'WHERE fMC.unit = fBC.unit AND egMC.fromunit = egBC.fromunit AND fMC.modelID = 37 AND fBC.modelID = 37 AND fMC.analyses_id = fBC.analyses_id']));
teBCb = cell2mat(teBCDC.Data(:,5));
teDC = cell2mat(teBCDC.Data(:,6));
bcunitBCDC = cell2mat(teBCDC.Data(:,7));
tuningtype = cell2mat(teBCDC.Data(:,8));
rotBCDC = (tuningtype == 1 | tuningtype == 3| tuningtype == 4);

corr(teMC, teBC)
corr(teMCb, teMC2)
corr(teMCc, teDC)
corr(teBCb, teDC)

colors = [1 0 0; 0 0 1];
colormap(colors)
subplot(2,2,1)
title('MC-BC')
cc = ones(size(bcunitMCBC));
cc(bcunitMCBC==1) = 2;
scatter(teMC, teBC, [], cc)
xlabel('MC')
ylabel('BC')
subplot(2,2,2)
title('MC-MC2')
cc = ones(size(bcunitMCMC2));
cc(bcunitMCMC2==1) = 2;
scatter(teMCb, teMC2, [], cc)
xlabel('MC')
ylabel('MC2')
subplot(2,2,3)
title('MC-DC')
cc = ones(size(bcunitMCDC));
cc(bcunitMCDC==1) = 2;
scatter(teMCc, teDC, [], cc)
xlabel('MC')
ylabel('DC')
subplot(2,2,4)
title('BC-DC')
cc = ones(size(bcunitBCDC));
cc(bcunitBCDC==1) = 2;
scatter(teBCb, teDC, [], cc)
xlabel('BC')
ylabel('DC')

%Save data for use in python
save('./worksheets/2016_06_10-resultsforpaper/teMCBCstats.mat')

diffteMCBC = abs(teMC-teBC);
diffteMCMC2 = abs(teMCb-teMC2);
diffteMCDC = abs(teMCc-teDC);
diffteBCDC = abs(teBCb-teDC);

figure
bar([mean(diffteMCMC2), mean(diffteMCBC), mean(diffteMCDC)]);

%Split into BC and non-BC units
diffteMCBCbci = abs(teMC(bcunitMCBC==1)-teBC(bcunitMCBC==1));
diffteMCBCnonbci = abs(teMC(bcunitMCBC~=1)-teBC(bcunitMCBC~=1));
diffteMCDCbci = abs(teMCc(bcunitMCDC==1)-teDC(bcunitMCDC==1));
diffteMCDCnonbci = abs(teMCc(bcunitMCDC~=1)-teDC(bcunitMCDC~=1));

figure
bar([mean(diffteMCBCbci), mean(diffteMCBCnonbci), mean(diffteMCDCbci), mean(diffteMCDCnonbci)]);
hold on 
errorbar([mean(diffteMCBCbci), mean(diffteMCBCnonbci), mean(diffteMCDCbci), mean(diffteMCDCnonbci)],...
	[std(diffteMCBCbci), std(diffteMCBCnonbci), std(diffteMCDCbci), std(diffteMCDCnonbci)]);

[h, p]=ttest2(diffteMCBCbci, diffteMCBCnonbci)
[h, p]=ttest2(diffteMCDCbci, diffteMCDCnonbci)

%Std dev. of change in TE between MC and MC2
std_diffteMCMC2 = std(teMCb-teMC2);
%== std_diffteMCMC2 =
%   5.4491e-04
