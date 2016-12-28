conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
grangerMCBC = fetch(exec(conn, ['SELECT fBC.analyses_id, fMC.`nev file`, fMC.unit tounit, '...
' egMC.fromunit, egMC.score egMCscore, egBC.score egBCscore FROM `experiment_tuning` et '...
'INNER JOIN `fits` fMC '...
'ON et.`manualrecording` = fMC.`nev file` '...
'INNER JOIN `estimates_granger` egMC '...
'ON egMC.`id` = fMC.`id` '...
'INNER JOIN `fits` fBC '...
'ON et.`1DBCrecording` = fBC.`nev file` '...
'INNER JOIN `estimates_granger` egBC '...
'ON egBC.`id` = fBC.`id` '...
'WHERE fMC.unit = fBC.unit AND egMC.fromunit = egBC.fromunit AND fMC.modelID = 29 AND fBC.modelID = 29 AND fMC.analyses_id = fBC.analyses_id']));
grangerMC = cell2mat(grangerMCBC.Data(:,5));
grangerBC = cell2mat(grangerMCBC.Data(:,6));

grangerMCMC2 = fetch(exec(conn, ['SELECT fBC.analyses_id, fMC.`nev file`, fMC.unit tounit, '...
' egMC.fromunit, egMC.score egMCscore, egBC.score egBCscore FROM `experiment_tuning` et '...
'INNER JOIN `fits` fMC '...
'ON et.`manualrecording` = fMC.`nev file` '...
'INNER JOIN `estimates_granger` egMC '...
'ON egMC.`id` = fMC.`id` '...
'INNER JOIN `fits` fBC '...
'ON et.`manualrecordingafter` = fBC.`nev file` '...
'INNER JOIN `estimates_granger` egBC '...
'ON egBC.`id` = fBC.`id` '...
'WHERE fMC.unit = fBC.unit AND egMC.fromunit = egBC.fromunit AND fMC.modelID = 29 AND fBC.modelID = 29 AND fMC.analyses_id = fBC.analyses_id']));
grangerMCb = cell2mat(grangerMCMC2.Data(:,5));
grangerMC2 = cell2mat(grangerMCMC2.Data(:,6));

grangerMCDC = fetch(exec(conn, ['SELECT fBC.analyses_id, fMC.`nev file`, fMC.unit tounit, '...
' egMC.fromunit, egMC.score egMCscore, egBC.score egBCscore FROM `experiment_tuning` et '...
'INNER JOIN `fits` fMC '...
'ON et.`manualrecording` = fMC.`nev file` '...
'INNER JOIN `estimates_granger` egMC '...
'ON egMC.`id` = fMC.`id` '...
'INNER JOIN `fits` fBC '...
'ON et.`dualrecording` = fBC.`nev file` '...
'INNER JOIN `estimates_granger` egBC '...
'ON egBC.`id` = fBC.`id` '...
'WHERE fMC.unit = fBC.unit AND egMC.fromunit = egBC.fromunit AND fMC.modelID = 29 AND fBC.modelID = 29 AND fMC.analyses_id = fBC.analyses_id']));
grangerMCc = cell2mat(grangerMCDC.Data(:,5));
grangerDC = cell2mat(grangerMCDC.Data(:,6));

grangerBCDC = fetch(exec(conn, ['SELECT fBC.analyses_id, fMC.`nev file`, fMC.unit tounit, '...
' egMC.fromunit, egMC.score egMCscore, egBC.score egBCscore FROM `experiment_tuning` et '...
'INNER JOIN `fits` fMC '...
'ON et.`1DBCrecording` = fMC.`nev file` '...
'INNER JOIN `estimates_granger` egMC '...
'ON egMC.`id` = fMC.`id` '...
'INNER JOIN `fits` fBC '...
'ON et.`dualrecording` = fBC.`nev file` '...
'INNER JOIN `estimates_granger` egBC '...
'ON egBC.`id` = fBC.`id` '...
'WHERE fMC.unit = fBC.unit AND egMC.fromunit = egBC.fromunit AND fMC.modelID = 29 AND fBC.modelID = 29 AND fMC.analyses_id = fBC.analyses_id']));
grangerBCb = cell2mat(grangerBCDC.Data(:,5));
grangerDC = cell2mat(grangerBCDC.Data(:,6));

corr(grangerMC, grangerBC)
corr(grangerMCb, grangerMC2)
corr(grangerMCc, grangerDC)
corr(grangerBCb, grangerDC)