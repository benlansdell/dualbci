%conn = database('',databaseuser,databasepwd,'com.mysql.jdbc.Driver', ...
%	databaseurl);
%
%%List of files
%files = fetch(exec(conn, ['SELECT et.`manualrecording` FROM experiment_tuning et']));
%files = files.Data;
%
%%Torque tuning angle (velocity)
%deltaBCI = [];
%deltacotuned = [];
%deltaother = [];
%
%dother = [];
%dcotuned = [];
%
%nR = 50;
%
%rng(15);
%
%for rep = 1:nR
%	display(['rep = ' num2str(rep)])
%	for idx = 1:length(files)
%		cotunedgranger = [];
%		othergranger = [];
%		mcfile = files{idx};
%		all_data = fetch(exec(conn, ['SELECT fl1.dir, fl2.dir, fl5.dir, et1.`tuning_type`, flin1.unit FROM '...
%		'`experiment_tuning` et1 '...
%		'INNER JOIN `fits` flin1 '...
%		'ON flin1.`nev file` = et1.`manualrecording`'...
%		'INNER JOIN `fits_linear` fl1 '...
%		'ON flin1.id = fl1.id '...
%		'INNER JOIN `fits` flin2 '...
%		'ON flin2.`nev file` = et1.`1DBCrecording`'...
%		'INNER JOIN `fits_linear` fl2 '...
%		'ON flin2.id = fl2.id '...
%		'INNER JOIN `fits` flin5 '...
%		'ON flin5.`nev file` = et1.`dualrecording` '...
%		'INNER JOIN `fits_linear` fl5 '...
%		'ON flin5.id = fl5.id '...
%		'WHERE flin1.modelID = 30 AND flin2.modelID = 30 AND flin5.modelID = 30 ' ...
%		'AND flin1.unit = flin2.unit AND flin2.unit = flin5.unit ' ...
%		'AND fl1.r2 > 0.0 '...
%		'AND NOT EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et1.`1DBCrecording` AND bci.unit = flin1.unit) '...
%		'AND et1.`manualrecording` = "' mcfile '" '...
%		'AND (et1.`tuning_type` = 1 OR et1.`tuning_type` = 3 OR et1.`tuning_type` = 4)']));
%		if strcmp(all_data.Data, 'No Data')
%			continue 
%		end
%		if size(all_data.Data, 1) < 4 
%			continue 
%		end
%	
%		all_theta = cell2mat(all_data.Data(:,1:3));
%		tuningtype = cell2mat(all_data.Data(:,4));
%		otherunits = all_data.Data(:,5);
%		
%		nU = size(otherunits, 1);
%		nonbciidx = randsample(nU,1);
%		otheridx = setdiff(1:nU, nonbciidx);
%
%		bci_theta = all_theta(nonbciidx,:);
%		bci_unit = otherunits(nonbciidx);
%
%		all_theta = all_theta(otheridx,:);
%		tuningtype = tuningtype(otheridx);
%		otherunits = otherunits(otheridx);
%
%		nU = nU - 1;
%		
%		%Pick cotuned units to BCI units in MC
%		diff_MC_theta = cos(bci_theta(1) - all_theta(:,1));
%		
%		%Pick top two
%		[l, cotunedidx] = sort(diff_MC_theta, 'descend'); 
%		cotunedidx = cotunedidx(1:2);
%		bcicotuned = diff_MC_theta(cotunedidx);
%		cotunedunits = otherunits(cotunedidx);
%	
%		%Pick random two 
%		otheridx = randsample(setdiff(1:nU, cotunedidx), 2);
%		otherunitnames = otherunits(otheridx);
%	
%		%Granger for cotuned units
%		skipthisone = 0;
%		for j = 1:length(cotunedunits)
%			unit = cotunedunits{j};
%			GC = fetch(exec(conn, ['SELECT fBC.analyses_id, fMC.`nev file`, fMC.unit tounit, '...
%			' egMC.fromunit, egMC.score egMCscore, egBC.score egBCscore, fconstMC.dev, fconstBC.dev, et.`tuning_type` FROM `experiment_tuning` et '...
%			'INNER JOIN `fits` fMC '...
%			'ON et.`manualrecording` = fMC.`nev file` '...
%			'INNER JOIN `estimates_te` egMC '...
%			'ON egMC.`id` = fMC.`id` '...
%			'INNER JOIN `fits` fBC '...
%			'ON et.`1DBCrecording` = fBC.`nev file` '...
%			'INNER JOIN `estimates_te` egBC '...
%			'ON egBC.`id` = fBC.`id` '...
%			'INNER JOIN `fits` fconstBC '...
%			'ON fconstBC.`nev file` = et.`1DBCrecording` AND fconstBC.unit = fMC.unit '...
%			'INNER JOIN `fits` fconstMC '...
%			'ON fconstMC.`nev file` = et.`manualrecording` AND fconstMC.unit = fMC.unit '...
%			'WHERE fMC.unit = fBC.unit AND egMC.fromunit = egBC.fromunit AND fMC.modelID = 37 AND fBC.modelID = 37 AND fMC.analyses_id = fBC.analyses_id '...
%			'AND fconstBC.modelID = 31 AND fconstMC.modelID = 31 '...
%			'AND NOT EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et.`1DBCrecording` AND bci.unit = fMC.unit) '...
%			'AND egMC.fromunit = "' unit '" '...
%			'AND fMC.unit = "' bci_unit{1} '" '...
%			'AND et.`manualrecording` = "' mcfile '"']));
%			if strcmp(GC.Data, 'No Data')
%				skipthisone = 1;
%				continue;
%			end		
%			cotunedgranger(j,1:2) = cell2mat(GC.Data(:,5:6));
%		end
%		if skipthisone == 1
%			continue;
%		end
%	
%		%Granger for other units
%		skipthisone = 0;
%		for j = 1:length(otherunitnames)
%			unit = otherunitnames{j};
%			GC = fetch(exec(conn, ['SELECT fBC.analyses_id, fMC.`nev file`, fMC.unit tounit, '...
%			' egMC.fromunit, egMC.score egMCscore, egBC.score egBCscore, fconstMC.dev, fconstBC.dev, et.`tuning_type` FROM `experiment_tuning` et '...
%			'INNER JOIN `fits` fMC '...
%			'ON et.`manualrecording` = fMC.`nev file` '...
%			'INNER JOIN `estimates_te` egMC '...
%			'ON egMC.`id` = fMC.`id` AND egMC.`fromunit` = "' unit '" ' ...
%			'INNER JOIN `fits` fBC '...
%			'ON et.`1DBCrecording` = fBC.`nev file` '...
%			'INNER JOIN `estimates_te` egBC '...
%			'ON egBC.`id` = fBC.`id` '...
%			'INNER JOIN `fits` fconstBC '...
%			'ON fconstBC.`nev file` = et.`1DBCrecording` AND fconstBC.unit = fMC.unit '...
%			'INNER JOIN `fits` fconstMC '...
%			'ON fconstMC.`nev file` = et.`manualrecording` AND fconstMC.unit = fMC.unit '...
%			'WHERE fMC.unit = fBC.unit AND egMC.fromunit = egBC.fromunit AND fMC.modelID = 37 AND fBC.modelID = 37 AND fMC.analyses_id = fBC.analyses_id '...
%			'AND fconstBC.modelID = 31 AND fconstMC.modelID = 31 '...
%			'AND NOT EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et.`1DBCrecording` AND bci.unit = fMC.unit) '...
%			'AND fMC.unit = "' bci_unit{1} '" '...
%			'AND et.`manualrecording` = "' mcfile '"']));
%			if strcmp(GC.Data, 'No Data')
%				skipthisone = 1;
%				continue;
%			end		
%			othergranger(j,1:2) = cell2mat(GC.Data(:,5:6));
%		end
%		if skipthisone == 1
%			continue;
%		end
%	
%		diffcotuned = diff(cotunedgranger, 1,2);
%		diffother = diff(othergranger, 1,2);
%		dcotuned = [dcotuned(:); diffcotuned(:)];
%		dother = [dother(:); diffother(:)];
%	end
%end
%
%dcotunedMCBC = dcotuned;
%dotherMCBC = dother; 
%
%%%%%%%%%%%%%%%%%%%%%
%%Do the same for DC%
%%%%%%%%%%%%%%%%%%%%%
%
%dother = [];
%dcotuned = [];
%
%for rep = 1:nR
%	display(['rep = ' num2str(rep)])
%	for idx = 1:length(files)
%		cotunedgranger = [];
%		othergranger = [];
%		mcfile = files{idx};
%		all_data = fetch(exec(conn, ['SELECT fl1.dir, fl2.dir, fl5.dir, et1.`tuning_type`, flin1.unit FROM '...
%		'`experiment_tuning` et1 '...
%		'INNER JOIN `fits` flin1 '...
%		'ON flin1.`nev file` = et1.`manualrecording`'...
%		'INNER JOIN `fits_linear` fl1 '...
%		'ON flin1.id = fl1.id '...
%		'INNER JOIN `fits` flin2 '...
%		'ON flin2.`nev file` = et1.`1DBCrecording`'...
%		'INNER JOIN `fits_linear` fl2 '...
%		'ON flin2.id = fl2.id '...
%		'INNER JOIN `fits` flin5 '...
%		'ON flin5.`nev file` = et1.`dualrecording` '...
%		'INNER JOIN `fits_linear` fl5 '...
%		'ON flin5.id = fl5.id '...
%		'WHERE flin1.modelID = 30 AND flin2.modelID = 30 AND flin5.modelID = 30 ' ...
%		'AND flin1.unit = flin2.unit AND flin2.unit = flin5.unit ' ...
%		'AND fl1.r2 > 0.0 '...
%		'AND NOT EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et1.`1DBCrecording` AND bci.unit = flin1.unit) '...
%		'AND et1.`manualrecording` = "' mcfile '" '...
%		'AND (et1.`tuning_type` = 1 OR et1.`tuning_type` = 3 OR et1.`tuning_type` = 4)']));
%		if strcmp(all_data.Data, 'No Data')
%			continue 
%		end
%		if size(all_data.Data, 1) < 4 
%			continue 
%		end
%	
%		all_theta = cell2mat(all_data.Data(:,1:3));
%		tuningtype = cell2mat(all_data.Data(:,4));
%		otherunits = all_data.Data(:,5);
%		
%		nU = size(otherunits, 1);
%		nonbciidx = randsample(nU,1);
%		otheridx = setdiff(1:nU, nonbciidx);
%
%		bci_theta = all_theta(nonbciidx,:);
%		bci_unit = otherunits(nonbciidx);
%
%		all_theta = all_theta(otheridx,:);
%		tuningtype = tuningtype(otheridx);
%		otherunits = otherunits(otheridx);
%
%		nU = nU - 1;
%		
%		%Pick cotuned units to BCI units in MC
%		diff_MC_theta = cos(bci_theta(1) - all_theta(:,1));
%		
%		%Pick top two
%		[l, cotunedidx] = sort(diff_MC_theta, 'descend'); 
%		cotunedidx = cotunedidx(1:2);
%		bcicotuned = diff_MC_theta(cotunedidx);
%		cotunedunits = otherunits(cotunedidx);
%	
%		%Pick random two 
%		otheridx = randsample(setdiff(1:nU, cotunedidx), 2);
%		otherunitnames = otherunits(otheridx);
%	
%		%Granger for cotuned units
%		skipthisone = 0;
%		for j = 1:length(cotunedunits)
%			unit = cotunedunits{j};
%			GC = fetch(exec(conn, ['SELECT fDC.analyses_id, fMC.`nev file`, fMC.unit tounit, '...
%			' egMC.fromunit, egMC.score egMCscore, egDC.score egDCscore, fconstMC.dev, fconstDC.dev, et.`tuning_type` FROM `experiment_tuning` et '...
%			'INNER JOIN `fits` fMC '...
%			'ON et.`manualrecording` = fMC.`nev file` '...
%			'INNER JOIN `estimates_te` egMC '...
%			'ON egMC.`id` = fMC.`id` '...
%			'INNER JOIN `fits` fDC '...
%			'ON et.`dualrecording` = fDC.`nev file` '...
%			'INNER JOIN `estimates_te` egDC '...
%			'ON egDC.`id` = fDC.`id` '...
%			'INNER JOIN `fits` fconstDC '...
%			'ON fconstDC.`nev file` = et.`dualrecording` AND fconstDC.unit = fMC.unit '...
%			'INNER JOIN `fits` fconstMC '...
%			'ON fconstMC.`nev file` = et.`manualrecording` AND fconstMC.unit = fMC.unit '...
%			'WHERE fMC.unit = fDC.unit AND egMC.fromunit = egDC.fromunit AND fMC.modelID = 37 AND fDC.modelID = 37 AND fMC.analyses_id = fDC.analyses_id '...
%			'AND fconstDC.modelID = 31 AND fconstMC.modelID = 31 '...
%			'AND NOT EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et.`dualrecording` AND bci.unit = fMC.unit) '...
%			'AND egMC.fromunit = "' unit '" '...
%			'AND fMC.unit = "' bci_unit{1} '" '...
%			'AND et.`manualrecording` = "' mcfile '"']));
%			if strcmp(GC.Data, 'No Data')
%				skipthisone = 1;
%				continue;
%			end		
%			cotunedgranger(j,1:2) = cell2mat(GC.Data(:,5:6));
%		end
%		if skipthisone == 1
%			continue;
%		end
%	
%		%Granger for other units
%		skipthisone = 0;
%		for j = 1:length(otherunitnames)
%			unit = otherunitnames{j};
%			GC = fetch(exec(conn, ['SELECT fBC.analyses_id, fMC.`nev file`, fMC.unit tounit, '...
%			' egMC.fromunit, egMC.score egMCscore, egBC.score egBCscore, fconstMC.dev, fconstBC.dev, et.`tuning_type` FROM `experiment_tuning` et '...
%			'INNER JOIN `fits` fMC '...
%			'ON et.`manualrecording` = fMC.`nev file` '...
%			'INNER JOIN `estimates_te` egMC '...
%			'ON egMC.`id` = fMC.`id` AND egMC.`fromunit` = "' unit '" ' ...
%			'INNER JOIN `fits` fBC '...
%			'ON et.`1DBCrecording` = fBC.`nev file` '...
%			'INNER JOIN `estimates_te` egBC '...
%			'ON egBC.`id` = fBC.`id` '...
%			'INNER JOIN `fits` fconstBC '...
%			'ON fconstBC.`nev file` = et.`1DBCrecording` AND fconstBC.unit = fMC.unit '...
%			'INNER JOIN `fits` fconstMC '...
%			'ON fconstMC.`nev file` = et.`manualrecording` AND fconstMC.unit = fMC.unit '...
%			'WHERE fMC.unit = fBC.unit AND egMC.fromunit = egBC.fromunit AND fMC.modelID = 37 AND fBC.modelID = 37 AND fMC.analyses_id = fBC.analyses_id '...
%			'AND fconstBC.modelID = 31 AND fconstMC.modelID = 31 '...
%			'AND NOT EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et.`1DBCrecording` AND bci.unit = fMC.unit) '...
%			'AND fMC.unit = "' bci_unit{1} '" '...
%			'AND et.`manualrecording` = "' mcfile '"']));
%			if strcmp(GC.Data, 'No Data')
%				skipthisone = 1;
%				continue;
%			end		
%			othergranger(j,1:2) = cell2mat(GC.Data(:,5:6));
%		end
%		if skipthisone == 1
%			continue;
%		end
%	
%	%'AND egMC.fromunit = "' unit '" '...
%		diffcotuned = diff(cotunedgranger, 1,2);
%		diffother = diff(othergranger, 1,2);
%		%pause 
%	
%		dcotuned = [dcotuned(:); diffcotuned(:)];
%		dother = [dother(:); diffother(:)];
%	end
%end
%
%save('./scripts/fig3b+c_part2.mat')
load('./scripts/fig3b+c_part2.mat')

dcotunedMCDC = dcotuned;
dotherMCDC = dother; 

[hMCBC, pMCBC] = ttest2(abs(dcotunedMCBC), abs(dotherMCBC))
[hMCDC, pMCDC] = ttest2(abs(dcotunedMCDC), abs(dotherMCDC))

[hMCBCsgn, pMCBCsgn] = ttest2((dcotunedMCBC), (dotherMCBC))
[hMCDCsgn, pMCDCsgn] = ttest2((dcotunedMCDC), (dotherMCDC))

[pMCBCrs, hMCBCrs] = ranksum((dcotunedMCBC), (dotherMCBC))
[pMCDCrs, hMCDCrs] = ranksum((dcotunedMCDC), (dotherMCDC))

[pMCBCsrC, hMCBCsrC] = signrank(dcotunedMCBC)
[pMCBCsrO, hMCBCsrO] = signrank(dotherMCBC)
[pMCDCsrC, hMCDCsrC] = signrank(dcotunedMCDC)
[pMCDCsrO, hMCDCsrO] = signrank(dotherMCDC)

[pBCDCrsC, hBCDCrsC] = ranksum((dcotunedMCBC), (dcotunedMCDC))
[pBCDCrsO, hBCDCrsO] = ranksum((dotherMCBC), (dotherMCDC))

figure
bar([mean((dcotunedMCBC)), mean((dotherMCBC)), mean((dcotunedMCDC)), mean((dotherMCDC))]);
hold on 
title(['MCBC: (dc)vs (do) p-value: ' num2str(pMCBCrs) ', MCDC: (dc)vs (do) p-value: ' num2str(pMCDCrs)])

errorbar([mean(dcotunedMCBC), mean(dotherMCBC), mean(dcotunedMCDC), mean(dotherMCDC)],...
	[std(dcotunedMCBC)/sqrt(length(dcotunedMCBC)), std(dotherMCBC)/sqrt(length(dotherMCBC)), std(dcotunedMCDC)/sqrt(length(dcotunedMCDC)), std(dotherMCDC)/sqrt(length(dotherMCDC))]);

saveplot(gcf, './figures/TE-decoupling-bargraph_bootstrap_rotated_signed_sem_noncontrol.eps')

save('bcitTEchanges_bootstrap_rotated_noncontrol.mat')