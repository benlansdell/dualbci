conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%List of files
files = fetch(exec(conn, ['SELECT et.`manualrecording` FROM experiment_tuning et']));
files = files.Data;

%Torque tuning angle (velocity)
deltaBCI = [];
deltacotuned = [];
deltaother = [];

dothercotuned = [];
dcotuned = [];

nR = 50;

for rep = 1:nR
	display(['rep = ' num2str(rep)])
	for idx = 1:length(files)
		cotunedgranger = [];
		othergranger = [];
		mcfile = files{idx};
		all_data = fetch(exec(conn, ['SELECT fl1.dir, fl2.dir, fl5.dir, et1.`tuning_type`, flin1.unit FROM '...
		'`experiment_tuning` et1 '...
		'INNER JOIN `fits` flin1 '...
		'ON flin1.`nev file` = et1.`manualrecording`'...
		'INNER JOIN `fits_linear` fl1 '...
		'ON flin1.id = fl1.id '...
		'INNER JOIN `fits` flin2 '...
		'ON flin2.`nev file` = et1.`1DBCrecording`'...
		'INNER JOIN `fits_linear` fl2 '...
		'ON flin2.id = fl2.id '...
		'INNER JOIN `fits` flin5 '...
		'ON flin5.`nev file` = et1.`dualrecording` '...
		'INNER JOIN `fits_linear` fl5 '...
		'ON flin5.id = fl5.id '...
		'WHERE flin1.modelID = 30 AND flin2.modelID = 30 AND flin5.modelID = 30 ' ...
		'AND flin1.unit = flin2.unit AND flin2.unit = flin5.unit ' ...
		'AND fl1.r2 > 0.0 '...
		'AND NOT EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et1.`1DBCrecording` AND bci.unit = flin1.unit) '...
		'AND et1.`manualrecording` = "' mcfile '" '...
		'AND (et1.`tuning_type` = 1 OR et1.`tuning_type` = 3 OR et1.`tuning_type` = 4)']));
		if strcmp(all_data.Data, 'No Data')
			continue 
		end
		if size(all_data.Data, 1) < 4 
			continue 
		end
	
		all_theta = cell2mat(all_data.Data(:,1:3));
		tuningtype = cell2mat(all_data.Data(:,4));
		otherunits = all_data.Data(:,5);
		
		%BCI unit
		all_data = fetch(exec(conn, ['SELECT fl1.dir, fl2.dir, fl5.dir, flin1.unit FROM '...
		'`experiment_tuning` et1 '...
		'INNER JOIN `fits` flin1 '...
		'ON flin1.`nev file` = et1.`manualrecording`'...
		'INNER JOIN `fits_linear` fl1 '...
		'ON flin1.id = fl1.id '...
		'INNER JOIN `fits` flin2 '...
		'ON flin2.`nev file` = et1.`1DBCrecording`'...
		'INNER JOIN `fits_linear` fl2 '...
		'ON flin2.id = fl2.id '...
		'INNER JOIN `fits` flin5 '...
		'ON flin5.`nev file` = et1.`dualrecording`'...
		'INNER JOIN `fits_linear` fl5 '...
		'ON flin5.id = fl5.id '...
		'WHERE flin1.modelID = 30 AND flin2.modelID = 30 AND flin5.modelID = 30 ' ...
		'AND flin1.unit = flin2.unit AND flin2.unit = flin5.unit ' ...
		'AND EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et1.`1DBCrecording` AND bci.unit = flin1.unit) '...
		'AND et1.`manualrecording` = "' mcfile '" '...
		'LIMIT 1']));
		%'AND et1.`tuning_type` = 5 '...
		bci_theta = cell2mat(all_data.Data(:,1:3));
		bci_unit = all_data.Data(:,4);
		nU = size(all_theta,1);
		
		%Pick cotuned units to BCI units in MC
		diff_MC_theta = cos(bci_theta(1) - all_theta(:,1));
		
		%Pick random two
		[l, cotunedidx] = sort(diff_MC_theta, 'descend'); 
		cotunedidx = randsample(1:nU,2);

		bcicotuned = diff_MC_theta(cotunedidx);
		cotunedunits = otherunits(cotunedidx);
	
		%Pick random unit
		otheridx = randsample(1:nU, 1);
		otherunitname = otherunits(otheridx);

		%Pick another random two units
		othercotunedidx = randsample(setdiff(1:nU,otheridx),2);

		othercotunedunits = otherunits(othercotunedidx);
	
		%Granger for cotuned units
		skipthisone = 0;
		for j = 1:length(cotunedunits)
			unit = cotunedunits{j};
			GC = fetch(exec(conn, ['SELECT fBC.analyses_id, fMC.`nev file`, fMC.unit tounit, '...
			' egMC.fromunit, egMC.score egMCscore, egBC.score egBCscore, fconstMC.dev, fconstBC.dev, et.`tuning_type` FROM `experiment_tuning` et '...
			'INNER JOIN `fits` fMC '...
			'ON et.`manualrecording` = fMC.`nev file` '...
			'INNER JOIN `estimates_granger` egMC '...
			'ON egMC.`id` = fMC.`id` '...
			'INNER JOIN `fits` fBC '...
			'ON et.`1DBCrecording` = fBC.`nev file` '...
			'INNER JOIN `estimates_granger` egBC '...
			'ON egBC.`id` = fBC.`id` '...
			'INNER JOIN `fits` fconstBC '...
			'ON fconstBC.`nev file` = et.`1DBCrecording` AND fconstBC.unit = fMC.unit '...
			'INNER JOIN `fits` fconstMC '...
			'ON fconstMC.`nev file` = et.`manualrecording` AND fconstMC.unit = fMC.unit '...
			'WHERE fMC.unit = fBC.unit AND egMC.fromunit = egBC.fromunit AND fMC.modelID = 29 AND fBC.modelID = 29 AND fMC.analyses_id = fBC.analyses_id '...
			'AND fconstBC.modelID = 31 AND fconstMC.modelID = 31 '...
			'AND EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et.`1DBCrecording` AND bci.unit = fMC.unit) '...
			'AND egMC.fromunit = "' unit '" '...
			'AND fMC.unit = "' bci_unit{1} '" '...
			'AND et.`manualrecording` = "' mcfile '"']));
			if strcmp(GC.Data, 'No Data')
				skipthisone = 1;
				continue;
			end		
			cotunedgranger(j,1:2) = cell2mat(GC.Data(:,5:6));
		end
		if skipthisone == 1
			continue;
		end
	
		%Granger for other units
		skipthisone = 0;
		for j = 1:length(othercotunedunits)
			unit = othercotunedunits{j};
			GC = fetch(exec(conn, ['SELECT fBC.analyses_id, fMC.`nev file`, fMC.unit tounit, '...
			' egMC.fromunit, egMC.score egMCscore, egBC.score egBCscore, fconstMC.dev, fconstBC.dev, et.`tuning_type` FROM `experiment_tuning` et '...
			'INNER JOIN `fits` fMC '...
			'ON et.`manualrecording` = fMC.`nev file` '...
			'INNER JOIN `estimates_granger` egMC '...
			'ON egMC.`id` = fMC.`id` AND egMC.`fromunit` = "' unit '" ' ...
			'INNER JOIN `fits` fBC '...
			'ON et.`1DBCrecording` = fBC.`nev file` '...
			'INNER JOIN `estimates_granger` egBC '...
			'ON egBC.`id` = fBC.`id` '...
			'INNER JOIN `fits` fconstBC '...
			'ON fconstBC.`nev file` = et.`1DBCrecording` AND fconstBC.unit = fMC.unit '...
			'INNER JOIN `fits` fconstMC '...
			'ON fconstMC.`nev file` = et.`manualrecording` AND fconstMC.unit = fMC.unit '...
			'WHERE fMC.unit = fBC.unit AND egMC.fromunit = egBC.fromunit AND fMC.modelID = 29 AND fBC.modelID = 29 AND fMC.analyses_id = fBC.analyses_id '...
			'AND fconstBC.modelID = 31 AND fconstMC.modelID = 31 '...
			'AND fMC.unit = "' otherunitname{1} '" '...
			'AND et.`manualrecording` = "' mcfile '"']));
			if strcmp(GC.Data, 'No Data')
				skipthisone = 1;
				continue;
			end		
			othergranger(j,1:2) = cell2mat(GC.Data(:,5:6));
		end
		if skipthisone == 1
			continue;
		end
	
		diffcotuned = diff(cotunedgranger, 1,2);
		diffother = diff(othergranger, 1,2);
		dcotuned = [dcotuned(:); diffcotuned(:)];
		dothercotuned = [dothercotuned(:); diffother(:)];
	end
end

figure
rng = linspace(-100, 100, 100);
histogram(dcotuned, rng, 'Normalization', 'probability')
hold on 
histogram(dothercotuned, rng, 'Normalization', 'probability')
%From teMCBCstats.m
MCBCpctile = 0.0011;
plot([MCBCpctile, MCBCpctile], [0 1], 'k--')
plot([-MCBCpctile, -MCBCpctile], [0 1], 'k--')
xlabel('MCGC - BCGC')
ylabel('Count')
%xlim([min(rng) max(rng)])
ylim([0, 0.3])
legend('With BC-unit', 'With randomly selected unit')
[h, p] = ttest2(abs(dcotuned), abs(dothercotuned))
100*(sum(dcotuned < -MCBCpctile) + sum(dcotuned > MCBCpctile))/length(dcotuned)
100*(sum(dothercotuned < -MCBCpctile) + sum(dothercotuned > MCBCpctile))/length(dothercotuned)

title(['abs(dbci)vs abs(dother) p-value: ' num2str(p)])
saveplot(gcf, './worksheets/2016_06_10-resultsforpaper/Granger-decoupling-MC-BC_histogram_bootstrap_othercontrol_rotated.eps')
%saveplot(gcf, './worksheets/2016_06_10-resultsforpaper/Granger-decoupling-MC-BC_histogram.png', 'png', [4 3])

figure 
qqplot(dcotuned, dothercotuned)
[h, p] = kstest2(dcotuned, dothercotuned)
title(['kstest2: p = ' num2str(p)])
saveplot(gcf, './worksheets/2016_06_10-resultsforpaper/Granger-decoupling-MC-BC_qqplot_bootstrap_othercontrol_rotated.eps')

%%%%%%%%%%%%%%%%%%%%
%Do the same for DC%
%%%%%%%%%%%%%%%%%%%%

dothercotuned = [];
dcotuned = [];

for rep = 1:nR
	display(['rep = ' num2str(rep)])
	for idx = 1:length(files)
		cotunedgranger = [];
		othergranger = [];
		mcfile = files{idx};
		all_data = fetch(exec(conn, ['SELECT fl1.dir, fl2.dir, fl5.dir, et1.`tuning_type`, flin1.unit FROM '...
		'`experiment_tuning` et1 '...
		'INNER JOIN `fits` flin1 '...
		'ON flin1.`nev file` = et1.`manualrecording`'...
		'INNER JOIN `fits_linear` fl1 '...
		'ON flin1.id = fl1.id '...
		'INNER JOIN `fits` flin2 '...
		'ON flin2.`nev file` = et1.`1DBCrecording`'...
		'INNER JOIN `fits_linear` fl2 '...
		'ON flin2.id = fl2.id '...
		'INNER JOIN `fits` flin5 '...
		'ON flin5.`nev file` = et1.`dualrecording` '...
		'INNER JOIN `fits_linear` fl5 '...
		'ON flin5.id = fl5.id '...
		'WHERE flin1.modelID = 30 AND flin2.modelID = 30 AND flin5.modelID = 30 ' ...
		'AND flin1.unit = flin2.unit AND flin2.unit = flin5.unit ' ...
		'AND fl1.r2 > 0.0 '...
		'AND NOT EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et1.`1DBCrecording` AND bci.unit = flin1.unit) '...
		'AND et1.`manualrecording` = "' mcfile '" '...
		'AND (et1.`tuning_type` = 1 OR et1.`tuning_type` = 3 OR et1.`tuning_type` = 4)']));
	%	'AND et1.`tuning_type` = 5']));
		if strcmp(all_data.Data, 'No Data')
			continue 
		end
		if size(all_data.Data, 1) < 4 
			continue 
		end
	
		all_theta = cell2mat(all_data.Data(:,1:3));
		tuningtype = cell2mat(all_data.Data(:,4));
		otherunits = all_data.Data(:,5);
		
		%BCI unit
		all_data = fetch(exec(conn, ['SELECT fl1.dir, fl2.dir, fl5.dir, flin1.unit FROM '...
		'`experiment_tuning` et1 '...
		'INNER JOIN `fits` flin1 '...
		'ON flin1.`nev file` = et1.`manualrecording`'...
		'INNER JOIN `fits_linear` fl1 '...
		'ON flin1.id = fl1.id '...
		'INNER JOIN `fits` flin2 '...
		'ON flin2.`nev file` = et1.`1DBCrecording`'...
		'INNER JOIN `fits_linear` fl2 '...
		'ON flin2.id = fl2.id '...
		'INNER JOIN `fits` flin5 '...
		'ON flin5.`nev file` = et1.`dualrecording`'...
		'INNER JOIN `fits_linear` fl5 '...
		'ON flin5.id = fl5.id '...
		'WHERE flin1.modelID = 30 AND flin2.modelID = 30 AND flin5.modelID = 30 ' ...
		'AND flin1.unit = flin2.unit AND flin2.unit = flin5.unit ' ...
		'AND EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et1.`1DBCrecording` AND bci.unit = flin1.unit) '...
		'AND et1.`manualrecording` = "' mcfile '" '...
		'LIMIT 1']));
		%'AND et1.`tuning_type` = 5 '...
		bci_theta = cell2mat(all_data.Data(:,1:3));
		bci_unit = all_data.Data(:,4);
		nU = size(all_theta,1);
		
		%Pick cotuned units to BCI units in MC
		diff_MC_theta = cos(bci_theta(1) - all_theta(:,1));
		
		%Pick random two
		[l, cotunedidx] = sort(diff_MC_theta, 'descend'); 
		cotunedidx = randsample(1:nU,2);

		bcicotuned = diff_MC_theta(cotunedidx);
		cotunedunits = otherunits(cotunedidx);
	
		%Pick random unit
		otheridx = randsample(1:nU, 1);
		otherunitname = otherunits(otheridx);

		%Pick another random two units
		othercotunedidx = randsample(setdiff(1:nU,otheridx),2);

		othercotunedunits = otherunits(othercotunedidx);
	
		%Granger for cotuned units
		skipthisone = 0;
		for j = 1:length(cotunedunits)
			unit = cotunedunits{j};
			GC = fetch(exec(conn, ['SELECT fBC.analyses_id, fMC.`nev file`, fMC.unit tounit, '...
			' egMC.fromunit, egMC.score egMCscore, egBC.score egBCscore, fconstMC.dev, fconstBC.dev, et.`tuning_type` FROM `experiment_tuning` et '...
			'INNER JOIN `fits` fMC '...
			'ON et.`manualrecording` = fMC.`nev file` '...
			'INNER JOIN `estimates_granger` egMC '...
			'ON egMC.`id` = fMC.`id` '...
			'INNER JOIN `fits` fBC '...
			'ON et.`dualrecording` = fBC.`nev file` '...
			'INNER JOIN `estimates_granger` egBC '...
			'ON egBC.`id` = fBC.`id` '...
			'INNER JOIN `fits` fconstBC '...
			'ON fconstBC.`nev file` = et.`dualrecording` AND fconstBC.unit = fMC.unit '...
			'INNER JOIN `fits` fconstMC '...
			'ON fconstMC.`nev file` = et.`manualrecording` AND fconstMC.unit = fMC.unit '...
			'WHERE fMC.unit = fBC.unit AND egMC.fromunit = egBC.fromunit AND fMC.modelID = 29 AND fBC.modelID = 29 AND fMC.analyses_id = fBC.analyses_id '...
			'AND fconstBC.modelID = 31 AND fconstMC.modelID = 31 '...
			'AND EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et.`dualrecording` AND bci.unit = fMC.unit) '...
			'AND egMC.fromunit = "' unit '" '...
			'AND fMC.unit = "' bci_unit{1} '" '...
			'AND et.`manualrecording` = "' mcfile '"']));
			if strcmp(GC.Data, 'No Data')
				skipthisone = 1;
				continue;
			end		
			cotunedgranger(j,1:2) = cell2mat(GC.Data(:,5:6));
		end
		if skipthisone == 1
			continue;
		end
	
		%Granger for other units
		skipthisone = 0;
		for j = 1:length(othercotunedunits)
			unit = othercotunedunits{j};
			GC = fetch(exec(conn, ['SELECT fBC.analyses_id, fMC.`nev file`, fMC.unit tounit, '...
			' egMC.fromunit, egMC.score egMCscore, egBC.score egBCscore, fconstMC.dev, fconstBC.dev, et.`tuning_type` FROM `experiment_tuning` et '...
			'INNER JOIN `fits` fMC '...
			'ON et.`manualrecording` = fMC.`nev file` '...
			'INNER JOIN `estimates_granger` egMC '...
			'ON egMC.`id` = fMC.`id` AND egMC.`fromunit` = "' unit '" ' ...
			'INNER JOIN `fits` fBC '...
			'ON et.`dualrecording` = fBC.`nev file` '...
			'INNER JOIN `estimates_granger` egBC '...
			'ON egBC.`id` = fBC.`id` '...
			'INNER JOIN `fits` fconstBC '...
			'ON fconstBC.`nev file` = et.`dualrecording` AND fconstBC.unit = fMC.unit '...
			'INNER JOIN `fits` fconstMC '...
			'ON fconstMC.`nev file` = et.`manualrecording` AND fconstMC.unit = fMC.unit '...
			'WHERE fMC.unit = fBC.unit AND egMC.fromunit = egBC.fromunit AND fMC.modelID = 29 AND fBC.modelID = 29 AND fMC.analyses_id = fBC.analyses_id '...
			'AND fconstBC.modelID = 31 AND fconstMC.modelID = 31 '...
			'AND fMC.unit = "' otherunitname{1} '" '...
			'AND et.`manualrecording` = "' mcfile '"']));
			if strcmp(GC.Data, 'No Data')
				skipthisone = 1;
				continue;
			end		
			othergranger(j,1:2) = cell2mat(GC.Data(:,5:6));
		end
		if skipthisone == 1
			continue;
		end
	
		diffcotuned = diff(cotunedgranger, 1,2);
		diffother = diff(othergranger, 1,2);
		dcotuned = [dcotuned(:); diffcotuned(:)];
		dothercotuned = [dothercotuned(:); diffother(:)];
	end
end

figure
MCDCpctile = 0.0011;
rng = linspace(-100, 100, 100);
histogram(dcotuned, rng, 'Normalization', 'probability')
hold on 
histogram(dothercotuned, rng, 'Normalization', 'probability')
plot([MCDCpctile, MCDCpctile], [0 1], 'k--')
plot([-MCDCpctile, -MCDCpctile], [0 1], 'k--')
xlabel('MCGC - DCGC')
ylabel('Count')
%xlim([min(rng) max(rng)])
ylim([0, 0.3])
legend('With BC-unit', 'With randomly selected unit')
[h, p] = ttest2(abs(dcotuned), abs(dothercotuned))
100*(sum(dcotuned < -MCDCpctile) + sum(dcotuned > MCDCpctile))/length(dcotuned)
100*(sum(dothercotuned < -MCDCpctile) + sum(dothercotuned > MCDCpctile))/length(dothercotuned)

title(['abs(dbci)vs abs(dother) p-value: ' num2str(p)])
saveplot(gcf, './worksheets/2016_06_10-resultsforpaper/Granger-decoupling-MC-DC_histogram_bootstrap_othercontrol_rotated.eps')
%saveplot(gcf, './worksheets/2016_06_10-resultsforpaper/Granger-decoupling-MC-DC_histogram.png', 'png', [4 3])

figure 
qqplot(dcotuned, dothercotuned);
[h, p] = kstest2(dcotuned, dothercotuned)
title(['kstest2: p = ' num2str(p)])
saveplot(gcf, './worksheets/2016_06_10-resultsforpaper/Granger-decoupling-MC-DC_qqplot_bootstrap_othercontrol_rotated.eps')

