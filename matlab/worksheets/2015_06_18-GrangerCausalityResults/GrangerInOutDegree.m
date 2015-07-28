%In degree
%Fetch paramcode to load
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/Spanky');
data = exec(conn, ['SELECT `f`.`nev file` AS `nev file`,`gcbc`.`unit` AS `unit`,avg(`ge`.`score`) AS mcgc, `gcbc`.`AVG(ge.score)` AS bcgc 
	FROM
	(((`Average in Granger score for BCI units in brain recordings` `gcbc`
	INNER JOIN `AnalysisLinear` `al` 
	ON ((`al`.`1DBCrecording` = `gcbc`.`nev file`))) 
	INNER JOIN `Fits` `f`
	ON(((`al`.`manualrecording` = `f`.`nev file`) AND (`f`.`unit` = `gcbc`.`unit`)))) 
	INNER JOIN `GrangerEstimates` `ge` 
	ON ((`ge`.`fitID` = `f`.`fitID`))) 
	GROUP BY `gcbc`.`unit`,`f`.`nev file`']);
data = fetch(data);
data = data.Data;

hold on 
x = [];
y = [];
%c = [];
for idx = 1:(size(data,1))
	x = [x, data{idx,3}];
	y = [y, data{idx,4}];
%		c = [c, (data{idx,2}+1)-(data{idx+1,2}+1)];
end
clf
scatter(x,y,'filled');
hold on
plot([0 30], [0 30], 'r')
xlabel('manual control')
ylabel('brain control')
xlim([0 30])
ylim([0 30])
title('Average in-degree GC for BCI units')
saveplot(gcf, './worksheets/2015_06_18/mcgcVsbcgc_indegree.eps')

%Out degree
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/Spanky');
data = exec(conn, ['select f.`nev file`, ge.fromunit, AVG(ge.score) mcgc, bcgc.`avg(ge.score)` bcgcs '...
'FROM '...
'`Average GC from BCI units in brain control recording` bcgc '...
'INNER JOIN `AnalysisLinear` al ON bcgc.`nev file` = al.`1DBCrecording` '...
'INNER JOIN Fits f ON f.`nev file` = al.manualrecording '...
'INNER JOIN GrangerEstimates ge ON ge.fromunit = bcgc.`fromunit` AND f.fitID = ge.fitID '...
'WHERE f.modelID = 2 '...
'GROUP BY f.`nev file`, ge.fromunit']);
data = fetch(data);
data = data.Data;

hold on 
x = [];
y = [];
%c = [];
for idx = 1:(size(data,1))
	x = [x, data{idx,3}];
	y = [y, data{idx,4}];
%		c = [c, (data{idx,2}+1)-(data{idx+1,2}+1)];
end
clf
scatter(x,y,'filled');
hold on
plot([0 30], [0 30], 'r')
xlabel('manual control')
ylabel('brain control')
xlim([0 30])
ylim([0 30])
title('Average out-degree GC for BCI units')
saveplot(gcf, './worksheets/2015_06_18/mcgcVsbcgc_outdegree.eps')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BC vs non-BC units in a BC recording%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%f.`nev file` AS `nev file`, f.`unit` AS `unit`, ge.`score`
bciin = exec(conn, ['SELECT f.`nev file`, f.unit, ge.score, ge.fromunit '...
'FROM '...
'GrangerEstimates ge '...
'INNER JOIN Fits f '...
'ON  f.id = ge.id '...
'INNER JOIN BCIUnits bci '...
'ON bci.ID = f.`nev file` AND bci.`unit` = f.`unit` '...
'AND f.modelID = 2']);
bciin = fetch(bciin);
bciin = bciin.Data;

nonbciin = exec(conn, ['SELECT f.`nev file`, f.unit, ge.score, ge.fromunit '...
'FROM '...
'GrangerEstimates ge '...
'INNER JOIN Fits f '...
'ON  f.id = ge.id '...
'WHERE f.unit NOT IN (SELECT bc.unit FROM BCIUnits bc WHERE bc.ID = f.`nev file`) '...
'AND f.`nev file` IN (SELECT ID FROM BCIUnits) '...
'AND f.`modelID` = 2'])
nonbciin = fetch(nonbciin);
nonbciin = nonbciin.Data;

%Group into the nev files, normalize by std
nonbcibynev = containers.Map;
nU = size(nonbciin,1);
for idx = 1:nU
	nevfile = nonbciin{idx,1};
	score = nonbciin{idx,3};
	dictstr = nevfile;
	if isKey(nonbcibynev, dictstr)
		vals = nonbcibynev(dictstr);
		nonbcibynev(dictstr) = [vals, score];
	else
		nonbcibynev(dictstr) = [];
	end
end

bcibynev = containers.Map;
nU = size(bciin,1);
for idx = 1:nU
	nevfile = bciin{idx,1};
	score = bciin{idx,3};
	dictstr = nevfile;
	if isKey(bcibynev, dictstr)
		vals = bcibynev(dictstr);
		bcibynev(dictstr) = [vals, score];
	else
		bcibynev(dictstr) = [];
	end
end

%Normalize by std
zscores = containers.Map;
keys = bcibynev.keys();

for idx = 1:length(bcibynev)
	nevfile = keys{idx};
	nonbcidata = nonbcibynev(nevfile);
	bcidata = bcibynev(nevfile);
	
	zscores(nevfile) = (bcidata-mu)./sigmaa;

end

%Sort by quantile...
quantiles = containers.Map;
for idx = 1:length(bcibynev)
	nevfile = keys{idx};
	nonbcidata = nonbcibynev(nevfile);
	bcidata = bcibynev(nevfile);
	tosort = [nonbcidata', zeros(size(nonbcidata')); bcidata', ones(size(bcidata'))];
	tosort = sortrows(tosort);
	quantiles(nevfile) = mean(find(tosort(:,2)))/size(tosort,1);
end