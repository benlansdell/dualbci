conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Proportion of GC score versus GLM tuning%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = exec(conn, ['SELECT ge.fromunit, ge.score, rec.successrate, f.`nev file`,'...
	' f1.`mse out`, f1.`dev`, f2.`mse out`, f2.`dev`, f3.`mse out`, f3.`dev`, f4.`mse out`, f4.`dev`, '...
	' f5.`mse out`, u.firingrate '...
	'FROM estimates_granger ge '...
	'INNER JOIN fits f ON f.id = ge.id '...
	'INNER JOIN bci_units bci ON bci.ID = f.`nev file` AND bci.unit = ge.fromunit '...
	'INNER JOIN recordings rec ON rec.`nev file` = f.`nev file` '...
	'INNER JOIN units u ON u.`nev file` = f.`nev file` AND u.`unit` = bci.`unit` '...
	'INNER JOIN experiment_tuning al ON al.`1DBCrecording` = f.`nev file` '...
	'INNER JOIN fits f1 ON f1.`nev file` = al.`manualrecording` AND f1.`modelID` = 14 and f1.unit = ge.fromunit '...
	'INNER JOIN fits f2 ON f2.`nev file` = al.`manualrecording` AND f2.`modelID` = 15 and f2.unit = ge.fromunit '...
	'INNER JOIN fits f3 ON f3.`nev file` = al.`manualrecording` AND f3.`modelID` = 1 and f3.unit = ge.fromunit '...
	'INNER JOIN fits f4 ON f4.`nev file` = al.`manualrecording` AND f4.`modelID` = 16 and f4.unit = ge.fromunit '...
	'INNER JOIN fits f5 ON f5.`nev file` = al.`manualrecording` AND f5.`modelID` = 19 and f5.unit = ge.fromunit '...
	'WHERE f.modelID = 3 ']);

data = fetch(data);
data = data.Data;

gc = [];
gc1 = [];
gc2 = [];
gcnorm = [];
performance = [];
torqmseout1 = [];
torqmseout2 = [];
torqdev1 = [];
torqdev2 = [];
spmseout1 = [];
spmseout2 = [];
spdev1 = [];
spdev2 = [];
deldev1 = [];
deldev2 = [];
linmseout1 = [];
linmseout2 = [];
lindev1 = [];
lindev2 = [];
targmseout1 = [];
targmseout2 = [];
targdev1 = [];
targdev2 = [];
mseconst1 = [];
mseconst2 = [];
firingrate1 = [];
firingrate2 = [];

for i = 1:(size(data,1)-1)
	if strcmp(data{i, 4}, data{i+1, 4})
		idx = length(gc)+1;
		%delangle(idx+1) = data{i+1,7};
		gc1(idx) = data{i,2};
		gc2(idx) = data{i+1,2};
		total = data{i,2}+data{i+1,2};
		gc(idx) = data{i,2}-data{i+1,2};
		gcnorm(idx) = gc(idx)/total;
		torqmseout1(idx) = data{i,5};
		torqmseout2(idx) = data{i+1,5};
		torqdev1(idx) = data{i,6};
		torqdev2(idx) = data{i+1,6};
		spmseout1(idx) = data{i,7};
		spmseout2(idx) = data{i+1,7};
		spdev1(idx) = data{i,8};
		spdev2(idx) = data{i+1,8};
		deldev1(idx) = spdev1(idx)-torqdev1(idx);
		deldev2(idx) = spdev2(idx)-torqdev2(idx);
		linmseout1(idx) = data{i,9};
		linmseout2(idx) = data{i+1,9};
		lindev1(idx) = data{i,10};
		lindev2(idx) = data{i+1,10};
		targmseout1(idx) = data{i, 11};
		targmseout2(idx) = data{i+1, 11};
		targdev1(idx) = data{i, 12};
		targdev2(idx) = data{i+1, 12};
		performance(idx) = data{i, 3};
		mseconst1(idx) = data{i, 13};
		mseconst2(idx) = data{i+1, 13};
		firingrate1(idx) = data{i, 14};
		firingrate2(idx) = data{i+1, 14};
	end
end

%%%%%%%%%%%
%MSE const%
%%%%%%%%%%%



figure
x = mseconst1; y = mseconst2; c = gcnorm; a = 100*performance+20;
scatter(x,y,[],c, 'filled');
hold on
%plot([0 .4], [0 .4], 'k')
ylim([0 2])
xlim([0 2])
xlabel('Variance Unit 1')
ylabel('Variance Unit 2')
title('\Delta Granger (Brain control)')
caxis([-.6 .6])
colorbar
saveplot(gcf, './worksheets/2015_07_09-WhatsDriverTheCursor/Const_mseoutVsGranger.eps')

%Ok then, so what is the relation between the unit's variance and other things????
%There's a vague relationship between firing and variance, as one would expect, of course...
plot(firingrate1/25, mseconst1, '.')

%But plotting firing vs Granger scores shows less of a relationship. Note that it is still there...
figure
x = firingrate1; y = firingrate2; c = gcnorm; a = 100*performance+20;
scatter(x,y,[],c, 'filled');
hold on
plot([0 .4], [0 .4], 'k')
%ylim([0 2])
%xlim([0 2])
xlabel('Firing rate Unit 1')
ylabel('Firing rate Unit 2')
title('\Delta Granger (Brain control)')
caxis([-.6 .6])
colorbar

%How stable is firing rate, and variance between MC and BC?
data = exec(conn, ['SELECT ge.fromunit, ge.score, rec.successrate, f.`nev file`,'...
	' f5.`mse out`, f6.`mse out`, u.firingrate, u2.firingrate '...
	'FROM estimates_granger ge '...
	'INNER JOIN fits f ON f.id = ge.id '...
	'INNER JOIN bci_units bci ON bci.ID = f.`nev file` AND bci.unit = ge.fromunit '...
	'INNER JOIN recordings rec ON rec.`nev file` = f.`nev file` '...
	'INNER JOIN units u ON u.`nev file` = f.`nev file` AND u.`unit` = bci.`unit` '...
	'INNER JOIN experiment_tuning al ON al.`1DBCrecording` = f.`nev file` '...
	'INNER JOIN fits f5 ON f5.`nev file` = al.`manualrecording` AND f5.`modelID` = 19 AND f5.unit = ge.fromunit '...
	'INNER JOIN fits f6 ON f6.`nev file` = al.`1DBCrecording` AND f6.`modelID` = 19 AND f6.unit = ge.fromunit '...
	'INNER JOIN units u2 ON u2.`nev file` = al.`manualrecording` AND u2.unit = f6.unit '...
	'WHERE f.modelID = 3 AND f5.modelID = 19 AND f6.modelID = 19 '...
	'AND u2.`nev file` IN (SELECT `manualrecording` FROM experiment_tuning)']);

data = fetch(data);
data = data.Data;

gc = [];
gc1 = [];
gc2 = [];
gcnorm = [];
performance = [];
mseconst1MC = [];
mseconst2MC = [];
firingrateBC1 = [];
firingrateBC2 = [];
mseconst1BC = [];
mseconst2BC = [];
firingrateMC1 = [];
firingrateMC2 = [];
msenormBC = [];

for i = 1:(size(data,1)-1)
	if strcmp(data{i, 4}, data{i+1, 4})
		idx = length(gc)+1;
		%delangle(idx+1) = data{i+1,7};
		gc1(idx) = data{i,2};
		gc2(idx) = data{i+1,2};
		total = data{i,2}+data{i+1,2};
		gc(idx) = data{i,2}-data{i+1,2};
		gcnorm(idx) = gc(idx)/total;
		performance(idx) = data{i, 3};
		mseconstMC1(idx) = data{i, 5};
		mseconstMC2(idx) = data{i+1, 5};
		mseconstBC1(idx) = data{i, 6};
		mseconstBC2(idx) = data{i+1, 6};
		total = data{i, 6}+data{i+1, 6};
		msen = data{i, 6}-data{i+1, 6};
		msenormBC(idx) = msen/total;
		firingrateBC1(idx) = data{i, 7};
		firingrateBC2(idx) = data{i+1, 7};
		firingrateMC1(idx) = data{i, 8};
		firingrateMC2(idx) = data{i+1, 8};
	end
end

clf
plot(msenormBC, gcnorm, '.')

clf
plot([mseconstMC1, mseconstMC2], [mseconstBC1, mseconstBC2], '.')

clf
plot([firingrateMC1, firingrateMC2], [firingrateBC1, firingrateBC2], '.')

figure
x = mseconstBC1; y = mseconstBC2; c = gcnorm; a = 100*performance+20;
scatter(x,y,[],c, 'filled');
hold on
plot([0 .4], [0 .4], 'k')
%ylim([0 2])
%xlim([0 2])
xlabel('Firing rate Unit 1')
ylabel('Firing rate Unit 2')
title('\Delta Granger (Brain control)')
caxis([-.6 .6])
colorbar

clf
plot([mseconstBC1, mseconstBC2], [gc1, gc2], '.')
