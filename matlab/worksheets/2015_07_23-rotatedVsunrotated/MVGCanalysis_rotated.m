conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

data2013 = exec(conn, ['SELECT ge.fromunit, ge.score, rec.successrate, f.`nev file`,'...
	' fl.dir, bci.angle, MOD(fl.dir-bci.angle, 2*PI()), fl.size, u.firingrate, f1.`mse out`, f1.`dev`, f2.`mse out`, f2.`dev`, '...
	'ai.`im tuning type` '...
	'FROM estimates_granger ge '...
	'INNER JOIN fits f ON f.id = ge.id '...
	'INNER JOIN recordings rec ON rec.`nev file` = f.`nev file` '...
	'INNER JOIN bci_units bci ON bci.ID = f.`nev file` AND bci.unit = ge.fromunit '...
	'INNER JOIN experiment_tuning al ON al.`1DBCrecording` = f.`nev file` '...
	'INNER JOIN Analysis1DBCIvana ai ON ai.`nev file` = f.`nev file` '...
	'INNER JOIN fits f1 ON f1.`nev file` = al.`manualrecording` AND f1.`modelID` = 1 and f1.unit = ge.fromunit '...
	'INNER JOIN fits_linear fl ON f1.`id` = fl.id '...
	'INNER JOIN units u ON u.`nev file` = f1.`nev file` AND u.`unit` = f1.unit '...
	'INNER JOIN fits f2 ON f2.`nev file` = al.`1DBCrecording` AND f2.`modelID` = 1 AND f2.unit = ge.fromunit WHERE f.modelID = 3 ']);

data2013 = fetch(data2013);
data2013 = data2013.Data;
for idx = 1:size(data2013,1)
	data2013{idx,14} = data2013{idx,14}~=5;
end

data2014 = exec(conn, ['SELECT ge.fromunit, ge.score, rec.successrate, f.`nev file`,'...
	' fl.dir, bci.angle, MOD(fl.dir-bci.angle, 2*PI()), fl.size, u.firingrate, f1.`mse out`, f1.`dev`, f2.`mse out`, f2.`dev`, '...
	'a2.`Type of trial`'...
	'FROM estimates_granger ge '...
	'INNER JOIN fits f ON f.id = ge.id '...
	'INNER JOIN recordings rec ON rec.`nev file` = f.`nev file` '...
	'INNER JOIN bci_units bci ON bci.ID = f.`nev file` AND bci.unit = ge.fromunit '...
	'INNER JOIN experiment_tuning al ON al.`1DBCrecording` = f.`nev file` '...
	'INNER JOIN `rotated_linear_analysis_2014` a2 ON a2.`File name BC` = f.`nev file` '...
	'INNER JOIN fits f1 ON f1.`nev file` = al.`manualrecording` AND f1.`modelID` = 1 and f1.unit = ge.fromunit '...
	'INNER JOIN fits_linear fl ON f1.`id` = fl.id '...
	'INNER JOIN units u ON u.`nev file` = f1.`nev file` AND u.`unit` = f1.unit '...
	'INNER JOIN fits f2 ON f2.`nev file` = al.`1DBCrecording` AND f2.`modelID` = 1 AND f2.unit = ge.fromunit WHERE f.modelID = 3 ']);

data2014 = fetch(data2014);
data2014 = data2014.Data;
for idx = 1:size(data2014,1)
	data2014{idx,14} = strcmp(data2014{idx,14}, 'rotated');
end

data = [data2013; data2014];

gc = [];
gc1 = [];
gc2 = [];
gcnorm = [];
delangle1 = [];
tuningsize1 = [];
firing1 = [];
delangle2 = [];
tuningsize2 = [];
firing2 = [];
performance = [];
mseout1 = [];
mseout2 = [];
dev1 = [];
dev2 = [];
dev1BC = [];
dev2BC = [];
mseout1BC = [];
mseout2BC = [];
rotated = [];

for i = 1:(size(data,1)-1)
	if strcmp(data{i, 4}, data{i+1, 4})
		idx = length(gc)+1;
		%delangle(idx+1) = data{i+1,7};
		gc1(idx) = data{i,2};
		gc2(idx) = data{i+1,2};
		total = data{i,2}+data{i+1,2};
		gc(idx) = data{i,2}-data{i+1,2};
		gcnorm(idx) = gc(idx)/total;
		delangle1(idx) = data{i,7};
		delangle2(idx) = data{i+1,7};
		tuningsize1(idx) = data{i,8};
		tuningsize2(idx) = data{i+1,8};
		firing1(idx) = data{i,9};
		firing2(idx) = data{i+1,9};
		mseout1(idx) = data{i,10};
		mseout2(idx) = data{i+1,10};
		dev1(idx) = data{i,11};
		dev2(idx) = data{i+1,11};
		mseout1BC(idx) = data{i,12};
		mseout2BC(idx) = data{i+1,12};
		dev1BC(idx) = data{i,13};
		dev2BC(idx) = data{i+1,13};
		performance(idx) = data{i, 3};
		rotated(idx) = data{i,14};
	end
end

delangle1 = mod(delangle1, 2*pi);
delangle1(delangle1>pi) = delangle1(delangle1>pi)-2*pi;
%Make between pi and -pi
delangle2 = mod(delangle2, 2*pi);
delangle2(delangle2>pi) = delangle2(delangle2>pi)-2*pi;

figure
subplot(1,2,1)
x = mseout1(rotated~=1); y = mseout2(rotated~=1); c = gcnorm(rotated~=1); a = 20*rotated+20;
scatter(x,y,[],c, 'filled');
hold on
plot([0 .4], [0 .4], 'k')
xlabel('MSE. Unit 1')
ylabel('MSE. Unit 2')
title('\Delta Granger (Brain control) (Unrotated)')
caxis([-.6 .6])
colorbar

subplot(1,2,2)
x = mseout1(rotated==1); y = mseout2(rotated==1); c = gcnorm(rotated==1); a = 20*rotated+20;
scatter(x,y,[],c, 'filled');
hold on
plot([0 .4], [0 .4], 'k')
xlabel('MSE. Unit 1')
ylabel('MSE. Unit 2')
title('\Delta Granger (Brain control) (Rotated)')
caxis([-.6 .6])
colorbar
saveplot(gcf, './worksheets/2015_07_23-rotatedVsunrotated/mseoutVsGranger_rotated.eps')

figure
subplot(121)
x = gc1(rotated~=1); y = gc2(rotated~=1); c = performance(rotated~=1);
scatter(x,y,[],c, 'filled');
hold on
%plot([0 4], [0 4])
xlabel('Granger score. Unit 1')
ylabel('Granger score. Unit 2')
title('Performance (successes/second) (unrotated)')
colorbar
subplot(122)
x = gc1(rotated==1); y = gc2(rotated==1); c = performance(rotated==1);
scatter(x,y,[],c, 'filled');
hold on
%plot([0 4], [0 4])
xlabel('Granger score. Unit 1')
ylabel('Granger score. Unit 2')
title('Performance (successes/second) (rotated)')
colorbar
saveplot(gcf, './worksheets/2015_07_23-rotatedVsunrotated/GrangervsPerformance_rotated.eps')

figure
subplot(121)
unrotatedGC = [gc1(rotated~=1), gc2(rotated~=1)]
rotatedGC = [gc1(rotated==1), gc2(rotated==1)]
hist(unrotatedGC);
set(get(gca,'child'),'FaceColor',[0 0.5 0.5],'EdgeColor','w');
xlabel('GC')
title('unrotated')
subplot(122)
hist(rotatedGC)
title('rotated')
saveplot(gcf, './worksheets/2015_07_23-rotatedVsunrotated/Granger_rotated.eps')

figure
unrotatedperf = [performance(rotated~=1)]
rotatedperf = [performance(rotated==1)]
subplot(121)
hist(unrotatedperf);
set(get(gca,'child'),'FaceColor',[0 0.5 0.5],'EdgeColor','w');
xlabel('Performance')
title('unrotated')
subplot(122)
hist(rotatedperf)
xlabel('Performance')
title('rotated')
saveplot(gcf, './worksheets/2015_07_23-rotatedVsunrotated/Performance_rotated.eps')

figure
subplot(121)
x = mseout1(rotated==1); y = mseout2(rotated==1); c = performance(rotated==1); a = min(gcnorm)+10+10*gcnorm;
scatter(x,y,[],c, 'filled');
hold on
%plot([0 4], [0 4])
xlabel('MSE. Unit 1')
ylabel('MSE. Unit 2')
title('Performance (successes/second) (rotated)')
colorbar
subplot(122)
x = mseout1(rotated~=1); y = mseout2(rotated~=1); c = performance(rotated~=1); a = min(gcnorm)+10+10*gcnorm;
scatter(x,y,[],c, 'filled');
hold on
%plot([0 4], [0 4])
xlabel('MSE. Unit 1')
ylabel('MSE. Unit 2')
title('Performance (successes/second) (unrotated)')
colorbar
saveplot(gcf, './worksheets/2015_07_23-rotatedVsunrotated/mseoutvsPerformance.eps')