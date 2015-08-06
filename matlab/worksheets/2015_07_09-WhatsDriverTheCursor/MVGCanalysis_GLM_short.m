
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/Spanky');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Proportion of GC score versus GLM tuning%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = exec(conn, ['SELECT ge.fromunit, ge.score, rec.successrate, f.`nev file`,'...
	' f1.`mse out`, f1.`dev`, f2.`mse out`, f2.`dev`, f3.`mse out`, f3.`dev`'...
	'FROM GrangerEstimates ge '...
	'INNER JOIN Fits f ON f.id = ge.id '...
	'INNER JOIN BCIUnits bci ON bci.ID = f.`nev file` AND bci.unit = ge.fromunit '...
	'INNER JOIN Recordings rec ON rec.`nev file` = f.`nev file` '...
	'INNER JOIN AnalysisLinear al ON al.`1DBCrecording` = f.`nev file` '...
	'INNER JOIN Fits f1 ON f1.`nev file` = al.`manualrecording` AND f1.`modelID` = 14 and f1.unit = ge.fromunit '...
	'INNER JOIN Fits f2 ON f2.`nev file` = al.`manualrecording` AND f2.`modelID` = 15 and f2.unit = ge.fromunit '...
	'INNER JOIN Fits f3 ON f3.`nev file` = al.`manualrecording` AND f3.`modelID` = 1 and f3.unit = ge.fromunit '...
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
		performance(idx) = data{i, 3};
	end
end

figure
x = torqmseout1; y = torqmseout2; c = gcnorm; a = 100*performance+20;
scatter(x,y,[],c, 'filled');
hold on
plot([0 .4], [0 .4], 'k')
ylim([0 1])
xlim([0 1])
xlabel('MSE. Unit 1')
ylabel('MSE. Unit 2')
title('\Delta Granger (Brain control)')
caxis([-.6 .6])
colorbar
saveplot(gcf, './worksheets/2015_07_09/GLM_torque_mseoutVsGranger.eps')

figure
plot([linmseout1, linmseout2], [torqmseout1, torqmseout2], '.')
xlabel('Linear MSE out')
ylabel('GLM(torque) MSE out')
saveplot(gcf, './worksheets/2015_07_09/GLM-LinMSEout.eps')

figure
x = deldev1; y = deldev2; c = gcnorm; a = 100*performance+20;
scatter(x,y,[],c, 'filled');
hold on
plot([0 .4], [0 .4], 'k')
%ylim([0 0.05])
%xlim([0 0.05])
xlabel('\Delta Dev Unit 1')
ylabel('\Delta Dev Unit 2')
title('\Delta Granger (Brain control)')
caxis([-.6 .6])
colorbar
saveplot(gcf, './worksheets/2015_07_09/GLM_torque-sp_devVsGranger.eps')

figure
x = torqmseout1; y = torqmseout2; c = performance; a = min(gcnorm)+10+10*gcnorm;
scatter(x,y,[],c, 'filled');
hold on
%plot([0 4], [0 4])
xlabel('MSE. Unit 1')
ylabel('MSE. Unit 2')
title('Performance (successes/second)')
colorbar
saveplot(gcf, './worksheets/2015_07_09/GLM_torque_mseoutvsPerformance.eps')

figure
x = mseout1; y = mseout2; c = performance;
plot([x y], [performance performance], '.');
%plot([0 4], [0 4])
xlabel('mse')
ylabel('performance')
saveplot(gcf, './worksheets/2015_07_09/mseoutvsPerformanceA.eps')

%Plot performance versus balance-ness:
figure
plot(abs(gcnorm), performance, '.');
xlabel('GC balance')
ylabel('Performance (successes/min)')
saveplot(gcf, './worksheets/2015_07_09/GCvsperformance.eps')