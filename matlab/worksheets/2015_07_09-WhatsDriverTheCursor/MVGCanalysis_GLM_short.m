
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Proportion of GC score versus GLM tuning%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = exec(conn, ['SELECT ge.fromunit, ge.score, rec.successrate, f.`nev file`,'...
	' f1.`mse out`, f1.`dev`, f2.`mse out`, f2.`dev`, f3.`mse out`, f3.`dev`, f4.`mse out`, f4.`dev`'...
	'FROM estimates_granger ge '...
	'INNER JOIN fits f ON f.id = ge.id '...
	'INNER JOIN bci_units bci ON bci.ID = f.`nev file` AND bci.unit = ge.fromunit '...
	'INNER JOIN recordings rec ON rec.`nev file` = f.`nev file` '...
	'INNER JOIN experiment_tuning al ON al.`1DBCrecording` = f.`nev file` '...
	'INNER JOIN fits f1 ON f1.`nev file` = al.`manualrecording` AND f1.`modelID` = 14 and f1.unit = ge.fromunit '...
	'INNER JOIN fits f2 ON f2.`nev file` = al.`manualrecording` AND f2.`modelID` = 15 and f2.unit = ge.fromunit '...
	'INNER JOIN fits f3 ON f3.`nev file` = al.`manualrecording` AND f3.`modelID` = 1 and f3.unit = ge.fromunit '...
	'INNER JOIN fits f4 ON f4.`nev file` = al.`manualrecording` AND f4.`modelID` = 16 and f4.unit = ge.fromunit '...
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
	end
end

figure
x = torqmseout1; y = torqmseout2; c = gcnorm; a = 100*performance+20;
scatter(x,y,[],c, 'filled');
hold on
plot([0 .4], [0 .4], 'k')
ylim([0 2])
xlim([0 2])
xlabel('MSE. Unit 1')
ylabel('MSE. Unit 2')
title('\Delta Granger (Brain control)')
caxis([-.6 .6])
colorbar
saveplot(gcf, './worksheets/2015_07_09-WhatsDriverTheCursor/GLM_torque_mseoutVsGranger.eps')

figure
x = spmseout1; y = spmseout2; c = gcnorm; a = 100*performance+20;
scatter(x,y,[],c, 'filled');
hold on
plot([0 .4], [0 .4], 'k')
ylim([0 2])
xlim([0 2])
xlabel('MSE. Unit 1')
ylabel('MSE. Unit 2')
title('\Delta Granger (Brain control)')
caxis([-.6 .6])
colorbar
saveplot(gcf, './worksheets/2015_07_09-WhatsDriverTheCursor/GLM_sp_mseoutVsGranger.eps')

figure
plot([linmseout1, linmseout2], [torqmseout1, torqmseout2], '.')
xlabel('Linear MSE out')
ylabel('GLM(torque) MSE out')
saveplot(gcf, './worksheets/2015_07_09-WhatsDriverTheCursor/GLM-LinMSEout.eps')

figure
plot([spmseout1, spmseout2], [torqmseout1, torqmseout2], '.')
xlabel('GLM(sp) MSE out')
ylabel('GLM(torque) MSE out')
title('53/204 have MSE_{sp} < MSE_{torq}')
saveplot(gcf, './worksheets/2015_07_09-WhatsDriverTheCursor/GLM-SpMSEoutVsGLM-TorqMSEOUT.eps')

figure
hist([spmseout1, spmseout2]-[torqmseout1, torqmseout2])
xlabel('GLM(sp)-GLM(torq) MSE out')
title('53/204 have MSE_{sp} < MSE_{torq}')
saveplot(gcf, './worksheets/2015_07_09-WhatsDriverTheCursor/GLM-SpMSEoutVsGLM-TorqMSEOUT_hist.eps')

figure
plot([spdev1, spdev2], [torqdev1, torqdev2], '.')
xlabel('GLM(sp) dev')
ylabel('GLM(torque) dev')
title('All dev_{sp} > dev{torq} (as expected)')
saveplot(gcf, './worksheets/2015_07_09-WhatsDriverTheCursor/GLM-LinMSEout.eps')

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
saveplot(gcf, './worksheets/2015_07_09-WhatsDriverTheCursor/GLM_torque-sp_devVsGranger.eps')

figure
x = torqmseout1; y = torqmseout2; c = performance; a = min(gcnorm)+10+10*gcnorm;
scatter(x,y,[],c, 'filled');
hold on
%plot([0 4], [0 4])
xlabel('MSE. Unit 1')
ylabel('MSE. Unit 2')
title('Performance (successes/second)')
colorbar
saveplot(gcf, './worksheets/2015_07_09-WhatsDriverTheCursor/GLM_torque_mseoutvsPerformance.eps')

figure
x = torqmseout1; y = torqmseout2; c = performance;
plot([x y], [performance performance], '.');
%plot([0 4], [0 4])
xlabel('mse')
ylabel('performance')
saveplot(gcf, './worksheets/2015_07_09-WhatsDriverTheCursor/mseoutvsPerformanceA.eps')

%%%%%%%%%%%%%%%%%%%%%
%Model 14 - model 15%
%%%%%%%%%%%%%%%%%%%%%

figure
x = deldev1; y = deldev2; c = performance;
plot([x y], [performance performance], '.');
%plot([0 4], [0 4])
xlabel('change in dev between w and w.out torque')
ylabel('performance')
saveplot(gcf, './worksheets/2015_07_09-WhatsDriverTheCursor/delDevvsPerformanceA.eps')

figure
x = deldev1; y = deldev2; c = performance; a = min(gcnorm)+10+10*gcnorm;
scatter(x,y,[],c, 'filled');
hold on
%plot([0 4], [0 4])
xlabel('MSE. Unit 1')
ylabel('MSE. Unit 2')
title('Performance (successes/second)')
colorbar
saveplot(gcf, './worksheets/2015_07_09-WhatsDriverTheCursor/delDevvsPerformance.eps')


%Plot performance versus balance-ness:
figure
plot(abs(gcnorm), performance, '.');
xlabel('GC balance')
ylabel('Performance (successes/min)')
saveplot(gcf, './worksheets/2015_07_09-WhatsDriverTheCursor/GCvsperformance.eps')

%%%%%%%%%%
%Model 16%
%%%%%%%%%%

figure
x = targmseout1; y = targmseout2; c = gcnorm; a = 100*performance+20;
scatter(x,y,[],c, 'filled');
hold on
plot([0 .4], [0 .4], 'k')
ylim([0 2])
xlim([0 2])
xlabel('MSE. Unit 1')
ylabel('MSE. Unit 2')
title('\Delta Granger (Brain control)')
caxis([-.6 .6])
colorbar
saveplot(gcf, './worksheets/2015_07_09-WhatsDriverTheCursor/GLM_target_mseoutVsGranger.eps')

figure
x = spdev1-targdev1; y = spdev2-targdev2; c = gcnorm; a = 100*performance+20;
scatter(x,y,[],c, 'filled');
hold on
plot([0 .4], [0 .4], 'k')
%ylim([0 2])
%xlim([0 2])
xlabel('MSE. Unit 1')
ylabel('MSE. Unit 2')
title('\Delta Granger (Brain control)')
caxis([-.6 .6])
colorbar
saveplot(gcf, './worksheets/2015_07_09-WhatsDriverTheCursor/GLM_target_dev-GLM_sp_dev_VsGranger.eps')

%Fit a surface to the data 
figure
f = fit([x', y'], c', 'cubicinterp')
plot(f)
hold on 
scatter(x,y,[],c, 'filled');
