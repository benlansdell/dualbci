conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Proportion of GC score versus GLM tuning%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = exec(conn, ['SELECT ge.fromunit, ge.score, rec.successrate, f.`nev file`,'...
	' f1.`mse out`, f1.`dev`, f2.`mse out`, f2.`dev`, f3.`mse out`, f3.`dev`, f4.`mse out`, f4.`dev`, '...
	' f5.`mse out`, u.firingrate, fl1.size '...
	'FROM estimates_granger ge '...
	'INNER JOIN fits f ON f.id = ge.id '...
	'INNER JOIN bci_units bci ON bci.ID = f.`nev file` AND bci.unit = ge.fromunit '...
	'INNER JOIN recordings rec ON rec.`nev file` = f.`nev file` '...
	'INNER JOIN units u ON u.`nev file` = f.`nev file` AND u.`unit` = bci.`unit` '...
	'INNER JOIN experiment_tuning al ON al.`dualrecording` = f.`nev file` '...
	'INNER JOIN fits f1 ON f1.`nev file` = al.`manualrecording` AND f1.`modelID` = 14 and f1.unit = ge.fromunit '...
	'INNER JOIN fits f2 ON f2.`nev file` = al.`manualrecording` AND f2.`modelID` = 15 and f2.unit = ge.fromunit '...
	'INNER JOIN fits f3 ON f3.`nev file` = al.`manualrecording` AND f3.`modelID` = 1 and f3.unit = ge.fromunit '...
	'INNER JOIN fits f4 ON f4.`nev file` = al.`manualrecording` AND f4.`modelID` = 16 and f4.unit = ge.fromunit '...
	'INNER JOIN fits f5 ON f5.`nev file` = al.`manualrecording` AND f5.`modelID` = 19 and f5.unit = ge.fromunit '...
	'INNER JOIN fits_linear fl1 ON f3.`id` = fl1.id '...
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
tuningsize1 = [];
tuningsize2 = [];

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
		tuningsize1(idx) = data{i, 15};
		tuningsize2(idx) = data{i+1, 15};
	end
end

%%%%%%%%%%%
%MSE const%
%%%%%%%%%%%

figure
colormap(parula)
x = gc1; y = gc2; c = performance;
scatter(x,y,[],c, 'filled');
hold on
flm = fitlm([x',y'], c')
%plot([0 .4], [0 .4], 'k')
ylim([0 1.2])
xlim([0 1.2])
xlabel('Granger unit 1')
ylabel('Granger unit 2')
title('Performance (successes/sec; dual control)')
%caxis([0 .4])
colorbar
saveplot(gcf, './figures/grangerVsPerformance_dual.eps')

figure
colormap(jet)
x = mseconst1; y = mseconst2; c = gcnorm; a = 100*performance+20;
%scatter(x,y,[],c, 'filled');
scatter(y-x,c, 'filled')
%Fit line
f = fit((y-x)',c','poly1')
r2 = corr((y-x)', c')^2
hold on
plot(f)
xlabel('MSE unit 2 - MSE unit 1')
ylabel('\Delta Granger (Brain control)')
title(['R^2: ' num2str(r2)])
%caxis([0 .4])
%colorbar
saveplot(gcf, ['./figures/Const_mseoutVsGranger_dual_trendline.eps'])


figure
colormap(parula)
x = mseconst1; y = mseconst2; c = performance; a = 100*performance+20;
z = performance;
save('./figures/varianceVsPerformance.mat', 'mseconst1', 'mseconst2', 'performance')
scatter(x,y,[],c, 'filled');
hold on
%plot([0 .4], [0 .4], 'k')
ylim([0 2])
xlim([0 2])
xlabel('MSE out Unit 1')
ylabel('MSE out Unit 2')
title('Performance (successes/sec)')
%caxis([0 .4])
colorbar
saveplot(gcf, './figures/varianceVsPerformance_dual.eps')
