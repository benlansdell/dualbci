
%%%%%%%%%%%%%%%%%%%%
%Summary MVGC stats%
%%%%%%%%%%%%%%%%%%%%

conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/Spanky');
data = exec(conn, ['SELECT ge.fromunit, ge.score, rec.successrate, f.`nev file`, bci.unit '...
	'FROM GrangerEstimates ge '...
	'INNER JOIN Fits f ON f.id = ge.id '...
	'INNER JOIN Recordings rec ON rec.`nev file` = f.`nev file`'...
    'LEFT JOIN BCIUnits bci ON bci.`id` = f.`nev file` AND bci.`unit` = ge.`fromunit`'...
    'WHERE f.`modelID` = 3']);
data = fetch(data);
data = data.Data;

%Histogram  
bcindices = ~cellfun(@(x) strcmp(x, 'null'), data(:,5)); 
mvgc = cell2mat(data(:,2));
clf
hist(log10(mvgc(bcindices)));
set(get(gca,'child'),'FaceColor',[0 0.5 0.5],'EdgeColor','w');
hold on
hist(log10(mvgc(~bcindices)));
legend('BCI units', 'Non-BCI units')
xlabel('log10(GC score)')
ylabel('frequency')
saveplot(gcf, './worksheets/2015_07_09/MVGCanalysis1.eps')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Proportion of GC score versus linear tuning%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = exec(conn, ['SELECT ge.fromunit, ge.score, rec.successrate, f.`nev file`, fl.dir, bci.angle, MOD(fl.dir-bci.angle, 2*PI()), fl.size, u.firingrate, f1.`mse out` '...
	'FROM GrangerEstimates ge '...
	'INNER JOIN Fits f ON f.id = ge.id '...
	'INNER JOIN Recordings rec ON rec.`nev file` = f.`nev file` '...
	'INNER JOIN BCIUnits bci ON bci.ID = f.`nev file` AND bci.unit = ge.fromunit '...
	'INNER JOIN AnalysisLinear al ON al.`1DBCrecording` = f.`nev file` '...
	'INNER JOIN Fits f1 ON f1.`nev file` = al.`manualrecording` AND f1.`modelID` = 1 and f1.unit = ge.fromunit '...
	'INNER JOIN FitsLinear fl ON f1.`id` = fl.id '...
	'INNER JOIN Units u ON u.`nev file` = f1.`nev file` AND u.`unit` = f1.unit WHERE f.modelID = 3 ']);

data = fetch(data);
data = data.Data;

gc = [];
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

for i = 1:(size(data,1)-1)
	if strcmp(data{i, 4}, data{i+1, 4})
		idx = length(gc)+1;
		%delangle(idx+1) = data{i+1,7};
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
		performance(idx) = data{i, 3};
	end
end
%Make between pi and -pi
delangle1 = mod(delangle1, 2*pi);
delangle1(delangle1>pi) = delangle1(delangle1>pi)-2*pi;
%Make between pi and -pi
delangle2 = mod(delangle2, 2*pi);
delangle2(delangle2>pi) = delangle2(delangle2>pi)-2*pi;

figure
x = delangle1; y = delangle2; c = gcnorm;
scatter(180*x/pi,180*y/pi,[],c, 'filled');
xlabel('Delta angle 1')
ylabel('Delta angle 2')
xlim([-180 180])
ylim([-180 180])
colorbar
%saveplot(gcf, './worksheets/2015_06_18/angleVsGranger_log.eps')
saveplot(gcf, './worksheets/2015_07_09/angleVsGranger.eps')

figure
x = tuningsize1; y = tuningsize2; c = gcnorm;
scatter(x,y,[],c, 'filled');
hold on
plot([0 4], [0 4])
xlabel('Tuning size1')
ylabel('Tuning size2')
colorbar
saveplot(gcf, './worksheets/2015_07_09/tuningsizeVsGranger.eps')

figure
x = mseout1; y = mseout2; c = gcnorm;
scatter(x,y,[],c, 'filled');
hold on
%plot([0 4], [0 4])
xlabel('MSE out of sample1')
ylabel('MSE out of sample2')
colorbar
saveplot(gcf, './worksheets/2015_07_09/mseoutVsGranger.eps')

figure
x = firing1; y = firing2; c = gcnorm;
%Looks like a slight pattern in this case, fit a linear model to confirm:
[bhat, dev, stats] = glmfit([x', y'], c, 'normal');
%Note in the stats struct that the coefficients are sig non-zero, and that they infact fit a plane 
%roughly oriented as we'd expect...!
[xx, yy] = meshgrid(0:50, 0:50);

scatter3(x,y,c, 'filled');
hold on
xlabel('Firing rate1')
ylabel('Firing rate2')
zlabel('Difference in Granger score')
plot3([0 50], [0 50], [0 0], 'r', 'LineWidth', 2)
h = surf(xx, yy, bhat(1)+bhat(2)*xx+bhat(3)*yy, 'FaceColor', 'none');
saveplot(gcf, './worksheets/2015_07_09/firingrateVsGranger3d.eps')

%Plot performance versus balance-ness:
figure
plot(abs(gcnorm), performance, '.');
xlabel('GC balance')
ylabel('Performance (successes/min)')
saveplot(gcf, './worksheets/2015_07_09/GCvsperformance.eps')

