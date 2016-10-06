conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%List of files
files = fetch(exec(conn, ['SELECT et.`manualrecording` FROM experiment_tuning et']));
files = files.Data;

%Torque tuning angle (velocity)
deltaBCI = [];
deltacotuned = [];
deltaother = [];

for idx = 1:length(files)
	mcfile = files{idx}
	all_data = fetch(exec(conn, ['SELECT fl1.dir, fl2.dir, fl5.dir, et1.`tuning_type` FROM '...
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
	'AND fl1.r2 > 0.01 '...
	'AND NOT EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et1.`1DBCrecording` AND bci.unit = flin1.unit) '...
	'AND et1.`manualrecording` = "' mcfile '" '...
	'AND et1.`tuning_type` = 5']));
	if strcmp(all_data.Data, 'No Data')
		continue 
	end
	if size(all_data.Data, 1) < 4 
		continue 
	end
	all_theta = cell2mat(all_data.Data(:,1:3));
	tuningtype = cell2mat(all_data.Data(:,4));
	
	%BCI unit
	all_data = fetch(exec(conn, ['SELECT fl1.dir, fl2.dir, fl5.dir FROM '...
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
	'AND et1.`tuning_type` = 5 '...
	'LIMIT 1']));
	bci_theta = cell2mat(all_data.Data(:,1:3));
	nU = size(all_theta,1);
	
	%Pick cotuned units to BCI units in MC
	diff_MC_theta = cos(bci_theta(1) - all_theta(:,1));
	
	%Pick top two
	[l, cotunedidx] = sort(diff_MC_theta, 'descend'); 
	cotunedidx = cotunedidx(1:2);
	bcicotuned = diff_MC_theta(cotunedidx);
	
	%Pick random two 
	otherunits = randsample(setdiff(1:nU, cotunedidx), 2);
	
	deltathetaBCIMCDC = (bci_theta(3)-bci_theta(1))*180/pi
	deltathetacotunedMCDC = (all_theta(cotunedidx,3) - all_theta(cotunedidx,1))*180/pi
	deltathetaotherMCDC = (all_theta(otherunits,3) - all_theta(otherunits,1))*180/pi

	deltaBCI = [deltaBCI(:); deltathetaBCIMCDC; deltathetaBCIMCDC];
	deltacotuned = [deltacotuned(:); deltathetacotunedMCDC];
	deltaother = [deltaother(:); deltathetaotherMCDC];	
	%pause 
end

deltaBCI = mod(deltaBCI, 2*pi);
deltacotuned = mod(deltacotuned, 2*pi);
deltaother = mod(deltaother, 2*pi);

figure
subplot(1,2,1)
%Linear regression
hold on
[f, gof] = fit(deltaBCI,deltacotuned,'poly1')
plot(f)%ylim([0 2])

%fit_deltacotuned = polyfit(deltaBCI, deltacotuned, 1)
%xf = 0:(pi/2):(2*pi);
%yf = polyval(fit_deltacotuned, xf)
scatter(deltaBCI, deltacotuned, '.r')
%hold on 
%plot(xf, yf, 'r', 'linewidth', 2)
xlim([0 2*pi])
ylim([0 2*pi])
xlabel('\Delta \theta BC unit')
ylabel('\Delta \theta Cotuned unit')

subplot(1,2,2)
scatter(deltaBCI, deltaother,'.b')

hold on
f = fit(deltaBCI,deltaother,'poly1')
plot(f)%ylim([0 2])

%fit_deltaother = polyfit(deltaBCI, deltaother, 1)
%yf = polyval(fit_deltaother, xf)
%hold on 
%plot(xf, yf, 'b', 'linewidth', 2)
xlim([0 2*pi])
ylim([0 2*pi])
xlabel('\Delta \theta BC unit')
ylabel('\Delta \theta Other unit')
saveplot(gcf, './worksheets/2015_11_03-bcireportplots/bcituningchanges_angle_trendline.eps', 'eps', [5 3])

%Bin into quadrants and plot trend
nBins = 4;
cts = {};
mean_deltacotuned = zeros(nBins,1);
std_deltacotuned = zeros(nBins,1);
for i = 1:nBins
	cts{i} = [];
end
for i = 1:length(deltaBCI)
	db = deltaBCI(i);
	b = ceil(db*nBins/(2*pi));
	cts{b} = [cts{b}, deltacotuned(i)]
end

%Compute mean and std 
for i = 1:nBins
	mean_deltacotuned(i) = mean(cts{i});
	std_deltacotuned(i) = std(cts{i});
end

nBins = 4;
cts = {};
mean_deltaother = zeros(nBins,1);
std_deltaother = zeros(nBins,1);
for i = 1:nBins
	cts{i} = [];
end
for i = 1:length(deltaBCI)
	db = deltaBCI(i);
	b = ceil(db*nBins/(2*pi));
	cts{b} = [cts{b}, deltaother(i)]
end

%Compute mean and std 
for i = 1:nBins
	mean_deltaother(i) = mean(cts{i});
	std_deltaother(i) = std(cts{i});
end

figure
subplot(1,2,1)
plot(0:90:270, mean_deltacotuned*180/pi)
ylim([0 360])
xlabel('\Delta \theta BC unit')
ylabel('\Delta \theta Cotuned unit')
subplot(1,2,2)
plot(0:90:270, mean_deltaother*180/pi)
ylim([0 360])
xlabel('\Delta \theta BC unit')
ylabel('\Delta \theta Other unit')