modelID = 34;
basename = './worksheets/2016_02_07-GPFApaired/run';

conn = db_conn();
paramcode = exec(conn, ['SELECT `description` FROM models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};
eval(paramcode);

%Extract the runs
conn = db_conn();
results = exec(conn, ['SELECT a.id, et.manualrecording, ' ...
' et.1DBCrecording, et.dualrecording, et.manualrecordingafter FROM analyses a ' ...
'INNER JOIN experiment_tuning et ' ...
'ON et.experiment_id = a.experiment_id WHERE modelID = ' num2str(modelID)]);
results = fetch(results);
results = results.Data;
nR = size(results,1);

ndim = 5;

paMCBC = [];
paMCDC = [];
paMCMC2 = [];
paBCDC = [];

paMCBCI = [];
paBCBCI = [];
paDCBCI = [];

pvalMCBC = [];
pvalMCDC = [];
pvalMCMC2 = [];
pvalBCDC = [];

nsigMCBC = [];
nsigMCDC = [];
nsigMCMC2 = [];
nsigBCDC = [];

pvalMCBCI = [];
pvalBCBCI = [];
pvalDCBCI = [];

allrates = [];

ratesMCBC = [];
ratesMCDC = [];
ratesBCDC = [];
ratesMCMC2 = [];

mu_dtheta_MCBC = [];
mu_dtheta_MCDC = [];
mu_dtheta_BCDC = [];
mu_dtheta_MCMC2 = [];

tuning_type = [];

bci_data = fetch(exec(conn, ['SELECT fl1.dir, fl2.dir, fl3.dir, fl5.dir, et1.`tuning_type`, rec1.`successrate`, rec2.`successrate`, '...
' IF(EXISTS (SELECT * FROM `bci_units` bci WHERE bci.`ID` = et1.`1DBCrecording` AND bci.unit = flin1.unit), '...
' 1, 0), flin1.`nev file` FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`manualrecording`'...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'INNER JOIN `fits` flin2 '...
'ON flin2.`nev file` = et1.`1DBCrecording`'...
'INNER JOIN `fits_linear` fl2 '...
'ON flin2.id = fl2.id '...
'INNER JOIN `fits` flin3 '...
'ON flin3.`nev file` = et1.`manualrecordingafter`'...
'INNER JOIN `fits_linear` fl3 '...
'ON flin3.id = fl3.id '...
'INNER JOIN `fits` flin5 '...
'ON flin5.`nev file` = et1.`dualrecording`'...
'INNER JOIN `fits_linear` fl5 '...
'ON flin5.id = fl5.id '...
'INNER JOIN `recordings` rec1 '...
'ON rec1.`nev file` = et1.`1DBCrecording` '...
'INNER JOIN `recordings` rec2 '...
'ON rec2.`nev file` = et1.`dualrecording` '...
'WHERE flin1.modelID = 30 AND flin2.modelID = 30 AND flin3.modelID = 30 AND flin5.modelID = 30 ' ...
'AND flin1.unit = flin2.unit AND flin2.unit = flin3.unit AND flin2.unit = flin5.unit '...
'AND fl1.r2 > .01 AND fl3.r2 > .01']));
all_r2 = cell2mat(bci_data.Data(:,1:4));
bcituningtype = cell2mat(bci_data.Data(:,5));
bciperformance = 50*cell2mat(bci_data.Data(:,6:7))+10;
bciunit = cell2mat(bci_data.Data(:,8));
nevfiles = bci_data.Data(:,9);

rot = (bcituningtype == 1 | bcituningtype == 3| bcituningtype == 4);
bci = bciunit == 1;
nonbci = bciunit == 0;

dof_bci = sum(bci)-1;
dof_nonbci = sum(nonbci)-1;

%1 = mc1 
%2 = bc
%3 = mc2 
%4 = dc 

mc1bc_difftheta_all = abs(center_angles(all_r2(:,1), all_r2(:,2)));
mc1mc2_difftheta_all = abs(center_angles(all_r2(:,1), all_r2(:,3)));
mc2bc_difftheta_all = abs(center_angles(all_r2(:,3), all_r2(:,2)));
mcdc_difftheta_all = abs(center_angles(all_r2(:,1), all_r2(:,4)));
bcdc_difftheta_all = abs(center_angles(all_r2(:,2), all_r2(:,4)));

mc1bc_difftheta = abs(center_angles(all_r2(rot,1), all_r2(rot,2)));
mc1mc2_difftheta = abs(center_angles(all_r2(rot,1), all_r2(rot,3)));
mc2bc_difftheta = abs(center_angles(all_r2(rot,3), all_r2(rot,2)));
mcdc_difftheta = abs(center_angles(all_r2(rot,1), all_r2(rot,4)));
bcdc_difftheta = abs(center_angles(all_r2(rot,2), all_r2(rot,4)));

bci_mc1bc_difftheta = abs(center_angles(all_r2(bci,1), all_r2(bci,2)));
bci_mc1mc2_difftheta = abs(center_angles(all_r2(bci,1), all_r2(bci,3)));
bci_mc2bc_difftheta = abs(center_angles(all_r2(bci,3), all_r2(bci,2)));
bci_mcdc_difftheta = abs(center_angles(all_r2(bci,1), all_r2(bci,4)));
bci_bcdc_difftheta = abs(center_angles(all_r2(bci,2), all_r2(bci,4)));

nonbci_mc1bc_difftheta = abs(center_angles(all_r2(nonbci,1), all_r2(nonbci,2)));
nonbci_mc1mc2_difftheta = abs(center_angles(all_r2(nonbci,1), all_r2(nonbci,3)));
nonbci_mc2bc_difftheta = abs(center_angles(all_r2(nonbci,3), all_r2(nonbci,2)));
nonbci_mcdc_difftheta = abs(center_angles(all_r2(nonbci,1), all_r2(nonbci,4)));
nonbci_bcdc_difftheta = abs(center_angles(all_r2(nonbci,2), all_r2(nonbci,4)));

%Open the file for each condition
for idx = 1:nR
	analysis_id = results{idx,1};
	mcnev = results{idx,2};
	bcnev = results{idx,3};
	dcnev = results{idx,4};
	mc2nev = results{idx,5};

	fn_out1 = [basename '_' num2str(analysis_id) '_' num2str(1) '.mat'];
	fn_out2 = [basename '_' num2str(analysis_id) '_' num2str(2) '.mat'];
	fn_out3 = [basename '_' num2str(analysis_id) '_' num2str(3) '.mat'];
	fn_out4 = [basename '_' num2str(analysis_id) '_' num2str(4) '.mat'];

	if exist(fn_out1, 'file')
		load(fn_out1);
		MCspace = estParams.Corth(:,1:ndim);
		load(fn_out2);
		BCspace = estParams.Corth(:,1:ndim);
		load(fn_out3);
		DCspace = estParams.Corth(:,1:ndim);
		load(fn_out4);
		MC2space = estParams.Corth(:,1:ndim);
	else 
		display(['Cant find all files. Continuing'])
		continue 
	end

	mcidx = find(cellfun(@(x) strcmp(mcnev, x), nevfiles));
	%Average change in tuning angle between conditions
	mu_dtheta_MCBC(idx) = mean(mc1bc_difftheta_all(mcidx));
	mu_dtheta_MCDC(idx) = mean(mcdc_difftheta_all(mcidx));
	mu_dtheta_MCMC2(idx) = mean(mc1mc2_difftheta_all(mcidx));
	mu_dtheta_BCDC(idx) = mean(bcdc_difftheta_all(mcidx));

	tuning_type(idx) = bcituningtype(mcidx(1));

	nU = size(MCspace, 1);

	%Figure out units to use... above 5hz in all recordings
	units = fetch(exec(conn, ['SELECT u1.unit FROM '...
	'`experiment_tuning` et1 '...
	'INNER JOIN `units` u1 '...
	'ON u1.`nev file` = et1.`manualrecording` '...
	'INNER JOIN `units` u2 '...
	'ON u2.`nev file` = et1.`1DBCrecording` '...
	'INNER JOIN `units` u3 '...
	'ON u3.`nev file` = et1.`manualrecordingafter` '...
	'INNER JOIN `units` u4 '...
	'ON u4.`nev file` = et1.`dualrecording` '...
	'WHERE u1.unit = u2.unit AND u1.unit = u4.unit AND u1.unit = u3.unit AND '...
	'u3.firingrate > ' num2str(threshold) ' AND ' ...
	'u1.firingrate > ' num2str(threshold) ' AND u2.firingrate > ' num2str(threshold) ' AND '...
	'u4.firingrate > ' num2str(threshold) ' AND et1.`manualrecording` = "' mcnev '"']));
	otherunits = units.Data;

	%Get the firing rate for the units
	rates = fetch(exec(conn, ['SELECT u1.firingrate, u2.firingrate, u4.firingrate, u3.firingrate FROM '...
	'`experiment_tuning` et1 '...
	'INNER JOIN `units` u1 '...
	'ON u1.`nev file` = et1.`manualrecording` '...
	'INNER JOIN `units` u2 '...
	'ON u2.`nev file` = et1.`1DBCrecording` '...
	'INNER JOIN `units` u3 '...
	'ON u3.`nev file` = et1.`manualrecordingafter` '...
	'INNER JOIN `units` u4 '...
	'ON u4.`nev file` = et1.`dualrecording` '...
	'WHERE u1.unit = u2.unit AND u1.unit = u4.unit AND u1.unit = u3.unit AND '...
	'u3.firingrate > ' num2str(threshold) ' AND ' ...
	'u1.firingrate > ' num2str(threshold) ' AND u2.firingrate > ' num2str(threshold) ' AND '...
	'u4.firingrate > ' num2str(threshold) ' AND et1.`manualrecording` = "' mcnev '"']));
	rates = cell2mat(rates.Data);

	ratesMCBC(idx) = corr(rates(:,1), rates(:,2));
	ratesMCDC(idx) = corr(rates(:,1), rates(:,3));
	ratesBCDC(idx) = corr(rates(:,2), rates(:,3));
	ratesMCMC2(idx) = corr(rates(:,1), rates(:,4));
	allrates = [allrates(:,:); rates];

	%Get brain control units for BC session
	bciunits = exec(conn, ['SELECT `unit` FROM bci_units WHERE `ID` = "' bcnev '"']);
	bciunits = fetch(bciunits);
	bciunits = bciunits.Data;
	%Get brain control units for DC session
	dualunits = exec(conn, ['SELECT `unit` FROM bci_units WHERE `ID` = "' dcnev '"']);
	dualunits = fetch(dualunits);
	dualunits = dualunits.Data;
	allunits = unique([otherunits; bciunits; dualunits]);
	a = sort(cellfun(@(x) str2num(x), allunits));
	nU = size(allunits,1);
	for j = 1:nU
			allunits{j} = num2str(a(j));
	end	

	%Get BCI space for BC/DC session
	BCIspace = zeros(nU, 2);
	ind1=find(ismember(allunits,bciunits{1}));
	ind2=find(ismember(allunits,bciunits{2}));
	BCIspace(ind1,1) = 1;
	BCIspace(ind2,2) = 1;

	%Compute the angles between each condition
	pasMCBC = subspacea(MCspace, BCspace);
	pasMCDC = subspacea(MCspace, DCspace);
	pasMCMC2 = subspacea(MCspace, MC2space);
	pasBCDC = subspacea(BCspace, DCspace);

	pvalsMCBC = 1-normcdf(sqrt(nU)*cos(pasMCBC));
	pvalsMCDC = 1-normcdf(sqrt(nU)*cos(pasMCDC));
	pvalsMCMC2 = 1-normcdf(sqrt(nU)*cos(pasMCMC2));
	pvalsBCDC = 1-normcdf(sqrt(nU)*cos(pasBCDC));

	nsigMCBC(idx) = sum(pvalsMCBC < 0.05);
	nsigMCDC(idx) = sum(pvalsMCDC < 0.05);
	nsigMCMC2(idx) = sum(pvalsMCMC2 < 0.05);
	nsigBCDC(idx) = sum(pvalsBCDC < 0.05);

	paMCBC(idx) = max(subspacea(MCspace, BCspace));
	paMCDC(idx) = max(subspacea(MCspace, DCspace));
	paMCMC2(idx) = max(subspacea(MCspace, MC2space));
	paBCDC(idx) = max(subspacea(BCspace, DCspace));

	pvalMCBC(idx) = 1-normcdf(sqrt(nU)*cos(paMCBC(idx)));
	pvalMCDC(idx) = 1-normcdf(sqrt(nU)*cos(paMCDC(idx)));
	pvalMCMC2(idx) = 1-normcdf(sqrt(nU)*cos(paMCMC2(idx)));
	pvalBCDC(idx) = 1-normcdf(sqrt(nU)*cos(paBCDC(idx)));

	%Compute the angles between neural space and the BC axes
	paMCBCI(idx) = max(subspacea(MCspace, BCIspace));
	paBCBCI(idx) = max(subspacea(BCspace, BCIspace));
	paDCBCI(idx) = max(subspacea(DCspace, BCIspace));

	pasMCBCI = subspacea(MCspace, BCIspace);
	pasBCBCI = subspacea(BCspace, BCIspace);
	pasDCBCI = subspacea(DCspace, BCIspace);

	pvalsMCBCI = 1-normcdf(sqrt(nU)*cos(pasMCBCI));
	pvalsBCBCI = 1-normcdf(sqrt(nU)*cos(pasBCBCI));
	pvalsDCBCI = 1-normcdf(sqrt(nU)*cos(pasDCBCI));

	nsigMCBCI(idx) = sum(pvalsMCBCI < 0.05);
	nsigBCBCI(idx) = sum(pvalsBCBCI < 0.05);
	nsigDCBCI(idx) = sum(pvalsDCBCI < 0.05);

	pvalMCBCI(idx) = 1-normcdf(sqrt(nU)*cos(paMCBCI(idx)));
	pvalBCBCI(idx) = 1-normcdf(sqrt(nU)*cos(paBCBCI(idx)));
	pvalDCBCI(idx) = 1-normcdf(sqrt(nU)*cos(paDCBCI(idx)));

	%Get whether was rotated or not...


	%Get peformance in BC and DC trials
	p = exec(conn, ['SELECT rec.successrate FROM `recordings` rec WHERE '...
		'rec.`nev file` = "' bcnev '"']);
	p = fetch(p);
	p = p.Data;
	BCperf(idx) = p{1};
	p = exec(conn, ['SELECT rec.successrate FROM `recordings` rec WHERE '...
		'rec.`nev file` = "' dcnev '"']);
	p = fetch(p);
	p = p.Data;
	DCperf(idx) = p{1};
end


%Scatter plot of angle between conditions over the days of recordings

%Scatter plot of angle between BCI axes and neural space (for all conditions) vs performance
figure 
subplot(2,2,1)
hold on 
scatter(log10(pvalMCBCI), BCperf)
plot([])
xlabel('pval BCI MC')
ylabel('BC success/sec')
subplot(2,2,2)
hold on 
scatter(log10(pvalBCBCI), BCperf)
xlabel('pval BCI BC')
ylabel('BC success/sec')
subplot(2,2,3)
hold on 
scatter(log10(pvalDCBCI), BCperf)
xlabel('pval BCI DC')
ylabel('BC success/sec')

figure 
subplot(2,2,1)
scatter3(pvalMCBCI, BCperf, DCperf)
xlabel('principal angle BCI MC')
ylabel('success/sec')
subplot(2,2,2)
scatter3(pvalBCBCI, BCperf, DCperf)
xlabel('principal angle BCI BC')
ylabel('success/sec')
subplot(2,2,3)
scatter3(pvalDCBCI, BCperf, DCperf)
xlabel('principal angle BCI DC')
ylabel('success/sec')

figure 
subplot(2,2,1)
scatter(DCperf, log10(pvalMCBCI))
ylabel('pval BCI MC')
xlabel('DC success/sec')
subplot(2,2,2)
scatter(DCperf, log10(pvalBCBCI))
ylabel('pval BCI BC')
xlabel('DC success/sec')
subplot(2,2,3)
scatter(DCperf, log10(pvalDCBCI))
ylabel('pval BCI DC')
xlabel('DC success/sec')

%%%%%%%%%%%%%%5
figure 
subplot(2,2,1)
scatter(paMCMC2, BCperf)
subplot(2,2,2)
scatter(paMCBC, BCperf)
subplot(2,2,3)
scatter(paMCDC, BCperf)
subplot(2,2,4)
scatter(paBCDC, BCperf)

figure 
subplot(2,2,1)
scatter(paMCMC2, DCperf)
subplot(2,2,2)
scatter(paMCBC, DCperf)
subplot(2,2,3)
scatter(paMCDC, DCperf)
subplot(2,2,4)
scatter(paBCDC, DCperf)

%Bin them into bins
nbins = 6;
MCDCperf = {};
for i = 1:nbins 
	MCDCperf{i} = [];
end

rmax = pi/2;
rmin = 0;
rnge = rmax - rmin;
for i = 1:length(DCperf)
	b = ceil(nbins*(paMCDC(i)-rmin)/rnge);
	b = max(b, 1);
	a = MCDCperf{b};
	MCDCperf{b} = [a; DCperf(i)];
end

mu_MCDCperf = zeros(1,nbins);
std_MCDCperf = zeros(1,nbins);
for i = 1:nbins
	m = MCDCperf{i};
	if length(m) > 0
		mu_MCDCperf(i) = mean(MCDCperf{i});
		std_MCDCperf(i) = std(MCDCperf{i});
	end
end

clf
hold on 
plot(mu_MCDCperf)
plot(mu_MCDCperf+std_MCDCperf)
plot(mu_MCDCperf-std_MCDCperf)

%Bin them into bins
nbins = 6;
MCBCperf = {};
for i = 1:nbins 
	MCBCperf{i} = [];
end

rmax = pi/2;
rmin = 0;
rnge = rmax - rmin;
for i = 1:length(BCperf)
	b = ceil(nbins*(paMCBC(i)-rmin)/rnge);
	b = max(b, 1);
	a = MCBCperf{b};
	MCBCperf{b} = [a; BCperf(i)];
end

mu_MCBCperf = zeros(1,nbins);
std_MCBCperf = zeros(1,nbins);
for i = 1:nbins
	m = MCBCperf{i};
	if length(m) > 0
		mu_MCBCperf(i) = mean(MCBCperf{i});
		std_MCBCperf(i) = std(MCBCperf{i});
	end
end

clf
hold on 
plot(mu_MCBCperf)
plot(mu_MCBCperf+std_MCBCperf)
plot(mu_MCBCperf-std_MCBCperf)

DCperfnonorth = DCperf(paMCBCI < 1.3);
DCperforth = DCperf(paMCBCI > 1.3);
[h1, p1] = ttest2(DCperfnonorth, DCperforth);

BCperfnonorth = BCperf(paMCBCI < 1.3);
BCperforth = BCperf(paMCBCI > 1.3);
[h2, p2] = ttest2(BCperfnonorth, BCperforth);

clf
hold on 
bar([60*mean(BCperforth), 60*mean(BCperfnonorth)])
errorbar([60*mean(BCperforth), 60*mean(BCperfnonorth)], [60*std(BCperforth), 60*std(BCperfnonorth)])
ylabel('performance (successes/minute)')
xlabel('orthogonal   non-orthogonal')
title(['2 sided t test p-val' num2str(p2)])
ylim([0 15])
saveplot(gcf, './worksheets/2016_02_07-GPFApaired/analyseGPFA_boxplotBCperf_paMCBCI.eps')
clf
hold on 
bar([60*mean(DCperforth), 60*mean(DCperfnonorth)])
errorbar([60*mean(DCperforth), 60*mean(DCperfnonorth)], [60*std(DCperforth), 60*std(DCperfnonorth)])
ylabel('performance (successes/minute)')
xlabel('orthogonal   non-orthogonal')
title(['2 sided t test p-val' num2str(p1)])
ylim([0 15])
saveplot(gcf, './worksheets/2016_02_07-GPFApaired/analyseGPFA_boxplotDCperf_paMCBCI.eps')

figure 
subplot(2,2,1)
scatter(paMCMC2, paMCBCI)
subplot(2,2,2)
scatter(paMCBC, paMCBCI)
subplot(2,2,3)
scatter(paMCDC, paMCBCI)
subplot(2,2,4)
scatter(paBCDC, paMCBCI)

figure 
subplot(2,2,1)
scatter(paMCMC2, paMCBCI)
subplot(2,2,2)
scatter(paMCBC, paMCBCI)
subplot(2,2,3)
scatter(paMCDC, paMCBCI)
subplot(2,2,4)
scatter(paBCDC, paMCBCI)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Num of significant non-orthogonal dims v performance

figure 
subplot(2,2,1)
scatter(nsigMCMC2, BCperf)
subplot(2,2,2)
scatter(nsigMCBC, BCperf)
subplot(2,2,3)
scatter(nsigMCDC, BCperf)
subplot(2,2,4)
scatter(nsigBCDC, BCperf)

figure 
subplot(2,2,1)
scatter(nsigMCMC2, DCperf)
subplot(2,2,2)
scatter(nsigMCBC, DCperf)
subplot(2,2,3)
scatter(nsigMCDC, DCperf)
subplot(2,2,4)
scatter(nsigBCDC, DCperf)

figure 
subplot(2,2,1)
boxplot(BCperf, nsigMCMC2)
subplot(2,2,2)
boxplot(BCperf, nsigMCBC)
subplot(2,2,3)
boxplot(BCperf, nsigMCDC)
subplot(2,2,4)
boxplot(BCperf, nsigBCDC)

figure 
subplot(2,2,1)
boxplot(DCperf, nsigMCMC2)
subplot(2,2,2)
boxplot(DCperf, nsigMCBC)
subplot(2,2,3)
boxplot(DCperf, nsigMCDC)
subplot(2,2,4)
boxplot(DCperf, nsigBCDC)

figure 
subplot(2,2,1)
boxplot(DCperf, nsigMCMC2)
subplot(2,2,2)
boxplot(DCperf, nsigMCBC)
subplot(2,2,3)
boxplot(DCperf, nsigMCDC)
subplot(2,2,4)
boxplot(DCperf, nsigDCBC)

figure 
subplot(2,2,1)
scatter(nsigMCBCI, BCperf)
subplot(2,2,2)
scatter(nsigBCBCI, BCperf)
subplot(2,2,3)
scatter(nsigDCBCI, BCperf)

figure 
subplot(2,2,1)
scatter(nsigMCBCI, DCperf)
subplot(2,2,2)
scatter(nsigBCBCI, DCperf)
subplot(2,2,3)
scatter(nsigDCBCI, DCperf)

figure 
subplot(2,2,1)
boxplot(BCperf, nsigMCBCI)
subplot(2,2,2)
boxplot(BCperf, nsigBCBCI)
subplot(2,2,3)
boxplot(BCperf, nsigDCBCI)

figure 
subplot(2,2,1)
boxplot(DCperf, nsigMCBCI)
subplot(2,2,2)
boxplot(DCperf, nsigBCBCI)
subplot(2,2,3)
boxplot(DCperf, nsigDCBCI)

figure
boxplot(DCperf*60, nsigMCBCI)
xlabel('Number non-orthogonal axes')
ylabel('Performance (successes/minute)')
saveplot(gcf, './worksheets/2016_02_07-GPFApaired/analyseGPFA_boxplotDCperf_nsigMCBCI.eps', 'eps', [3 3])

figure
boxplot(BCperf*60, nsigMCBCI)
xlabel('Number non-orthogonal axes')
ylabel('Performance (successes/minute)')
saveplot(gcf, './worksheets/2016_02_07-GPFApaired/analyseGPFA_boxplotBCperf_nsigMCBCI.eps', 'eps', [3 3])

%Figure out the proportion of DC performances above 5/min (5/60 = .0833/s)
DCperfabove0 = DCperf(DCperf > 0/12 & DCperf < 1/12);
nsigMCBCIabove0 = nsigMCBCI(DCperf > 0/12 & DCperf < 1/12);
DCperfabove10 = DCperf(DCperf > 1/12 & DCperf < 1/6);
nsigMCBCIabove10 = nsigMCBCI(DCperf > 1/12 & DCperf < 1/6);
DCperfabove20 = DCperf(DCperf > 1/6 & DCperf < 1/3);
nsigMCBCIabove20 = nsigMCBCI(DCperf > 1/6 & DCperf < 1/3);
DCperfabove30 = DCperf(DCperf > 1/3);
nsigMCBCIabove30 = nsigMCBCI(DCperf > 1/3);

tot0 = length(nsigMCBCIabove0);
tot10 = length(nsigMCBCIabove10);
tot20 = length(nsigMCBCIabove20);
tot30 = length(nsigMCBCIabove30);

data = [sum(nsigMCBCIabove0 == 0)/tot0,  sum(nsigMCBCIabove0 == 1)/tot0, sum(nsigMCBCIabove0 == 2)/tot0;
sum(nsigMCBCIabove10 == 0)/tot10, sum(nsigMCBCIabove10 == 1)/tot10, sum(nsigMCBCIabove10 == 2)/tot10;
sum(nsigMCBCIabove20 == 0)/tot20, sum(nsigMCBCIabove20 == 1)/tot20, sum(nsigMCBCIabove20 == 2)/tot20;
sum(nsigMCBCIabove30 == 0)/tot30, sum(nsigMCBCIabove30 == 1)/tot30, sum(nsigMCBCIabove30 == 2)/tot30];

area(data);

BCperfabove0 = BCperf(BCperf > 0/12 & BCperf < 1/12);
nsigMCBCIabove0 = nsigMCBCI(BCperf > 0/12 & BCperf < 1/12);
BCperfabove10 = BCperf(BCperf > 1/12 & BCperf < 1/6);
nsigMCBCIabove10 = nsigMCBCI(BCperf > 1/12 & BCperf < 1/6);
BCperfabove20 = BCperf(BCperf > 1/6 & BCperf < 1/3);
nsigMCBCIabove20 = nsigMCBCI(BCperf > 1/6 & BCperf < 1/3);
BCperfabove30 = BCperf(BCperf > 1/3);
nsigMCBCIabove30 = nsigMCBCI(BCperf > 1/3);

tot0 = length(nsigMCBCIabove0);
tot10 = length(nsigMCBCIabove10);
tot20 = length(nsigMCBCIabove20);
tot30 = length(nsigMCBCIabove30);

data = [sum(nsigMCBCIabove0 == 0)/tot0,  sum(nsigMCBCIabove0 == 1)/tot0, sum(nsigMCBCIabove0 == 2)/tot0;
sum(nsigMCBCIabove10 == 0)/tot10, sum(nsigMCBCIabove10 == 1)/tot10, sum(nsigMCBCIabove10 == 2)/tot10;
sum(nsigMCBCIabove20 == 0)/tot20, sum(nsigMCBCIabove20 == 1)/tot20, sum(nsigMCBCIabove20 == 2)/tot20;
sum(nsigMCBCIabove30 == 0)/tot30, sum(nsigMCBCIabove30 == 1)/tot30, sum(nsigMCBCIabove30 == 2)/tot30];
figure 
area(data);

figure 
hist(nsigBCDC)

%Plots of firing rates and performance:
%Don't see anything striking between firing rate and performance
%Probably have to look at variance information if want to see some relation
%between performance
figure 
subplot(2,2,1)
scatter(ratesMCBC, DCperf)
subplot(2,2,2)
scatter(ratesMCDC, DCperf)
subplot(2,2,3)
scatter(ratesBCDC, DCperf)
subplot(2,2,4)
scatter(ratesMCMC2, DCperf)

figure 
subplot(2,2,1)
scatter(ratesMCBC, BCperf)
subplot(2,2,2)
scatter(ratesMCDC, BCperf)
subplot(2,2,3)
scatter(ratesBCDC, BCperf)
subplot(2,2,4)
scatter(ratesMCMC2, BCperf)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Average change in tuning angle vs principal angle between spaces

figure 
subplot(2,2,1)
scatter(paMCMC2, mu_dtheta_MCMC2)
xlabel('principal angle')
ylabel('average |\Delta \theta|')
subplot(2,2,2)
scatter(paMCBC, mu_dtheta_MCBC)
xlabel('principal angle')
ylabel('average |\Delta \theta|')
subplot(2,2,3)
scatter(paMCDC, mu_dtheta_MCDC)
xlabel('principal angle')
ylabel('average |\Delta \theta|')
subplot(2,2,4)
scatter(paBCDC, mu_dtheta_BCDC)
xlabel('principal angle')
ylabel('average |\Delta \theta|')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rotated vs unrotated 

rot = (tuning_type == 1 | tuning_type == 3| tuning_type == 4);
cc = ones(size(rot));
cc(rot == 1) = 2;
colors = [1 0 0; 0 0 1];
c = [];
for idx = 1:size(cc, 1);
	c(idx,1:3) = colors(cc(idx),:);
end

figure
subplot(2,2,1)
scatter(paMCBCI, BCperf, [], cc, '.')
xlabel('principal angle')
ylabel('performance (successes/sec)')
subplot(2,2,2)
scatter(paBCBCI, BCperf, [], cc, '.')
xlabel('principal angle')
ylabel('performance (successes/sec)')
subplot(2,2,3)
scatter(paDCBCI, BCperf, [], cc, '.')
xlabel('principal angle')
ylabel('performance (successes/sec)')
%subplot(2,2,4)
%scatter(paBCDC, DCperf, [], cc, '.')
%xlabel('principal angle')
%ylabel('average |\Delta \theta|')

DCperfnonorth = DCperf(paMCBCI < 1.3);
DCperforth = DCperf(paMCBCI > 1.3);
[h1, p1] = ttest2(DCperfnonorth, DCperforth);

BCperfnonorth = BCperf(paMCBCI < 1.3);
BCperforth = BCperf(paMCBCI > 1.3);
[h2, p2] = ttest2(BCperfnonorth, BCperforth);


DCperfnonorth_rot = DCperf(paMCBCI < 1.3 & rot);
DCperforth_rot = DCperf(paMCBCI > 1.3 & rot);
[h3, p3] = ttest2(DCperfnonorth_rot, DCperforth_rot);

DCperfnonorth_unrot = DCperf(paMCBCI < 1.3 & ~rot);
DCperforth_unrot = DCperf(paMCBCI > 1.3 & ~rot);
[h4, p4] = ttest2(DCperfnonorth_unrot, DCperforth_unrot);

BCperfnonorth_rot = BCperf(paMCBCI < 1.3 & rot);
BCperforth_rot = BCperf(paMCBCI > 1.3 & rot);
[h5, p5] = ttest2(BCperfnonorth_rot, BCperforth_rot);

BCperfnonorth_unrot = BCperf(paMCBCI < 1.3 & ~rot);
BCperforth_unrot = BCperf(paMCBCI > 1.3 & ~rot);
[h6, p6] = ttest2(BCperfnonorth_unrot, BCperforth_unrot);

clf
hold on 
bar([60*mean(BCperforth_rot), 60*mean(BCperfnonorth_rot)])
errorbar([60*mean(BCperforth_rot), 60*mean(BCperfnonorth_rot)], [60*std(BCperforth_rot), 60*std(BCperfnonorth_rot)])
ylabel('performance (successes/minute)')
xlabel('orthogonal   non-orthogonal')
title(['Rotated. 2 sided t test p-val' num2str(p2)])
ylim([0 15])
saveplot(gcf, './worksheets/2016_02_07-GPFApaired/analyseGPFA_boxplotBCperf_paMCBCI_rotated.eps')

%clf
figure 
hold on 
bar([60*mean(BCperforth_unrot), 60*mean(BCperfnonorth_unrot)])
errorbar([60*mean(BCperforth_unrot), 60*mean(BCperfnonorth_unrot)], [60*std(BCperforth_unrot), 60*std(BCperfnonorth_unrot)])
ylabel('performance (successes/minute)')
xlabel('orthogonal   non-orthogonal')
title(['Unrotated. 2 sided t test p-val' num2str(p2)])
ylim([0 15])
saveplot(gcf, './worksheets/2016_02_07-GPFApaired/analyseGPFA_boxplotBCperf_paMCBCI_unrotated.eps')

figure 
hold on 
bar([60*mean(DCperforth_rot), 60*mean(DCperfnonorth_rot)])
errorbar([60*mean(DCperforth_rot), 60*mean(DCperfnonorth_rot)], [60*std(DCperforth_rot), 60*std(DCperfnonorth_rot)])
ylabel('performance (successes/minute)')
xlabel('orthogonal   non-orthogonal')
title(['Rotated. 2 sided t test p-val' num2str(p1)])
ylim([0 15])
saveplot(gcf, './worksheets/2016_02_07-GPFApaired/analyseGPFA_boxplotDCperf_paMCBCI_rotated.eps')

figure 
hold on 
bar([60*mean(DCperforth_unrot), 60*mean(DCperfnonorth_unrot)])
errorbar([60*mean(DCperforth_unrot), 60*mean(DCperfnonorth_unrot)], [60*std(DCperforth_unrot), 60*std(DCperfnonorth_unrot)])
ylabel('performance (successes/minute)')
xlabel('orthogonal   non-orthogonal')
title(['Unrotated. 2 sided t test p-val' num2str(p1)])
ylim([0 15])
saveplot(gcf, './worksheets/2016_02_07-GPFApaired/analyseGPFA_boxplotDCperf_paMCBCI_unrotated.eps')
