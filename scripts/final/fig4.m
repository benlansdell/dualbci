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
