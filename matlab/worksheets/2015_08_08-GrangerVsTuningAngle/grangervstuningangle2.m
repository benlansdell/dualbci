conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%%%%%%%%%%%%%%%%
%Manual control%
%%%%%%%%%%%%%%%%

data = fetch(exec(conn, ['SELECT fl.`dir`, fl2.`dir`, eg.`score` FROM fits f '...
' INNER JOIN fits_linear fl '...
' ON f.id = fl.id '...
' INNER JOIN fits f2 '...
' ON f.`nev file` = f2.`nev file` AND f.`unit` = f2.`unit`'...
' INNER JOIN estimates_granger eg'...
' ON f2.id = eg.id'...
' INNER JOIN fits f3 '...
' ON f3.`nev file` = f.`nev file` AND f3.`unit` = eg.fromunit'...
' INNER JOIN fits_linear fl2 '...
' ON fl2.id = f3.id '...
' INNER JOIN experiment_tuning et '...
' ON et.`manualrecording` = f.`nev file`'...
' WHERE f.modelID = 1 AND f2.modelID = 2 AND f3.modelID = 1']));

d = cell2mat(data.Data(:,1:3));

gcthresholds = 0.1:5:200;
anglethreshs = [10]/180*pi;
proptuned = zeros(length(anglethreshs),length(gcthresholds));
diffangle = mod(d(:,1)-d(:,2), 2*pi);
diffangle(diffangle < 0) = diffangle(diffangle < 0) + 2*pi;
diffangle(diffangle > pi) = diffangle(diffangle > pi) - 2*pi;
diffangle = abs(diffangle);

nabove = [];
nsimtuned = [];
stdproptuned = [];
sem =[];
corrangles = [];
overlap = [];
stdoverlap = [];

nC = 200;
N = size(d,1);
cmap = jet(nC);

nabovecmap = zeros(size(cmap));
for i = 1:length(anglethreshs)
	for idx = 1:length(gcthresholds)
		%Find those above threshold
		gcthresh = gcthresholds(idx);
		gcabove = d(:,3)>gcthresh;
		nabove(i,idx) = sum(gcabove);
		a1 = d(gcabove,1);
		a2 = d(gcabove,2);
		overlap(idx) = mean(cos(a1-a2));
		stdoverlap(idx) = std(cos(a1-a2))/sqrt(nabove(i,idx));
		corrangles(idx) = corr(a1,a2);
		ii = ceil(nC*nabove(i,idx)/N);
		nabovecmap(idx,:) = cmap(ii,:);
		nsimtuned(i,idx) = sum(diffangle(gcabove)<anglethreshs(i));
		%Find the prop of those whose difference in tuning angle is less than anglethresh
		p = nsimtuned(i,idx)/nabove(i,idx);
		proptuned(i,idx) = p;
		stdproptuned(i,idx) = sqrt(p*(1-p));
		sem(i,idx) = stdproptuned(i,idx)/sqrt(nabove(i,idx));
	end
end
colormap(nabovecmap);

nbetween = diff(nabove);
nsimtunedbetween = diff(nsimtuned);

figure
[gcthresholds; nabove; nsimtuned]
for i = 1:length(anglethreshs)
	ii = 1:find(nabove(i,:)<30, 1);
	hold on
	plot(gcthresholds(ii), proptuned(i,ii), gcthresholds(ii), proptuned(i,ii)+sem(i,ii), '--r', gcthresholds(ii), proptuned(i,ii)-sem(i,ii), '--r', gcthresholds, ones(size(gcthresholds))*(10/180), ':k')
end
colorbar('south')
xlabel('GC score')
ylabel('Proportion tuned within +/- 10 degrees')
title('Manual control. Ipsi+contra')
ylim([0 .4])
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/gcscorevproptuned_MC.eps')

ii = 1:find(nabove(i,:)<400, 1);
plot(gcthresholds(ii), corrangles(ii))
xlabel('GC score')
ylabel('Correlation between tuning angles')
title('Manual control. Ipsi+contra')
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/gcscorevtuningcorrelation_MC.eps')

ii = 1:find(nabove(i,:)<400, 1);
plot(gcthresholds(ii), overlap(ii), gcthresholds(ii), overlap(ii)-stdoverlap(ii), 'r--', gcthresholds(ii), overlap(ii)+stdoverlap(ii), 'r--')
xlabel('GC score')
ylabel('Overlap between tuning angles')
title('Manual control. Ipsi+contra')
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/gcscorevtuningoverlap_MC.eps')

%Note the non-uniform distribution of preferred angles (even with the linear model)
hist(d(:,2), 50) 

absangle = [];
stdabsangle = [];

figure
%cmap = jet(length(gcthresholds))
for idx = 1:length(gcthresholds)
%for idx = 1:20
	%Find those above threshold
	gcthresh = gcthresholds(idx);
	gcabove = d(:,3)>gcthresh;
	nabove(idx) = sum(gcabove);
	%absangle(idx) = mean(abs(diffangle(gcabove)));
	i90 = diffangle < pi/2;
	absangle(idx) = mean(abs(diffangle(gcabove&i90)));
	stdabsangle(idx) = std(abs(diffangle(gcabove)));
	%hist(d(gcabove,2), 50) 
	%title(['GC threshold: ' num2str(gcthresh)])
	%pause
	%[f,xi] = ksdensity(d(gcabove,2), 'width', 0.1);
	%hold on
	%plot(xi, f, 'Color', cmap(idx,:));
end

ii = 1:find(nabove<30, 1);
%plot(gcthresholds(ii), absangle(ii), gcthresholds(ii), absangle(ii)-stdabsangle(ii), gcthresholds(ii), absangle(ii)+stdabsangle(ii), gcthresholds, ones(size(gcthresholds))*(10/180), ':k')
plot(gcthresholds, absangle, gcthresholds, ones(size(gcthresholds))*(10/180), ':k')
xlabel('GC score')
ylabel('Mean |diff in angle|')
title('Manual control. Ipsi+contra')
ylim([0 pi])
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/gcscorevmeantuned_MC.eps')


%Plot on a log scale

gcthresholds = 0.1:.1:log(200);
anglethreshs = [10]/180*pi;
proptuned = zeros(length(anglethreshs),length(gcthresholds));
diffangle = mod(d(:,1)-d(:,2), 2*pi);
diffangle(diffangle < 0) = diffangle(diffangle < 0) + 2*pi;
diffangle(diffangle > pi) = diffangle(diffangle > pi) - 2*pi;
diffangle = abs(diffangle);

nabove = [];
nsimtuned = [];
stdproptuned = [];
sem =[];

nC = length(gcthresholds);
N = size(d,1);
cmap = bone(nC);

nabovecmap = zeros(length(gcthresholds),3);
for i = 1:length(anglethreshs)
	for idx = 1:length(gcthresholds)
		%Find those above threshold
		gcthresh = gcthresholds(idx);
		gcabove = log(d(:,3))>gcthresh;
		nabove(i,idx) = sum(gcabove);
		ii = ceil(nC*nabove(i,idx)/N);
		nabovecmap(idx,:) = cmap(ii,:);
		nsimtuned(i,idx) = sum(diffangle(gcabove)<anglethreshs(i));
		%Find the prop of those whose difference in tuning angle is less than anglethresh
		p = nsimtuned(i,idx)/nabove(i,idx);
		proptuned(i,idx) = p;
		stdproptuned(i,idx) = sqrt(p*(1-p));
		sem(i,idx) = stdproptuned(i,idx)/sqrt(nabove(i,idx));
	end
end

nbetween = diff(nabove);
nsimtunedbetween = diff(nsimtuned);

figure
[gcthresholds; nabove; nsimtuned]
for i = 1:length(anglethreshs)
	ii = 1:find(nabove(i,:)<30, 1);
	hold on
	plot(gcthresholds(ii), proptuned(i,ii), gcthresholds(ii), proptuned(i,ii)+sem(i,ii), '--r', gcthresholds(ii), proptuned(i,ii)-sem(i,ii), '--r', gcthresholds, ones(size(gcthresholds))*(10/180), ':k')
end
colorbar('southoutside')
colormap(cmap);
caxis([0, 100])
xlabel('GC score. Percent units above.')
ylabel('Proportion tuned within +/- 10 degrees')
title('Manual control. Ipsi+contra')
ylim([0 .4])
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/loggcscorevproptuned_MC.eps')


%%%%%%%%%%%%%%%
%Brain control%
%%%%%%%%%%%%%%%

data = fetch(exec(conn, ['SELECT fl.`dir`, fl2.`dir`, eg.`score` FROM fits f '...
' INNER JOIN fits_linear fl '...
' ON f.id = fl.id '...
' INNER JOIN fits f2 '...
' ON f.`nev file` = f2.`nev file` AND f.`unit` = f2.`unit`'...
' INNER JOIN estimates_granger eg'...
' ON f2.id = eg.id'...
' INNER JOIN fits f3 '...
' ON f3.`nev file` = f.`nev file` AND f3.`unit` = eg.fromunit'...
' INNER JOIN fits_linear fl2 '...
' ON fl2.id = f3.id '...
' INNER JOIN experiment_tuning et '...
' ON et.`1DBCrecording` = f.`nev file`'...
' WHERE f.modelID = 1 AND f2.modelID = 2 AND f3.modelID = 1']));

d = cell2mat(data.Data(:,1:3));

gcthresholds = 0:5:200;
proptuned = zeros(size(gcthresholds));
anglethresh = 10/180*pi;
diffangle = mod(d(:,1)-d(:,2), 2*pi);
diffangle(diffangle < 0) = diffangle(diffangle < 0) + 2*pi;
diffangle(diffangle > pi) = diffangle(diffangle > pi) - 2*pi;
diffangle = abs(diffangle);

nabove = [];
nsimtuned = [];

for idx = 1:length(gcthresholds)
	%Find those above threshold
	gcthresh = gcthresholds(idx);
	gcabove = d(:,3)>gcthresh;
	nabove(idx) = sum(gcabove);
	nsimtuned(idx) = sum(diffangle(gcabove)<anglethresh);
	%Find the prop of those whose difference in tuning angle is less than anglethresh
	proptuned(idx) = nsimtuned(idx)/nabove(idx);
end

[gcthresholds; nabove; nsimtuned]

figure
ii = 1:find(nabove<30, 1);
plot(gcthresholds(ii), proptuned(ii), gcthresholds, ones(size(gcthresholds))*(10/180), ':k')
xlabel('GC score threshold')
ylabel('Proportion tuned within +/- 10 degrees')
title('Brain control. Ipsi+contra')
ylim([0 .4])
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/gcscorevproptuned_BC_all.eps')

%Note the non-uniform distribution of preferred angles (even with the linear model)
hist(d(:,2), 50) 

absangle = [];
stdabsangle = [];

figure
cmap = jet(length(gcthresholds))
%for idx = 1:length(gcthresholds)
for idx = 1:20
	%Find those above threshold
	gcthresh = gcthresholds(idx);
	gcabove = d(:,3)>gcthresh;
	nabove(idx) = sum(gcabove);
	absangle(idx) = mean(abs(diffangle(gcabove)));
	stdabsangle(idx) = std(abs(diffangle(gcabove)));
	%hist(d(gcabove,2), 50) 
	%title(['GC threshold: ' num2str(gcthresh)])
	%pause
	[f,xi] = ksdensity(d(gcabove,2), 'width', 0.1);
	hold on
	plot(xi, f, 'Color', cmap(idx,:));
end

ii = 1:find(nabove<30, 1);
plot(gcthresholds(ii), absangle(ii), gcthresholds(ii), absangle(ii)-stdabsangle(ii), gcthresholds(ii), absangle(ii)+stdabsangle(ii), gcthresholds, ones(size(gcthresholds))*(10/180), ':k')
xlabel('GC score')
ylabel('Mean |diff in angle|')
title('Manual control. Ipsi+contra')
ylim([0 pi])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%All data, split into hemispheres%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Connections from ipsi- to contralateral

data = fetch(exec(conn, ['SELECT fl.`dir`, fl2.`dir`, eg.`score`, MOD(fl.dir- fl2.dir, 2*PI()) FROM fits f '...
' INNER JOIN fits_linear fl '...
' ON f.id = fl.id '...
' INNER JOIN fits f2 '...
' ON f.`nev file` = f2.`nev file` AND f.`unit` = f2.`unit`'...
' INNER JOIN estimates_granger eg'...
' ON f2.id = eg.id'...
' INNER JOIN fits f3 '...
' ON f3.`nev file` = f.`nev file` AND f3.`unit` = eg.fromunit'...
' INNER JOIN fits_linear fl2 '...
' ON fl2.id = f3.id '...
' WHERE f.modelID = 1 AND f2.modelID = 2 AND f3.modelID = 1 AND f.unit < 97 AND f3.unit > 97']));

%' INNER JOIN experiment_tuning et '...
%' ON et.`manualrecording` = f.`nev file`'...


d = cell2mat(data.Data(:,1:4));

nbins = 100;
bins = cell(nbins,1);
%Set scores to zero
for idx = 1:nbins
	bins{idx} = [];
end

for idx = 1:size(d,1)
	%Figure out which bin to put score in...
	angle = d(idx,4);
	if angle < 0 
		angle = angle + 2*pi;
	end
	if angle > pi
		angle = angle - 2*pi;
	end
	score = d(idx,3);
	bin = ceil(nbins*(angle+pi)/(2*pi));
	bins{bin} = [bins{bin},score];
end

nK_sp = 6;
sigthresh = 22.35; %Because: 1-chi2cdf(12.59, nK_sp) = 0.05;

%Average scores, and variance of scores
scores = [];
stdscores = [];
catted = [];
groups = [];
for idx = 1:nbins
	scores(idx) = mean(bins{idx});
	stdscores(idx) = std(bins{idx});
	propabove(idx) = sum(bins{idx}>sigthresh)/length(bins{idx});
	maxscores(idx) = max(bins{idx});
	catted = [catted, bins{idx}];
	groups = [groups, idx*ones(size(bins{idx}))];
end
angles = linspace(-pi, pi, nbins)/pi*180;
figure 
plot(angles, scores);
hold on
plot(angles, scores-stdscores, '--r')
plot(angles, scores+stdscores, '--r')
plot([-180, 180], [sigthresh, sigthresh], ':k')
xlabel('Difference in tuning angles between units')
ylabel('Granger score between units')
title('Brain+manual control. Ipsi->contra GC')
xlim([-180 180])
ylim([0 40])
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/GCvtuning_alldata_ipsi-contra.eps')

figure
plot(angles, propabove)
ylabel('Proportion of significant connections. \alpha = 0.001')
ylim([0 .4])
xlabel('Difference in tuning angles between units')
title('Brain+manual control. Ipsi->contra GC')
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/tuningvpropsig_alldata_ipsi-contra.eps')


gcthresholds = 0:5:200;
proptuned = zeros(size(gcthresholds));
anglethresh = 10/180*pi;
diffangle = mod(d(:,1)-d(:,2), 2*pi);
diffangle(diffangle < 0) = diffangle(diffangle < 0) + 2*pi;
diffangle(diffangle > pi) = diffangle(diffangle > pi) - 2*pi;
diffangle = abs(diffangle);

nabove = [];
nsimtuned = [];

for idx = 1:length(gcthresholds)
	%Find those above threshold
	gcthresh = gcthresholds(idx);
	gcabove = d(:,3)>gcthresh;
	nabove(idx) = sum(gcabove);
	nsimtuned(idx) = sum(diffangle(gcabove)<anglethresh);
	%Find the prop of those whose difference in tuning angle is less than anglethresh
	proptuned(idx) = nsimtuned(idx)/nabove(idx);
end

[gcthresholds; nabove; nsimtuned]

figure
ii = 1:find(nabove<30, 1);
plot(gcthresholds(ii), proptuned(ii), gcthresholds, ones(size(gcthresholds))*(10/180), ':k')
xlabel('GC score threshold')
ylabel('Proportion of similarly tuned')
title('Brain+manual. Ipsi->contra GC')
ylim([0 .4])
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/GCscorevproptuned_alldata_ipsi-contra.eps')

%Compare with intra-hemisphere connections

data = fetch(exec(conn, ['SELECT fl.`dir`, fl2.`dir`, eg.`score`, MOD(fl.dir- fl2.dir, 2*PI()) FROM fits f '...
' INNER JOIN fits_linear fl '...
' ON f.id = fl.id '...
' INNER JOIN fits f2 '...
' ON f.`nev file` = f2.`nev file` AND f.`unit` = f2.`unit`'...
' INNER JOIN estimates_granger eg'...
' ON f2.id = eg.id'...
' INNER JOIN fits f3 '...
' ON f3.`nev file` = f.`nev file` AND f3.`unit` = eg.fromunit'...
' INNER JOIN fits_linear fl2 '...
' ON fl2.id = f3.id '...
' WHERE f.modelID = 1 AND f2.modelID = 2 AND f3.modelID = 1 AND f.unit < 97 AND f3.unit < 97']));

%' INNER JOIN experiment_tuning et '...
%' ON et.`manualrecording` = f.`nev file`'...


d = cell2mat(data.Data(:,1:4));

nbins = 100;
bins = cell(nbins,1);
%Set scores to zero
for idx = 1:nbins
	bins{idx} = [];
end

for idx = 1:size(d,1)
	%Figure out which bin to put score in...
	angle = d(idx,4);
	if angle < 0 
		angle = angle + 2*pi;
	end
	if angle > pi
		angle = angle - 2*pi;
	end
	score = d(idx,3);
	bin = ceil(nbins*(angle+pi)/(2*pi));
	bins{bin} = [bins{bin},score];
end

nK_sp = 6;
sigthresh = 22.35; %Because: 1-chi2cdf(12.59, nK_sp) = 0.05;

%Average scores, and variance of scores
scores = [];
stdscores = [];
catted = [];
groups = [];
for idx = 1:nbins
	scores(idx) = mean(bins{idx});
	stdscores(idx) = std(bins{idx});
	propabove(idx) = sum(bins{idx}>sigthresh)/length(bins{idx});
	maxscores(idx) = max(bins{idx});
	catted = [catted, bins{idx}];
	groups = [groups, idx*ones(size(bins{idx}))];
end
angles = linspace(-pi, pi, nbins)/pi*180;
figure
plot(angles, scores);
hold on
plot(angles, scores-stdscores, '--r')
plot(angles, scores+stdscores, '--r')
plot([-180, 180], [sigthresh, sigthresh], ':k')
xlabel('Difference in tuning angles between units')
ylabel('Granger score between units')
title('Brain+manual control. Contra->contra GC')
xlim([-180 180])
ylim([0 40])

saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/GCvtuning_alldata_contra-contra.eps')

figure
plot(angles, propabove)
ylabel('Proportion of significant connections. \alpha = 0.001')
ylim([0 .4])
xlabel('Difference in tuning angles between units')
title('Brain+manual control. Contra->contra GC')
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/tuningvpropsig_alldata_contra-contra.eps')

gcthresholds = 0:5:200;
proptuned = zeros(size(gcthresholds));
anglethresh = 10/180*pi;
diffangle = mod(d(:,1)-d(:,2), 2*pi);
diffangle(diffangle < 0) = diffangle(diffangle < 0) + 2*pi;
diffangle(diffangle > pi) = diffangle(diffangle > pi) - 2*pi;
diffangle = abs(diffangle);

nabove = [];
nsimtuned = [];

for idx = 1:length(gcthresholds)
	%Find those above threshold
	gcthresh = gcthresholds(idx);
	gcabove = d(:,3)>gcthresh;
	nabove(idx) = sum(gcabove);
	nsimtuned(idx) = sum(diffangle(gcabove)<anglethresh);
	%Find the prop of those whose difference in tuning angle is less than anglethresh
	proptuned(idx) = nsimtuned(idx)/nabove(idx);
end

[gcthresholds; nabove; nsimtuned]

figure
ii = 1:find(nabove<30, 1);
plot(gcthresholds(ii), proptuned(ii), gcthresholds, ones(size(gcthresholds))*(10/180), ':k')
xlabel('GC score')
ylabel('Proportion of similarly tuned')
title('Brain+manual control. Contra->contra GC')
ylim([0 .4])
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/GCscorevproptuned_alldata_contra-contra.eps')

%Connections from contra- to ipsilateral

data = fetch(exec(conn, ['SELECT fl.`dir`, fl2.`dir`, eg.`score`, MOD(fl.dir- fl2.dir, 2*PI()) FROM fits f '...
' INNER JOIN fits_linear fl '...
' ON f.id = fl.id '...
' INNER JOIN fits f2 '...
' ON f.`nev file` = f2.`nev file` AND f.`unit` = f2.`unit`'...
' INNER JOIN estimates_granger eg'...
' ON f2.id = eg.id'...
' INNER JOIN fits f3 '...
' ON f3.`nev file` = f.`nev file` AND f3.`unit` = eg.fromunit'...
' INNER JOIN fits_linear fl2 '...
' ON fl2.id = f3.id '...
' WHERE f.modelID = 1 AND f2.modelID = 2 AND f3.modelID = 1 AND f.unit > 97 AND f3.unit < 97']));

%' INNER JOIN experiment_tuning et '...
%' ON et.`manualrecording` = f.`nev file`'...

d = cell2mat(data.Data(:,1:4));

nbins = 100;
bins = cell(nbins,1);
%Set scores to zero
for idx = 1:nbins
	bins{idx} = [];
end

for idx = 1:size(d,1)
	%Figure out which bin to put score in...
	angle = d(idx,4);
	if angle < 0 
		angle = angle + 2*pi;
	end
	if angle > pi
		angle = angle - 2*pi;
	end
	score = d(idx,3);
	bin = ceil(nbins*(angle+pi)/(2*pi));
	bins{bin} = [bins{bin},score];
end

nK_sp = 6;
sigthresh = 22.35; %Because: 1-chi2cdf(12.59, nK_sp) = 0.05;

%Average scores, and variance of scores
scores = [];
stdscores = [];
catted = [];
groups = [];
for idx = 1:nbins
	scores(idx) = mean(bins{idx});
	stdscores(idx) = std(bins{idx});
	propabove(idx) = sum(bins{idx}>sigthresh)/length(bins{idx});
	maxscores(idx) = max(bins{idx});
	catted = [catted, bins{idx}];
	groups = [groups, idx*ones(size(bins{idx}))];
end
angles = linspace(-pi, pi, nbins)/pi*180;
figure 
plot(angles, scores);
hold on
plot(angles, scores-stdscores, '--r')
plot(angles, scores+stdscores, '--r')
plot([-180, 180], [sigthresh, sigthresh], ':k')
xlabel('Difference in tuning angles between units')
ylabel('Granger score between units')
title('Brain+manual control. Contra->ipsi GC')
xlim([-180 180])
ylim([0 40])

saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/GCvtuning_alldata_contra-ipsi.eps')

figure
plot(angles, propabove)
ylabel('Proportion of significant connections. \alpha = 0.001')
ylim([0 .4])
xlabel('Difference in tuning angles between units')
title('Brain+manual control. Contra->ipsi GC')

saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/tuningvpropsig_alldata_contra-ipsi.eps')

gcthresholds = 0:5:200;
proptuned = zeros(size(gcthresholds));
anglethresh = 10/180*pi;
diffangle = mod(d(:,1)-d(:,2), 2*pi);
diffangle(diffangle < 0) = diffangle(diffangle < 0) + 2*pi;
diffangle(diffangle > pi) = diffangle(diffangle > pi) - 2*pi;
diffangle = abs(diffangle);

nabove = [];
nsimtuned = [];

for idx = 1:length(gcthresholds)
	%Find those above threshold
	gcthresh = gcthresholds(idx);
	gcabove = d(:,3)>gcthresh;
	nabove(idx) = sum(gcabove);
	nsimtuned(idx) = sum(diffangle(gcabove)<anglethresh);
	%Find the prop of those whose difference in tuning angle is less than anglethresh
	proptuned(idx) = nsimtuned(idx)/nabove(idx);
end

[gcthresholds; nabove; nsimtuned]

figure
ii = 1:find(nabove<30, 1);
plot(gcthresholds(ii), proptuned(ii), gcthresholds, ones(size(gcthresholds))*(10/180), ':k')
xlabel('GC score threshold')
ylabel('Proportion of similarly tuned')
title('Brain+manual control. Contra->ipsi GC')
ylim([0 .4])

saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/GCscorevproptuned_alldata_contra-ipsi.eps')
