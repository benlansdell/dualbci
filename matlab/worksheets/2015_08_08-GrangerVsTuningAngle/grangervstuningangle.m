conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

data = fetch(exec(conn, ['SELECT ge.score, fl1.dir, fl2.dir, MOD(fl1.dir- fl2.dir, 2*PI()), flin1.`nev file`, ge.`fromunit`, flin1.`unit` FROM `fits` fgranger '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = fgranger.`nev file` AND flin1.`unit` = fgranger.`unit` '...
'INNER JOIN `estimates_granger` ge '...
'ON fgranger.id = ge.id '...
'INNER JOIN `fits` flin2 '...
'ON flin2.`nev file` = fgranger.`nev file` AND flin2.`unit` = ge.`fromunit` '...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'INNER JOIN `fits_linear` fl2 '...
'ON flin2.id = fl2.id '...
'WHERE flin1.modelID = 1 AND flin2.modelID = 1 AND fgranger.modelID = 2']));

d = cell2mat(data.Data(:,1:4));

%Bin them into bins...
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
	score = d(idx,1);
	bin = ceil(nbins*(angle+pi)/(2*pi));
	bins{bin} = [bins{bin},score];
end

nK_sp = 6;
sigthresh = 22.35; %Because: 1-chi2cdf(22.35, nK_sp) = 0.001;
%sigthresh = 12.59; %Because: 1-chi2cdf(12.59, nK_sp) = 0.05;

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
title('Brain+manual. Ipsi+contra')
xlim([-180 180])
ylim([0 40])
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/GCvtuning_alldata.eps')

figure
plot(angles, propabove)
ylabel('Proportion of significant connections. \alpha = 0.001')
ylim([0 .4])
xlabel('Difference in tuning angles between units')
title('Brain+manual. Ipsi+contra')
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/tuningvpropsig_alldata.eps')

figure
plot(angles, maxscores)

figure
boxplot(catted, groups)


%%%%%%%%%%%%%%%%
%Manual control%
%%%%%%%%%%%%%%%%


data = fetch(exec(conn, ['SELECT ge.score, fl1.dir, fl2.dir, MOD(fl1.dir- fl2.dir, 2*PI()), flin1.`nev file`, ge.`fromunit`, flin1.`unit` FROM `fits` fgranger '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = fgranger.`nev file` AND flin1.`unit` = fgranger.`unit` '...
'INNER JOIN `estimates_granger` ge '...
'ON fgranger.id = ge.id '...
'INNER JOIN `fits` flin2 '...
'ON flin2.`nev file` = fgranger.`nev file` AND flin2.`unit` = ge.`fromunit` '...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'INNER JOIN `fits_linear` fl2 '...
'ON flin2.id = fl2.id '...
' INNER JOIN experiment_tuning et '...
' ON et.`manualrecording` = flin1.`nev file` '...
'WHERE flin1.modelID = 1 AND flin2.modelID = 1 AND fgranger.modelID = 2']));

d = cell2mat(data.Data(:,1:4));

%Bin them into bins...
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
	score = d(idx,1);
	bin = ceil(nbins*(angle+pi)/(2*pi));
	bins{bin} = [bins{bin},score];
end

nK_sp = 6;
sigthresh = 22.35; %Because: 1-chi2cdf(22.35, nK_sp) = 0.001;
%sigthresh = 12.59; %Because: 1-chi2cdf(12.59, nK_sp) = 0.05;

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
title('manual. Ipsi+contra')
xlim([-180 180])
ylim([0 40])
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/GCvtuning_MCdata.eps')

figure
plot(angles, propabove)
ylabel('Proportion of significant connections. \alpha = 0.001')
ylim([0 .4])
xlabel('Difference in tuning angles between units')
title('manual. Ipsi+contra')
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/tuningvpropsig_MCdata.eps')


%%%%%%%%%%%%%%%%
%Brain control%
%%%%%%%%%%%%%%%%


data = fetch(exec(conn, ['SELECT ge.score, fl1.dir, fl2.dir, MOD(fl1.dir- fl2.dir, 2*PI()), flin1.`nev file`, ge.`fromunit`, flin1.`unit` FROM `fits` fgranger '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = fgranger.`nev file` AND flin1.`unit` = fgranger.`unit` '...
'INNER JOIN `estimates_granger` ge '...
'ON fgranger.id = ge.id '...
'INNER JOIN `fits` flin2 '...
'ON flin2.`nev file` = fgranger.`nev file` AND flin2.`unit` = ge.`fromunit` '...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'INNER JOIN `fits_linear` fl2 '...
'ON flin2.id = fl2.id '...
' INNER JOIN experiment_tuning et '...
' ON et.`1DBCrecording` = flin1.`nev file` '...
'WHERE flin1.modelID = 1 AND flin2.modelID = 1 AND fgranger.modelID = 2']));

d = cell2mat(data.Data(:,1:4));

%Bin them into bins...
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
	score = d(idx,1);
	bin = ceil(nbins*(angle+pi)/(2*pi));
	bins{bin} = [bins{bin},score];
end

nK_sp = 6;
sigthresh = 22.35; %Because: 1-chi2cdf(22.35, nK_sp) = 0.001;
%sigthresh = 12.59; %Because: 1-chi2cdf(12.59, nK_sp) = 0.05;

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
title('Brain. Ipsi+contra')
xlim([-180 180])
ylim([0 40])
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/GCvtuning_BCdata.eps')

figure
plot(angles, propabove)
ylabel('Proportion of significant connections. \alpha = 0.001')
ylim([0 .4])
xlabel('Difference in tuning angles between units')
title('Brain. Ipsi+contra')
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/tuningvpropsig_BCdata.eps')
