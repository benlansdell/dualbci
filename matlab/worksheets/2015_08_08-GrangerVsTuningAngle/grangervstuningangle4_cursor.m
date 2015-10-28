conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%%%%%%%%%%%%%%%%
%Manual control%
%%%%%%%%%%%%%%%%

data = fetch(exec(conn, ['SELECT fl.`dir` '...
'	FROM `fits` f '...
' INNER JOIN `fits_linear` fl'...
' ON f.id = fl.id '...
' INNER JOIN recordings rec '...
' ON rec.`nev file` = f.`nev file` '...
' WHERE f.modelID = 7 AND rec.tasktype = "brain"']));

d = cell2mat(data.Data(:,1));

gcthresholds = 0:5:200;
anglethreshs = [10 20 30 40 50 60 70 80 90]/180*pi;
proptuned = zeros(length(anglethreshs),length(gcthresholds));
diffangle = mod(d(:,1)-d(:,2), 2*pi);
diffangle(diffangle < 0) = diffangle(diffangle < 0) + 2*pi;
diffangle(diffangle > pi) = diffangle(diffangle > pi) - 2*pi;
diffangle = abs(diffangle);

nabove = [];
nsimtuned = [];

for i = 1:length(anglethreshs)
	for idx = 1:length(gcthresholds)
		%Find those above threshold
		gcthresh = gcthresholds(idx);
		gcabove = d(:,3)>gcthresh;
		nabove(i,idx) = sum(gcabove);
		nsimtuned(i,idx) = sum(diffangle(gcabove)<anglethreshs(i));
		%Find the prop of those whose difference in tuning angle is less than anglethresh
		proptuned(i,idx) = nsimtuned(i,idx)/nabove(i,idx);
	end
end

nbetween = diff(nabove);
nsimtunedbetween = diff(nsimtuned);

figure
[gcthresholds; nabove; nsimtuned]
for i = 1:length(anglethreshs)
	ii = 1:find(nabove(i,:)<30, 1);
	hold on
	plot(gcthresholds(ii), proptuned(i,ii), gcthresholds, ones(size(gcthresholds))*(10/180), ':k')
end
xlabel('GC score')
ylabel('Proportion tuned within +/- 10 degrees')
title('Manual control. Ipsi+contra')
ylim([0 .4])
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/gcscorevproptuned_MC.eps')

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
	hist(d(gcabove,2), 50) 
	title(['GC threshold: ' num2str(gcthresh)])
	pause
	%[f,xi] = ksdensity(d(gcabove,2), 'width', 0.1);
	%hold on
	%plot(xi, f, 'Color', cmap(idx,:));
end

ii = 1:find(nabove<30, 1);
plot(gcthresholds(ii), absangle(ii), gcthresholds(ii), absangle(ii)-stdabsangle(ii), gcthresholds(ii), absangle(ii)+stdabsangle(ii), gcthresholds, ones(size(gcthresholds))*(10/180), ':k')
xlabel('GC score')
ylabel('Mean |diff in angle|')
title('Manual control. Ipsi+contra')
ylim([0 pi])


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
' WHERE f.modelID = 7 AND f2.modelID = 2 AND f3.modelID = 7']));

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

