%%%%%%%%%%%%%%%%
%Manual control%
%%%%%%%%%%%%%%%%
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

data = fetch(exec(conn, ['SELECT fl.`dir`, fl2.`dir`, eg.`score`, f.unit, eg.fromunit FROM fits f '...
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
' WHERE f.modelID = 1 AND f2.modelID = 2 AND f3.modelID = 1 AND fl.size > .2 AND fl2.size > .2']));

d = cell2mat(data.Data(:,1:3));
units = floor(cellfun(@str2num,data.Data(:,4:5)));
include = units(:,1)~=units(:,2);
d = d(include,:);

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

nbetween = diff(nabove);
nsimtunedbetween = diff(nsimtuned);

[gcthresholds; nabove; nsimtuned]
ii = 1:find(nabove<30, 1);
plot(gcthresholds(ii), proptuned(ii), gcthresholds, ones(size(gcthresholds))*(10/180), ':k')
xlabel('GC score')
ylabel('Proportion tuned within +/- 10 degrees')
title('Manual control. Ipsi+contra')
ylim([0 .4])
saveplot(gcf, './worksheets/2015_08_08-GrangerVsTuningAngle/gcscorevproptuned_MC.eps')

%Note the non-uniform distribution of preferred angles (even with the linear model)
hist(d(:,2), 50) 