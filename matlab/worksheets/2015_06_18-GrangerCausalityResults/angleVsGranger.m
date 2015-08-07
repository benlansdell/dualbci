%Fetch paramcode to load
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
data = exec(conn, ['SELECT ge.fromunit, ge.score, rec.successrate, f.`nev file`, fl.dir, bci.angle, MOD(fl.dir-bci.angle, 2*PI()), fl.size, u.firingrate '...
	'FROM estimates_granger ge '...
	'INNER JOIN fits f ON f.id = ge.id '...
	'INNER JOIN recordings rec ON rec.`nev file` = f.`nev file` '...
	'INNER JOIN bci_units bci ON bci.ID = f.`nev file` AND bci.unit = ge.fromunit '...
	'INNER JOIN experiment_tuning al ON al.`1DBCrecording` = f.`nev file` '...
	'INNER JOIN fits f1 ON f1.`nev file` = al.`manualrecording` AND f1.`modelID` = 1 and f1.unit = ge.fromunit '...
	'INNER JOIN fits_linear fl ON f1.`id` = fl.id '...
	'INNER JOIN units u ON u.`nev file` = f1.`nev file` AND u.`unit` = f1.unit WHERE f.modelID = 3 '...
	'INNER JOIN estimates_parameters fl1 ON fl1.`id` = ']);

data = fetch(data);
data = data.Data;

hold on 
delangle = [];
gc = [];
tuningsize = [];
firing = [];

for i = 1:(size(data,1)-1)
	if strcmp(data{i, 4}, data{i+1, 4})
		idx = length(gc)+1;
		%delangle(idx+1) = data{idx+1,7};
		total = data{idx,2}+data{idx+1,2};
		[gc(idx),maxidx] = max([data{idx,2},data{idx+1,2}]/total);
		%gc(idx+1) = data{idx+1,2}/total;
		delangle(idx) = data{idx-1+maxidx,7};

		tuningsize(idx) = data{idx,8}-data{idx+1,8};
		%tuningsize(idx+1) = data{idx+1,8};
		firing(idx) = data{idx,9}-data{idx+1,9};
		%firing(idx+1) = data{idx+1,9};
	end
end
%Make between pi and -pi
delangle = mod(delangle, 2*pi);
delangle(delangle>pi) = delangle(delangle>pi)-2*pi;
figure
plot(delangle, gc, '.')
f=smooth(delangle, gc, 0.5, 'loess');
[xx,ind] = sort(delangle);
hold on
plot(xx, f(ind)', '-');


xlabel('Change in angle')
ylabel('Proportion of granger score')
figure
plot(tuningsize, gc, '.')
xlabel('Difference in tuning')
ylabel('Proportion of granger score')
figure
plot(firing, gc, '.')
xlabel('Difference in firing')
ylabel('Proportion of granger score')



scatter(180*x/pi,180*y/pi,[],c, 'filled');
xlabel('Delta angle 1')
ylabel('Delta angle 2')
xlim([-180 180])
ylim([-180 180])
colorbar
%saveplot(gcf, './worksheets/2015_06_18/angleVsGranger_log.eps')
saveplot(gcf, './worksheets/2015_06_18/angleVsGranger.eps')