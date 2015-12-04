conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

all_data = fetch(exec(conn, ['SELECT flin1.`mse out`, flin1.`dev`, fl1.`r2`, fl1.`r2out`, fl1.`size` FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`manualrecording` '...
'INNER JOIN `fits_linear` fl1 '...
'ON fl1.id = flin1.id '...
'WHERE flin1.modelID = 1']));
all_d = cell2mat(all_data.Data(:,1:5));

%Plot of deviance vs MSE out of sample error:
plot(all_d(:,1), all_d(:,2), '.')
xlabel('MSE out');
ylabel('Deviance')
saveplot(gcf, './worksheets/2015_10_22-lineartuningsignificance/mseoutVsDevMC.eps')

%Plot of deviance vs MSE out of sample error:
plot(all_d(:,3), all_d(:,2), '.')
xlabel('R2 within sample')
ylabel('Deviance')
saveplot(gcf, './worksheets/2015_10_22-lineartuningsignificance/R2withinVsDevMC.eps')

plot(all_d(:,1), all_d(:,4), '.')
xlabel('MSE out of sample')
ylabel('R2 out of sample')
saveplot(gcf, './worksheets/2015_10_22-lineartuningsignificance/mseoutVsR2out.eps')

plot(all_d(:,3), all_d(:,4), '.')
xlabel('R2 within sample')
ylabel('R2 out of sample')
saveplot(gcf, './worksheets/2015_10_22-lineartuningsignificance/r2withinVsr2out.eps')

plot(all_d(:,3), all_d(:,5), '.')
xlabel('R2 within sample')
ylabel('"Tuning size"')
saveplot(gcf, './worksheets/2015_10_22-lineartuningsignificance/r2withinVsTuningSize.eps')

plot(all_d(:,2), all_d(:,5), '.')
xlabel('Dev')
ylabel('"Tuning size"')
saveplot(gcf, './worksheets/2015_10_22-lineartuningsignificance/DevVsTuningSize.eps')

plot(all_d(:,1), all_d(:,5), '.')
xlabel('MSE out')
ylabel('"Tuning size"')
saveplot(gcf, './worksheets/2015_10_22-lineartuningsignificance/MSEoutVsTuningSize.eps')

plot(all_d(:,4), all_d(:,5), '.')
xlabel('R2 out of sample')
ylabel('"Tuning size"')
saveplot(gcf, './worksheets/2015_10_22-lineartuningsignificance/r2outVsTuningSize.eps')

%Conclusion: R2 within and R2 without are fairly correlated
%Deviance and MSE out are fairly correlated
%But comparing a normalized to non-normalized quantity (R2within or without to Deviance, for example) are fairly
%uncorrelated. 'Tuning size' is most correlated with R2 values (normalized), not with deviance or MSE. --> interesting!
%Suggests tuning size is an OK indication of tuning to task
%Going forward: Use normalized quantities... R2 and tuning size.