%Compute sum of granger scores for out degree for given unit and nev file
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

data = exec(conn, ['SELECT f.`nev file`, f.`unit`, sum(ge.score), f3.`mse out`, fglmtrqsp.`mse out`, fglmtrqsp.`dev`, '...
	'fglmsp.`mse out`,fglmsp.`dev`,fglmtrqtarsp.`mse out`, fglmtrqtarsp.`dev` '...
	'FROM '...
	'fits f '...
	'INNER JOIN fits f2 '...
    'ON f2.`nev file` = f.`nev file` AND f2.modelID = 2 '...
	'INNER JOIN estimates_granger ge '...
	'ON f2.id = ge.id AND ge.`fromunit` = f.unit '...
	'INNER JOIN experiment_tuning et '...
	'ON f.`nev file` = et.`manualrecording` '...
    'INNER JOIN fits f3 '...
    'ON f3.`nev file` = f.`nev file` AND f3.`unit` = f.`unit` AND f3.`modelID` = 1 '...
    'INNER JOIN fits fglmtrqsp ON fglmtrqsp.`nev file` = et.`manualrecording` AND fglmtrqsp.`modelID` = 14 and fglmtrqsp.unit = ge.fromunit '...
	'INNER JOIN fits fglmsp ON fglmsp.`nev file` = et.`manualrecording` AND fglmsp.`modelID` = 15 and fglmsp.unit = ge.fromunit '...
	'INNER JOIN fits fglmtrqtarsp ON fglmtrqtarsp.`nev file` = et.`manualrecording` AND fglmtrqtarsp.`modelID` = 16 and fglmtrqtarsp.unit = ge.fromunit '...
	'WHERE f.modelID = 2 '...
    'GROUP BY f.`nev file`, f.`unit`'])

data = fetch(data);
data = data.Data;

a = cell2mat(data(:,3:4));
smoothhist2D(a, 2, [200 200]); xlim([40 500]); ylim([0.01 .4]); axis xy; xlabel('Outward sum Granger scores'); ylabel('MSE out (linear model [1])' )
f = LinearModel.fit(a(:,1), a(:,2));
hold on 
plot(f)
xlabel('Outward sum Granger scores'); ylabel('MSE out (linear [1]])' )
saveplot(gcf, './worksheets/2015_08_07-GrangerAndMSEOUT/GrangerOutvsMSE.eps')
%Returns:
%mdl = LinearModel.fit(a(:,1), a(:,2))
%
%mdl = 
%
%
%Linear regression model:
%    y ~ 1 + x1
%
%Estimated Coefficients:
%                   Estimate      SE            tStat     pValue   
%    (Intercept)         0.011     0.0040266    2.7319    0.0063903
%    x1             0.00043869    2.5196e-05    17.411     1.04e-60
%
%
%Number of observations: 1202, Error degrees of freedom: 1200
%Root Mean Squared Error: 0.0647
%R-squared: 0.202,  Adjusted R-Squared 0.201
%F-statistic vs. constant model: 303, p-value = 1.04e-60

%...sig value, but tiny effect....

%Of course, when we look at the mse out for other models... we see something similar:
%And even, disappointingly, we see the same pattern in a model with just the spike history component:

d = cell2mat(data(:,3:end));
clf
%smoothhist2D([d(:,5), d(:,1)], 2, [200 200]); xlim([40 500]); ylim([0.01 .4]); axis xy; xlabel('Outward sum Granger scores'); ylabel('MSE out (linear model [1])' )
f = LinearModel.fit(d(:,1), d(:,5));
hold on 
plot(f)
saveplot(gcf, './worksheets/2015_08_07-GrangerAndMSEOUT/GrangerOutvsMSESp.eps')
xlabel('Outward sum Granger scores'); ylabel('MSE out (linear model [1])' )
clf
%smoothhist2D([d(:,1), d(:,3)], 2, [200 200]); 
%xlim([40 500]); ylim([0.01 .4]); axis xy; 
f = LinearModel.fit(d(:,1), d(:,3));
hold on 
plot(f)
xlabel('Outward sum Granger scores'); ylabel('MSE out (GLM trq sp [14])' )
saveplot(gcf, './worksheets/2015_08_07-GrangerAndMSEOUT/GrangerOutvsMSETrqSp.eps')

clf
hold off 
plot(d(:,1), d(:,6)-d(:,4), '.')
xlabel('Outward sum Granger scores')
ylabel('Change in dev between GLM(sp) and GLM(sptrq)')
title('Note: in all cases dev_{sp} > dev_{trq}')
saveplot(gcf, './worksheets/2015_08_07-GrangerAndMSEOUT/GrangerOutvsDevSp-DevTrqSp.eps')

plot(d(:,1), d(:,5)-d(:,3), '.')
xlabel('Outward sum Granger scores')
ylabel('Change in MSE between GLM(sp) and GLM(sptrq)')
title('Note: in 810/1202 cases, MSE_{sp} > MSE_{sptrq}')
saveplot(gcf, './worksheets/2015_08_07-GrangerAndMSEOUT/GrangerOutvsMSEoutSp-MSETrqSp.eps')

plot(d(:,1), d(:,5)-d(:,7), '.')
xlabel('Outward sum Granger scores')
ylabel('Change in MSE between GLM(sp) and GLM(sptrqtar)')
title('Note: in 916/1202 cases, MSE_{sp} > MSE_{sptrqtar}')
saveplot(gcf, './worksheets/2015_08_07-GrangerAndMSEOUT/GrangerOutvsMSEoutSp-MSETrqTarSp.eps')

plot(d(:,1), d(:,6)-d(:,8), '.')
xlabel('Outward sum Granger scores')
ylabel('Change in dev between GLM(sp) and GLM(sptrqtar)')
title('Note: in all cases dev_{sp} > dev_{tartrq}')
saveplot(gcf, './worksheets/2015_08_07-GrangerAndMSEOUT/GrangerOutvsDevoutSp-TrqTarSp.eps')

%plot(d(:,6)-d(:,8), d(:,6)-d(:,4), '.')


%Compute sum of granger scores for in degree for given unit and nev file
data = exec(conn, ['SELECT f.`nev file`, f.`unit`, sum(ge.score), f2.`mse out`, fglmtrqsp.`mse out`, fglmtrqsp.`dev`, '...
	'fglmsp.`mse out`,fglmsp.`dev`,fglmtrqtarsp.`mse out`, fglmtrqtarsp.`dev` '...
	'FROM '...
	'fits f '...
	'INNER JOIN estimates_granger ge '...
	'ON f.id = ge.id '...
	'INNER JOIN experiment_tuning et '...
	'ON f.`nev file` = et.`manualrecording` '...
	'INNER JOIN fits f2 '...
    'ON f2.`nev file` = f.`nev file` AND f2.`unit` = f.`unit` AND f2.modelID = 1 '...
    'INNER JOIN fits fglmtrqsp ON fglmtrqsp.`nev file` = et.`manualrecording` AND fglmtrqsp.`modelID` = 14 and fglmtrqsp.unit = f.unit '...
	'INNER JOIN fits fglmsp ON fglmsp.`nev file` = et.`manualrecording` AND fglmsp.`modelID` = 15 and fglmsp.unit = f.unit '...
	'INNER JOIN fits fglmtrqtarsp ON fglmtrqtarsp.`nev file` = et.`manualrecording` AND fglmtrqtarsp.`modelID` = 16 and fglmtrqtarsp.unit = f.unit '...
	'WHERE f.modelID = 2 '...
    'GROUP BY f.`nev file`, f.`unit` ' ...
'ORDER BY `f`.`nev file` ASC, `f`.`unit` ASC'])

data = fetch(data);
data = data.Data;

a = cell2mat(data(:,3:4));
smoothhist2D(a, 2, [200 200]); xlim([40 500]); ylim([0.01 .4]); axis xy;
f = LinearModel.fit(a(:,1), a(:,2));
hold on 
plot(f)
xlabel('Inward sum Granger scores'); ylabel('MSE out (linear model [1])' )
saveplot(gcf, './worksheets/2015_08_07-GrangerAndMSEOUT/GrangerInvsMSE.eps')

d = cell2mat(data(:,3:end));
%smoothhist2D([d(:,5), d(:,1)], 2, [200 200]); xlim([40 500]); ylim([0.01 .4]); axis xy; xlabel('Outward sum Granger scores'); ylabel('MSE out (linear model [1])' )
f = LinearModel.fit(d(:,1), d(:,5));
clf
hold on 
plot(f)
xlabel('Inward sum Granger scores'); ylabel('MSE out (GLM sp [15])' )
saveplot(gcf, './worksheets/2015_08_07-GrangerAndMSEOUT/GrangerInvsMSESp.eps')
clf
%smoothhist2D([d(:,3), d(:,1)], 2, [200 200]); xlim([40 500]); ylim([0.01 .4]); axis xy; xlabel('Outward sum Granger scores'); ylabel('MSE out (linear model [1])' )
f = LinearModel.fit(d(:,1), d(:,3));
hold on 
plot(f)
xlabel('Inward sum Granger scores'); ylabel('MSE out (GLM sp trq [14])' )
saveplot(gcf, './worksheets/2015_08_07-GrangerAndMSEOUT/GrangerInvsMSETrqSp.eps')

clf
hold off 
plot(d(:,1), d(:,6)-d(:,4), '.')
xlabel('Inward sum Granger scores')
ylabel('Change in dev between GLM(sp) and GLM(sptrq)')
title('Note: in all cases dev_{sp} > dev_{trq}')
saveplot(gcf, './worksheets/2015_08_07-GrangerAndMSEOUT/GrangerInvsDevSp-DevTrqSp.eps')

plot(d(:,1), d(:,5)-d(:,3), '.')
xlabel('Inward sum Granger scores')
ylabel('Change in MSE between GLM(sp) and GLM(sptrq)')
title('Note: in 810/1202 cases, MSE_{sp} > MSE_{sptrq}')
saveplot(gcf, './worksheets/2015_08_07-GrangerAndMSEOUT/GrangerInvsMSEoutSp-MSETrqSp.eps')

plot(d(:,1), d(:,5)-d(:,7), '.')
xlabel('Inward sum Granger scores')
ylabel('Change in MSE between GLM(sp) and GLM(sptrqtar)')
title('Note: in 916/1202 cases, MSE_{sp} > MSE_{sptrqtar}')
saveplot(gcf, './worksheets/2015_08_07-GrangerAndMSEOUT/GrangerInvsMSEoutSp-MSETrqTarSp.eps')

plot(d(:,1), d(:,6)-d(:,8), '.')
xlabel('Inward sum Granger scores')
ylabel('Change in dev between GLM(sp) and GLM(sptrqtar)')
title('Note: in all cases dev_{sp} > dev_{tartrq}')
saveplot(gcf, './worksheets/2015_08_07-GrangerAndMSEOUT/GrangerInvsDevoutSp-TrqTarSp.eps')







%Compute sum of granger scores for out degree for given unit and nev file
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

data = exec(conn, ['SELECT f.`nev file`, f.`unit`, sum(ge.score), f3.`mse out`, fglmtrqsp.`mse out`, fglmtrqsp.`dev`, '...
	'fglmsp.`mse out`,fglmsp.`dev` '...
	'FROM '...
	'fits f '...
	'INNER JOIN fits f2 '...
    'ON f2.`nev file` = f.`nev file` AND f2.modelID = 2 '...
	'INNER JOIN estimates_granger ge '...
	'ON f2.id = ge.id AND ge.`fromunit` = f.unit '...
	'INNER JOIN experiment_tuning et '...
	'ON f.`nev file` = et.`manualrecording` '...
    'INNER JOIN fits f3 '...
    'ON f3.`nev file` = f.`nev file` AND f3.`unit` = f.`unit` AND f3.`modelID` = 1 '...
    'INNER JOIN fits fglmtrqsp ON fglmtrqsp.`nev file` = et.`manualrecording` AND fglmtrqsp.`modelID` = 4 and fglmtrqsp.unit = ge.fromunit '...
	'INNER JOIN fits fglmsp ON fglmsp.`nev file` = et.`manualrecording` AND fglmsp.`modelID` = 5 and fglmsp.unit = ge.fromunit '...
	'WHERE f.modelID = 2 '...
    'GROUP BY f.`nev file`, f.`unit`'])

data = fetch(data);
data = data.Data;

d = cell2mat(data(:,3:end));
clf
f = LinearModel.fit(d(:,1), d(:,3));
hold on 
plot(f)
xlabel('Outward sum Granger scores'); ylabel('MSE out (GLM trq sp [4])' )
title('Intercept not significantly non-zero')
saveplot(gcf, './worksheets/2015_08_07-GrangerAndMSEOUT/GrangerOutvsMSETrqSp_long.eps')

clf
f = LinearModel.fit(d(:,1), d(:,5));
hold on 
plot(f)
xlabel('Outward sum Granger scores'); ylabel('MSE out (GLM sp [5])' )
title('Intercept not significantly non-zero')
saveplot(gcf, './worksheets/2015_08_07-GrangerAndMSEOUT/GrangerOutvsMSESp_long.eps')

clf
plot(d(:,5), d(:,3), '.'); 
xlabel('MSE out GLM sp');
ylabel('MSE out GLM trq sp');
title('582/1202 cases have MSE_{sp} < MSE_{trqsp}');
saveplot(gcf, './worksheets/2015_08_07-GrangerAndMSEOUT/MSEsp_MSEtrqsp_long.eps')
