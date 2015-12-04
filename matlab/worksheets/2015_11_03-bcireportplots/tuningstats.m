conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gather data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MC
%%%
data = fetch(exec(conn, ['SELECT fl1.r2, fl1.dir, fl1.size, fl1.r2out, flin1.dev FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`manualrecording`'...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'WHERE flin1.modelID = 1']));
MCtorq = cell2mat(data.Data(:,1:5));

data = fetch(exec(conn, ['SELECT fl1.r2, fl1.dir, fl1.size, fl1.r2out, flin1.dev FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`manualrecording`'...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'WHERE flin1.modelID = 30']));
MCtorqV = cell2mat(data.Data(:,1:5));

data = fetch(exec(conn, ['SELECT fl1.r2, fl1.dir, fl1.size, fl1.r2out, flin1.dev FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`manualrecording`'...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'WHERE flin1.modelID = 7']));
MCcursor = cell2mat(data.Data(:,1:5));

data = fetch(exec(conn, ['SELECT fl1.r2, fl1.dir, fl1.size, fl1.r2out, flin1.dev FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`manualrecording`'...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'WHERE flin1.modelID = 28']));
MCcursorV = cell2mat(data.Data(:,1:5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BC
%%%
data = fetch(exec(conn, ['SELECT fl1.r2, fl1.dir, fl1.size, fl1.r2out, flin1.dev FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`1DBCrecording`'...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'WHERE flin1.modelID = 1']));
BCtorq = cell2mat(data.Data(:,1:5));

data = fetch(exec(conn, ['SELECT fl1.r2, fl1.dir, fl1.size, fl1.r2out, flin1.dev FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`1DBCrecording`'...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'WHERE flin1.modelID = 30']));
BCtorqV = cell2mat(data.Data(:,1:5));

data = fetch(exec(conn, ['SELECT fl1.r2, fl1.dir, fl1.size, fl1.r2out, flin1.dev FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`1DBCrecording`'...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'WHERE flin1.modelID = 7']));
BCcursor = cell2mat(data.Data(:,1:5));

data = fetch(exec(conn, ['SELECT fl1.r2, fl1.dir, fl1.size, fl1.r2out, flin1.dev FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`1DBCrecording`'...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'WHERE flin1.modelID = 28']));
BCcursorV = cell2mat(data.Data(:,1:5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DC
%%%
data = fetch(exec(conn, ['SELECT fl1.r2, fl1.dir, fl1.size, fl1.r2out, flin1.dev FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`dualrecording`'...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'WHERE flin1.modelID = 1']));
DCtorq = cell2mat(data.Data(:,1:5));

data = fetch(exec(conn, ['SELECT fl1.r2, fl1.dir, fl1.size, fl1.r2out, flin1.dev FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`dualrecording`'...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'WHERE flin1.modelID = 30']));
DCtorqV = cell2mat(data.Data(:,1:5));

data = fetch(exec(conn, ['SELECT fl1.r2, fl1.dir, fl1.size, fl1.r2out, flin1.dev FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`dualrecording`'...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'WHERE flin1.modelID = 7']));
DCcursor = cell2mat(data.Data(:,1:5));

data = fetch(exec(conn, ['SELECT fl1.r2, fl1.dir, fl1.size, fl1.r2out, flin1.dev FROM '...
'`experiment_tuning` et1 '...
'INNER JOIN `fits` flin1 '...
'ON flin1.`nev file` = et1.`dualrecording`'...
'INNER JOIN `fits_linear` fl1 '...
'ON flin1.id = fl1.id '...
'WHERE flin1.modelID = 28']));
DCcursorV = cell2mat(data.Data(:,1:5));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tuning histograms%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Tuning angle

%Note the non-uniform distribution of preferred angles (even with the linear model)
clf
subplot(3,4,1)
hist(MCtorq(:,2)/pi*180, 50) 
title('MC torq')
xlabel('angle')
ylabel('frequency')
subplot(3,4,2)
hist(MCtorqV(:,2)/pi*180, 50)
title('MC torq vel')
xlabel('angle')
ylabel('frequency')
subplot(3,4,3)
hist(MCcursor(:,2)/pi*180, 50)
title('MC cursor')
xlabel('angle')
ylabel('frequency')
subplot(3,4,4)
hist(MCcursorV(:,2)/pi*180, 50) 
title('MC cursor vel')
xlabel('angle')
ylabel('frequency')
subplot(3,4,5)
hist(BCtorq(:,2)/pi*180, 50)
title('BC torq')
xlabel('angle')
ylabel('frequency')
subplot(3,4,6)
hist(BCtorqV(:,2)/pi*180, 50)
title('BC torq vel')
xlabel('angle')
ylabel('frequency')
subplot(3,4,7)
hist(BCcursor(:,2)/pi*180, 50) 
title('BC cursor')
xlabel('angle')
ylabel('frequency')
subplot(3,4,8)
hist(BCcursorV(:,2)/pi*180, 50)
title('BC cursor vel')
xlabel('angle')
ylabel('frequency')
subplot(3,4,9)
hist(DCtorq(:,2)/pi*180, 50)
title('DC torq')
xlabel('angle')
ylabel('frequency')
subplot(3,4,10)
hist(DCtorqV(:,2)/pi*180, 50) 
title('DC torq vel')
xlabel('angle')
ylabel('frequency')
subplot(3,4,11)
hist(DCcursor(:,2)/pi*180, 50)
title('DC cursor')
xlabel('angle')
ylabel('frequency')
subplot(3,4,12)
hist(DCcursorV(:,2)/pi*180, 50)
title('DC cursor vel')
xlabel('angle')
ylabel('frequency')
saveplot(gcf, './worksheets/2015_11_03-bcireportplots/hist_tuning_direction.eps', 'eps', [10 10])


%Tuning size

%Note the non-uniform distribution of preferred angles (even with the linear model)
clf
subplot(3,4,1)
hist(MCtorq(:,1), 50) 
title('MC torq')
xlabel('R^2'); xlim([0 0.15]);
ylabel('frequency')
subplot(3,4,2)
hist(MCtorqV(:,1), 50)
title('MC torq vel')
xlabel('R^2'); xlim([0 0.15]);
ylabel('frequency')
subplot(3,4,3)
hist(MCcursor(:,1), 50)
title('MC cursor')
xlabel('R^2'); xlim([0 0.15]);
ylabel('frequency')
subplot(3,4,4)
hist(MCcursorV(:,1), 50) 
title('MC cursor vel')
xlabel('R^2'); xlim([0 0.15]);
ylabel('frequency')
subplot(3,4,5)
hist(BCtorq(:,1), 50)
title('BC torq')
xlabel('R^2'); xlim([0 0.15]);
ylabel('frequency')
subplot(3,4,6)
hist(BCtorqV(:,1), 50)
title('BC torq vel')
xlabel('R^2'); xlim([0 0.15]);
ylabel('frequency')
subplot(3,4,7)
hist(BCcursor(:,1), 50) 
title('BC cursor')
xlabel('R^2'); xlim([0 0.15]);
ylabel('frequency')
subplot(3,4,8)
hist(BCcursorV(:,1), 50)
title('BC cursor vel')
xlabel('R^2'); xlim([0 0.15]);
ylabel('frequency')
subplot(3,4,9)
hist(DCtorq(:,1), 50)
title('DC torq')
xlabel('R^2'); xlim([0 0.15]);
ylabel('frequency')
subplot(3,4,10)
hist(DCtorqV(:,1), 50) 
title('DC torq vel')
xlabel('R^2'); xlim([0 0.15]);
ylabel('frequency')
subplot(3,4,11)
hist(DCcursor(:,1), 50)
title('DC cursor')
xlabel('R^2'); xlim([0 0.15]);
ylabel('frequency')
subplot(3,4,12)
hist(DCcursorV(:,1), 50)
title('DC cursor vel')
xlabel('R^2'); xlim([0 0.15]);
ylabel('frequency')
saveplot(gcf, './worksheets/2015_11_03-bcireportplots/hist_tuning_R^2.eps', 'eps', [10 10])


%Box plot of data
MCalldata = [MCtorq(:,1); MCtorqV(:,1); MCcursor(:,1); MCcursorV(:,1)];
MClengths = [size(MCtorq(:,1)); size(MCtorqV(:,1)); size(MCcursor(:,1)); size(MCcursorV(:,1))];
MCnames = {'MC torq', 'MC torq vel', 'MC curs', 'MC curs vel'};
MCconds = {};
for idx = 1:length(MCnames)
	nE = MClengths(idx);
	name = MCnames{idx};
	MCconds = [MCconds; repmat({name}, nE,1)];
end

BCalldata = [BCtorq(:,1); BCtorqV(:,1); BCcursor(:,1); BCcursorV(:,1)];
BClengths = [size(BCtorq(:,1)); size(BCtorqV(:,1)); size(BCcursor(:,1)); size(BCcursorV(:,1))];
BCnames = {'BC torq', 'BC torq vel', 'BC curs', 'BC curs vel'};
BCconds = {};
for idx = 1:length(BCnames)
	nE = BClengths(idx);
	name = BCnames{idx};
	BCconds = [BCconds; repmat({name}, nE,1)];
end

DCalldata = [DCtorq(:,1); DCtorqV(:,1); DCcursor(:,1); DCcursorV(:,1)];
DClengths = [size(DCtorq(:,1)); size(DCtorqV(:,1)); size(DCcursor(:,1)); size(DCcursorV(:,1))];
DCnames = {'DC torq', 'DC torq vel', 'DC curs', 'DC curs vel'};
DCconds = {};
for idx = 1:length(DCnames)
	nE = DClengths(idx);
	name = DCnames{idx};
	DCconds = [DCconds; repmat({name}, nE,1)];
end

clf
subplot(3,1,1)
boxplot(MCalldata, MCconds);
ylim([0 .06])
ylabel('R^2')
subplot(3,1,2)
boxplot(BCalldata, BCconds);
ylim([0 .06])
ylabel('R^2')
subplot(3,1,3)
boxplot(DCalldata, DCconds);
ylim([0 .06])
ylabel('R^2')
saveplot(gcf, './worksheets/2015_11_03-bcireportplots/hist_tuning_R^2_boxplot.eps', 'eps', [10 10])
