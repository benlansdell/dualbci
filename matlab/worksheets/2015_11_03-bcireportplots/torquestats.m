conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');

manualcontrol = fetch(exec(conn, ['SELECT `File name MCP 1` FROM `rotated_linear_analysis_2014` WHERE `Type of trial` = "control"']))
manualcontrol = manualcontrol.Data;
braincontrol_rotated = fetch(exec(conn, ['SELECT `nev file` FROM `rotated_linear_analysis_2013` WHERE `im tuning type` = 1 OR `im tuning type` = 4']))
braincontrol_rotated = braincontrol_rotated.Data;
braincontrol_control = fetch(exec(conn, ['SELECT `File name BC` FROM `rotated_linear_analysis_2014` WHERE `Type of trial` = "control"']));
braincontrol_control = braincontrol_control.Data;

dualcontrol = fetch(exec(conn, ['SELECT `dualrecording` FROM `experiment_tuning` LIMIT 20']));
dualcontrol = dualcontrol.Data;

braincontrol = [braincontrol_rotated(1:20); braincontrol_control(1:20)];
condition = [ones(size(braincontrol_rotated(1:20))); 2*ones(size(braincontrol_control(1:20)))];
conds = {'rotated', 'control'};

binsize = 0.002;
samplerate = 1/binsize;
angles = {};
torque = [];
dtorque = [];
ddtorque = [];
%Load torque data
for idx = 1:length(braincontrol)
	cnd = conds{condition(idx)};
	nevfile = ['./blackrock/' braincontrol{idx}]
	ns3file = [nevfile(1:end-3) 'ns3'];
	chans = findTorqueChannels(nevfile);
	NEV = openNEV(nevfile, 'nosave');
	dur = NEV.MetaTags.DataDurationSec;
	nB = dur/binsize;
	NS3 = openNSx(ns3file, 'read', ['c:' num2str(chans(1)) ':' num2str(chans(2))]);
	nsxtorque = double(NS3.Data);
	nsxsamplerate = double(NS3.MetaTags.SamplingFreq);
	%Switch sign of FE axis for coordinate consistency
	nsxtorque(1,:)=-nsxtorque(1,:);
	nsxpts = ((1:size(nsxtorque,2))-1)/nsxsamplerate;
	pts = ((1:nB)-1)/samplerate;
	%Smoothing parameter
	%p = 1/(1+binsize^3/0.001);
	p = 0.9999;
	clear torqueii; clear dtorqueii; clear ddtorqueii;
	for j=1:2
		%Scale from uint16 value to proportion
		nsxtorque(j,:) = nsxtorque(j,:)/(2^15);
		%Subtract mean
		nsxtorque(j,:) = nsxtorque(j,:)-mean(nsxtorque(j,:));
		%Smooth w spline
		sp = csaps(nsxpts, nsxtorque(j,:), p);
		torqueii(:,j) = fnval(sp, pts);
		%Compute velocity of spline
		vel = fnder(sp);
		dtorqueii(:,j) = fnval(vel, pts);
		%Compute accel of spline
		accel = fnder(vel);
		ddtorqueii(:,j) = fnval(accel, pts);
	end
	torque = [torque; torqueii];
	dtorque = [dtorque; dtorqueii];
	ddtorque = [ddtorque; ddtorqueii];
end


torqueMC = [];
dtorqueMC = [];
ddtorqueMC = [];
%Load torque data
for idx = 1:length(manualcontrol)
	if strcmp(manualcontrol{idx}, 'null')
		continue
	end
	nevfile = ['./blackrock/' manualcontrol{idx}]
	ns3file = [nevfile(1:end-3) 'ns3'];
	chans = findTorqueChannels(nevfile);
	NEV = openNEV(nevfile, 'nosave');
	dur = NEV.MetaTags.DataDurationSec;
	nB = dur/binsize;
	NS3 = openNSx(ns3file, 'read', ['c:' num2str(chans(1)) ':' num2str(chans(2))]);
	nsxtorque = double(NS3.Data);
	nsxsamplerate = double(NS3.MetaTags.SamplingFreq);
	%Switch sign of FE axis for coordinate consistency
	nsxtorque(1,:)=-nsxtorque(1,:);
	nsxpts = ((1:size(nsxtorque,2))-1)/nsxsamplerate;
	pts = ((1:nB)-1)/samplerate;
	%Smoothing parameter
	%p = 1/(1+binsize^3/0.001);
	p = 0.9999;
	clear torqueii; clear dtorqueii; clear ddtorqueii;
	for j=1:2
		%Scale from uint16 value to proportion
		nsxtorque(j,:) = nsxtorque(j,:)/(2^15);
		%Subtract mean
		nsxtorque(j,:) = nsxtorque(j,:)-mean(nsxtorque(j,:));
		%Smooth w spline
		sp = csaps(nsxpts, nsxtorque(j,:), p);
		torqueii(:,j) = fnval(sp, pts);
		%Compute velocity of spline
		vel = fnder(sp);
		dtorqueii(:,j) = fnval(vel, pts);
		%Compute accel of spline
		accel = fnder(vel);
		ddtorqueii(:,j) = fnval(accel, pts);
	end
	torqueMC = [torqueMC; torqueii];
	dtorqueMC = [dtorqueMC; dtorqueii];
	ddtorqueMC = [ddtorqueMC; ddtorqueii];
end


torqueDC = [];
dtorqueDC = [];
ddtorqueDC = [];
%Load torque data
for idx = 1:length(dualcontrol)
	if strcmp(dualcontrol{idx}, 'null')
		continue
	end
	nevfile = ['./blackrock/' dualcontrol{idx}]
	ns3file = [nevfile(1:end-3) 'ns3'];
	chans = findTorqueChannels(nevfile);
	NEV = openNEV(nevfile, 'nosave');
	dur = NEV.MetaTags.DataDurationSec;
	nB = dur/binsize;
	NS3 = openNSx(ns3file, 'read', ['c:' num2str(chans(1)) ':' num2str(chans(2))]);
	nsxtorque = double(NS3.Data);
	nsxsamplerate = double(NS3.MetaTags.SamplingFreq);
	%Switch sign of FE axis for coordinate consistency
	nsxtorque(1,:)=-nsxtorque(1,:);
	nsxpts = ((1:size(nsxtorque,2))-1)/nsxsamplerate;
	pts = ((1:nB)-1)/samplerate;
	%Smoothing parameter
	%p = 1/(1+binsize^3/0.001);
	p = 0.9999;
	clear torqueii; clear dtorqueii; clear ddtorqueii;
	for j=1:2
		%Scale from uint16 value to proportion
		nsxtorque(j,:) = nsxtorque(j,:)/(2^15);
		%Subtract mean
		nsxtorque(j,:) = nsxtorque(j,:)-mean(nsxtorque(j,:));
		%Smooth w spline
		sp = csaps(nsxpts, nsxtorque(j,:), p);
		torqueii(:,j) = fnval(sp, pts);
		%Compute velocity of spline
		vel = fnder(sp);
		dtorqueii(:,j) = fnval(vel, pts);
		%Compute accel of spline
		accel = fnder(vel);
		ddtorqueii(:,j) = fnval(accel, pts);
	end
	torqueDC = [torqueDC; torqueii];
	dtorqueDC = [dtorqueDC; dtorqueii];
	ddtorqueDC = [ddtorqueDC; ddtorqueii];
end


save('./worksheets/2015_11_03-bcireportplots/torquestats.mat');
load('./worksheets/2015_11_03-bcireportplots/torquestats.mat');

%Plot the PCA components over the distribution
[coeffMC, scoreMC, latentMC] = pca(torqueMC);
[coeff, score, latent] = pca(torque);
[coeffDC, scoreDC, latentDC] = pca(torqueDC);

meanBC = mean(torque);
meanMC = mean(torqueMC);
meanDC = mean(torqueDC);

h = figure
smoothhist2D(torqueMC, 5, [100 100], .05)
hold on
plot(meanMC(1), meanMC(2), 'b.')
plot([meanMC(1), meanMC(1)+sqrt(latentMC(1))*coeffMC(1,1)], [meanMC(2), meanMC(2)+sqrt(latentMC(1))*coeffMC(2,1)], 'r')
plot([meanMC(1), meanMC(1)+sqrt(latentMC(2))*coeffMC(1,2)], [meanMC(2), meanMC(2)+sqrt(latentMC(2))*coeffMC(2,2)], 'r')
xlim([-.6 .6])
ylim([-.6 .6])
xlabel('x')
ylabel('y')
title('Manual control')
saveplot(gcf, './worksheets/2015_11_03-bcireportplots/MCtorque.eps')

g = figure
smoothhist2D(torque, 5, [100 100], .05)
hold on 
plot(meanBC(1), meanBC(2), 'b.')
plot([meanBC(1), meanBC(1)+sqrt(latent(1))*coeff(1,1)], [meanBC(2), meanBC(2)+sqrt(latent(1))*coeff(2,1)], 'r')
plot([meanBC(1), meanBC(1)+sqrt(latent(2))*coeff(1,2)], [meanBC(2), meanBC(2)+sqrt(latent(2))*coeff(2,2)], 'r')
xlim([-1 1])
ylim([-1 1])
xlabel('x')
ylabel('y')
title('Brain control')
saveplot(gcf, './worksheets/2015_11_03-bcireportplots/BCtorque.eps')

i = figure
smoothhist2D(torqueDC, 5, [100 100], .05)
hold on 
plot(meanDC(1), meanDC(2), 'b.')
plot([meanDC(1), meanDC(1)+sqrt(latentDC(1))*coeffDC(1,1)], [meanDC(2), meanDC(2)+sqrt(latentDC(1))*coeffDC(2,1)], 'r')
plot([meanDC(1), meanDC(1)+sqrt(latentDC(2))*coeffDC(1,2)], [meanDC(2), meanDC(2)+sqrt(latentDC(2))*coeffDC(2,2)], 'r')
xlim([-1 1])
ylim([-1 1])
xlabel('x')
ylabel('y')
title('Dual control')
saveplot(gcf, './worksheets/2015_11_03-bcireportplots/DCtorque.eps')

%Compute speed
speedMC = sqrt(dtorqueMC(:,1).^2 + dtorqueMC(:,2).^2);
speedBC = sqrt(dtorque(:,1).^2 + dtorque(:,2).^2);
speedDC = sqrt(dtorqueDC(:,1).^2 + dtorqueDC(:,2).^2);
figure;
subplot(1,3,1)
hist(speedBC, 500); xlim([0 5])
xlabel('Torque speed (BC)')
title(['Median:' num2str(median(speedBC))])
subplot(1,3,2)
hist(speedMC, 100); xlim([0 5])
xlabel('Torque speed (MC)')
title(['Median:' num2str(median(speedMC))])
subplot(1,3,3)
hist(speedDC, 100); xlim([0 5])
xlabel('Torque speed (DC)')
title(['Median:' num2str(median(speedDC))])
saveplot(gcf, './worksheets/2015_11_03-bcireportplots/torqueSpeed.eps', 'eps', [8 4])