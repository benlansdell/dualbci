%Optimal lag on position and velocity models
%Preprocess data (timebis of 2ms, smooth trajectories)
%Apply no offset
nevfile = './testdata/20130117SpankyUtah001.nev';
binsize = 0.002;
dt_sp = binsize;
dt_pos = 0.05;
dt_vel = 0.05;
nK_pos = 1;
nK_vel = 1;
nK_sp = 100;
threshold = 5;
const = 'on';
offsets = linspace(-0.5, 0.5, 10);

nK_pos = 0;
processed = preprocess_spline(nevfile, binsize, threshold, offsets(1));
data = filters_sp_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
NMS = size(data.y,2);
model = MLE_glmfit(data, const);
dev_MS = deviance(model, data);

%Fit different filters for different timebin sizes
nK_pos = 1;
nU = length(processed.unitnames);
devs = zeros(length(offsets), nU);
for j = 1:length(offsets);
	offset = offsets(j);
	processed = preprocess_spline(nevfile, binsize, threshold, offset);
	data = filters_sp_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
	N = size(data.y, 2);
	model = MLE_glmfit(data, const);
	devs(j,:) = deviance(model, data)/N*NMS;
end

dev_MS1 = repmat(dev_MS, length(offsets),1);

%Plot deviance as a function of lag for each unit 
plot(offsets, dev_MS1-devs);
xlabel('lag(s)')
ylabel('\Delta D')
title('Position tuning')
saveplot(gcf, './worksheets/11_15_2014/plots/lags_pos.eps')

%Find max lag for each

%Redo for velocity based tuning
%Fit different filters for different timebin sizes
devs_vel = zeros(length(offsets), nU);
nK_vel = 1;
for j = 1:length(offsets);
	offset = offsets(j);
	processed = preprocess_spline(nevfile, binsize, threshold, offset);
	data = filters_sp_vel(processed, nK_sp, nK_vel, dt_sp, dt_vel);
	N = size(data.y, 2);
	model = MLE_glmfit(data, const);
	devs_vel(j,:) = deviance(model, data)/N*NMS;
end

%Plot deviance as a function of lag for each unit 
plot(offsets, dev_MS1 - devs_vel);
xlabel('lag(s)')
ylabel('\Delta D')
title('Velocity tuning')
saveplot(gcf, './worksheets/11_15_2014/plots/lags_vel.eps')
save('./worksheets/11_15_2014/data.mat', 'devs', 'dev_MS', 'dev_MS1', 'devs_vel', 'offsets');
