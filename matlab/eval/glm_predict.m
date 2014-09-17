function [corrs, simTrain] = glm_predict(processed, data, model, fn_out)
	%Predict spike trains for GLM given a model and cursor data. Compare the spike trains to the actual spike train generated
	%
	%Usage:
	%	y = glm_predict(processed, model, data, fn_out)
	%     
	%Input:
	%	processed = a structure output from a preprocess function
	%	data = a structure of stimulus and spike history data from ./models
	%	model = a structure of fit coefficients from MLE_glmfit
	%	fn_out = base file name to write plots to
	%   
	%Output:
	%	corrs = a vector of correlation coefficients between simulated firing rate, and actual firing rate
	%		 (a smoothed version of spike train)
	%  
	%Test code:
	%	%Load test preprocessed data
	%	fn_out = './worksheets/09_12_2014/plots/test';
	%	pre = load('./testdata/test_preprocess_spline_short.mat');
	%	const = 'on';
	%	nK_sp = 50; 
	%	nK_pos = 10;
	%	dt_sp = 0.002;
	%	dt_pos = 0.05;
	%	data = filters_sp_pos(pre.processed, nK_sp, nK_pos, dt_sp, dt_pos);
	%	model = MLE_glmfit(data, const);
	%	corrs = glm_predict(pre.processed, data, model, fn_out);

	nTrains = 20;
	nU = size(data.X,1);
	N = size(data.X,2);
	nK_sp = length(data.sp_hist);
	corrs = zeros(nU, 1);
	actualTrain = data.y';
	simTrain = zeros(N, nU);
	binsize = processed.binsize;

	%Make a Gaussian filter to smooth estimates
	sigma = 0.25;
	%sigma = 0.001;
	sigma = sigma/binsize;
	sz = sigma*3*3;
	x = linspace(-sz/2, sz/2, sz);
	gaussFilter = exp(-x.^2/(2*sigma^2));
	gaussFilter = gaussFilter/sum(gaussFilter);

	%Plot ten seconds worth of data
	t_i = 30;
	%t_f = 33;
	t_f = 50; 
	ii = 1:N;
	tt = ii*binsize;

	for i = 1:nTrains
		i
		simTrain = simTrain + glmsim(processed, model, data);
	end
	simTrain = simTrain/nTrains;

	%for each unit
	for i = 1:nU
		%Smooth this average train, and the original train
		smthfittedrates = conv(simTrain(:,i), gaussFilter, 'same')/processed.binsize;
		smthrates = conv(actualTrain(:,i), gaussFilter, 'same')/processed.binsize;
		%Compute correlation coefficient between the two trains
		actualSimCorr = corrcoef(smthfittedrates, smthrates);
		corrs(i) = actualSimCorr(1,2);
		%Plot the simulated trains, the actual train, and the cursor data used to generate it
		clf;
		%Plot firing rate data
		subplot(2,1,1);
		hold on
		plot(tt, smthrates(ii), tt, smthfittedrates(ii))
		%plot(tt, data.y(idx,ii)/processed.binsize, tt, rho_hat(ii)/processed.binsize);
		xlim([t_i, t_f])
		legend('Actual', 'GLM')
		xlabel('time (s)')
		ylabel('estimated firing rate (Hz)')
		title(['Unit: ' processed.unitnames{i} ' correlation: ' num2str(corrs(i))]);

		%Plot cursor data
		subplot(2,1,2);
		plot(tt, data.torque(ii,1), tt, data.torque(ii,2))
		xlim([t_i, t_f])
		legend('RU', 'FE')
		xlabel('time (s)')		

		%Save plot
		saveplot(gcf, [fn_out '_unit_' processed.unitnames{i} '_pred.eps'], 'eps', [6 6]);
		%save fig
		saveas(gcf, [fn_out '_unit_' processed.unitnames{i} '_pred.fig'])
	end