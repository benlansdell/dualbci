function [maxR2 sumR2] = regression_nev(nevfile, fn_out, kernellength, binsize, sigma_fr, sigma_trq, offset, verbosity)
	%regression_nev	Function to fit the following model to spike and torque data:
	%			lambda_i(t) = \lambda_0 + \sum_j^N k_j^1 . x_i^1(t+\tau+jh) + \sum_j^N k_j^2 . x_i^2(t+\tau+jh)
	%		That is, it fits a linear filter to the torque data. In the above formula, time-step size h
	%		is given by parameter binsize; max value of N is given by kernellength; \tau is given by offset.
	%		Program fits all models with kernels of length between 1 and N, meaning that for each single-unit N models
	%		will be fit. Both filters are applied before resampling
	%
	% 		Input:
	%			nevfile = file to process
	%			fn_out = base name for output plots
	%			kernellength = (optional, default = 6) max 
	%			binsize = (optional, default = 0.05) size of timebins over which to compute regression
	%			sigma_fr = (optional, default = 0.25) width of gaussian filter to apply to spikes for firing rate. If 0 then no filter applied
	%			sigma_trq = (optional, default = 0.25) width of gaussian filter to apply to torque. If 0 then no filter applied
	%				Note: Note: both sigmas are in units of seconds, and then are scaled according to binsize
	%			offset = (optional, default = 0) number of seconds to add to spike data before comparing with torque
	%			versbose = (optional, default = 0) verbosity level. If above 0 then plot/print extra info
	%		
	%		Output:
	%			maxR2 = maximum R^2 value for linear fit over all units and all fitted kernel lengths
	%			sumR2 = sum of all R^2 values for every linear fit, over all units and all fitted kernel lengths
	%
	%		Test code:
	%			nevfile = './testdata/20130117SpankyUtah001.nev';
	%			binsize = 0.05;
	%			kernellength = 6;
	%			sigma_fr = 0.25;
	%			sigma_trq = 0.25;
	%			offset = 0;
	%			verbosity = 1;
	%			fn_out = './worksheets/tuning/01172013/20130117SpankyUtah001';
	%			regression_nev(nevfile, fn_out, kernellength, binsize, sigma_fr, sigma_trq, offset, verbosity);
	
	%Optional arguments
	if (nargin < 3)	kernellength = 6; end
	if (nargin < 4) binsize = 0.05; end
	if (nargin < 5) sigma_fr = 0.25; end
	if (nargin < 6) sigma_trq = 0.25; end
	if (nargin < 7) offset = 0; end
	if (nargin < 8) verbosity = 0; end
	%Set this to above 0 if want debug info/plots
	%verbosity = 1;
	%Threshold firing reate below which we ignore that unit
	threshold = 5;
	samplerate = 1/binsize;
	%Preprocess spike and torque data
	[binnedspikes rates torque unitnames] = preprocess_nev(nevfile, fn_out, binsize, sigma_fr, sigma_trq, threshold, offset, verbosity);
	nU = length(unitnames);

    %%%%%%%%%%%%
    %Fit models%
    %%%%%%%%%%%%
    %Produce shifted torque data to fit model to
    laggedtorque = [];
    sumR2 = 0;
    for i=1:2
		for j=1:kernellength
	    	laggedtorque = [laggedtorque, torque(j:end-kernellength+j,i)];
    	end
    end
    rates = rates(1:end-kernellength+1,:);
    szRates = size(rates,1);
    fittedrates = zeros(kernellength, nU, szRates);
    %Do fitting for a range of kernel lengths
    models = {};
    lambda_0 = zeros(kernellength, nU);
    kernelRU = zeros(kernellength, nU, kernellength);
    kernelFE = zeros(kernellength, nU, kernellength);
    r2 = zeros(kernellength, nU);
    sumr2 = zeros(kernellength, 1);
    for i=1:kernellength
		display(['Fitting model with kernel length: ' num2str(i)])
		for j=1:nU
			%Estimate parameters as a linear model
			X = [laggedtorque(:,1:i), laggedtorque(:,kernellength+1:kernellength+i)];
			Y = rates(:,j);
			model = LinearModel.fit(X,Y);
			models{i,j} = model;
			r2(i,j) = model.Rsquared.Ordinary;
			sumR2 = sumR2 + r2(i,j);
			coeffs = model.Coefficients.Estimate;
			lambda_0(i,j) = coeffs(1);
			kernelRU(i,j,1:i) = coeffs(2:i+1);
			kernelFE(i,j,1:i) = coeffs(i+2:end);
			%Do the same thing but do the least squares fit ourselves...
			%%Add a constant intercept term
			X = [ones(szRates, 1), X];
			%%Solve normal equations
			beta = (X'*X)\(X'*Y);
	    	%%...these give the same result, but it's easy to get the predicted Y values from here
	    	fittedrates(i,j,:) = X*beta;
	    end
    end
    maxR2 = max(max(r2));
    display(['Done fitting. Max R^2 value: ' num2str(maxR2) ' Sum of all R^2 values: '  num2str(sumR2)])

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   	%Plot kernels for each neuron%
   	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   	colors = get(gca, 'colororder');
    for i=1:nU
    	figure
    	keys1 = {};
    	keys2 = {};
    	subplot(1,2,1)
		for j=1:kernellength
			hold on
			plot((1:j)*binsize, squeeze(kernelRU(j,i,1:j)), 'Color', colors(j,:));
			keys1 = {keys1{:}, ['R^2=' num2str(r2(j,i))]};
			keys2 = {keys2{:}, ['lambda_0=' num2str(lambda_0(j,i))]};
		end
		title(['Unit: ' unitnames(i)])
		xlabel('time (s)')
		ylabel('k [Radial-Ulnar]')
		legend(keys1, 'Location', 'NorthOutside')
		subplot(1,2,2)
		for j=1:kernellength
			hold on
			plot((1:j)*binsize, squeeze(kernelFE(j,i,1:j)), 'Color', colors(j,:));
		end
		legend(keys2, 'Location', 'NorthOutside')
		xlabel('time (s)')
		ylabel('k [Flexion-Extension]')
		saveplot(gcf, [fn_out '_kernel_unit_' unitnames{i} '_maxR2_' num2str(max(r2(:,i))) '.eps'])
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Heat map of goodness-of-fit for all units%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for i=1:kernellength
		subplot(kernellength,1,i)
		gridx = floor(sqrt(nU));
		if (gridx < sqrt(nU)) 
			gridx = gridx + 1;
			gridy = gridx; 
			r2_imputed = [squeeze(r2(i,:)), ones(1,gridx*gridy-nU)];
		else
			r2_imputed = [squeeze(r2(i,:))];
		end
		image(reshape(r2_imputed,gridx,gridy), 'CDataMapping', 'scaled');
		zaxis = [0 1];
		caxis(zaxis);
		set(gca,'Zlim',zaxis,'Ztick',zaxis);
		xlabel('channel');
  		title(['R^2 values for kernels of length ' num2str(i)])
   		colorbar
	end
	saveplot(gcf, [fn_out '_regress_R2_heatmap_maxR2_' num2str(maxR2) '_sumR2_' num2str(sumR2) '.eps'], 'eps', [3 1.5*kernellength]);
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Plot filtered data%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for i=1:kernellength
		for j=1:nU
			subplot(1,3,1)
			%Plot histogram of filtered data
			hist(squeeze(fittedrates(i,j,:)));
			title(['Unit: ' unitnames{j} ' kernel length: ' num2str(i)])
			xlabel('predicted lambda')
			%Plot scatter plot of filtered data vs firing rate
			subplot(1,3,2)
			%smoothhist2D([reshape(squeeze(fittedrates(i,j,:)),szRates,1), -rates(:,j)], 5, [100 100], 0.05);
			smoothhist2D([reshape(squeeze(fittedrates(i,j,:)),szRates,1), rates(:,j)], 5, [100 100], 0.05);
			axis xy;
			xlabel('linear prediction (lambda)')
			ylabel('estimated instantaneous rate (hat lambda)')
			subplot(1,3,3)
	    	t = 50; unit = 18;
    		times = (1:(t*samplerate))*binsize;
			plot(times, rates(1:(t*samplerate), j), times, squeeze(fittedrates(i,j,1:(t*samplerate))))
			xlabel('time (s)')
			ylabel('spikes/s')
			legend('Actual', 'Estimation')
			saveplot(gcf, [fn_out '_regress_filtered_unit_' unitnames{j} '_klength_' num2str(i) '.eps'], 'eps', [9 3]);
		end
	end
end
