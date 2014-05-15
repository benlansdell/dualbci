function [maxR2 sumR2] = regression_shoham_nev(nevfile, fn_out, kernellength, polydeg)
	%regression_shoham_nev	Function to fit the following model to spike and torque data:
	%			lambda_i(t) = \lambda_0 + \sum_j^N k_j^1 . x_i^1(t+\tau+jh) + \sum_j^N k_j^2 . x_i^2(t+\tau+jh)
	%		That is, it fits a linear filter to the torque data. In the above formula, time-step size h
	%		is given by parameter binsize; max value of N is given by kernellength; \tau is given by offset.
	%		Program fits all models with kernels of length between 1 and N, meaning that for each single-unit N models
	%		will be fit. Both smoothing filters are applied before resampling
	%
	%		Usage:
	%			[maxR2, sumR2] = regression_shoham_nev(nevfile, fn_out, kernlelength. polydeg)
	%
	% 		Input:
	%			nevfile = file to process
	%			fn_out = base name for output plots
	%			kernellength = number of points in linear filter to fit (recommend 4)
	%			polydeg = degree of polynomial for N part of model (recommend 4)
	%		
	%		Output:
	%			maxR2 = maximum R^2 value for linear fit over all units
	%			sumR2 = sum of all R^2 values for every linear fit, over all units
	%
	%		Test code:
	%			nevfile = './testdata/20130117SpankyUtah001.nev';
	%			kernellength = 4;
	%			polydeg = 4;
	%			fn_out = './worksheets/LNmodel/plots/20130117SpankyUtah001';
	%			regression_shoham_nev(nevfile, fn_out, kernellength. polydeg);
	
	%Threshold firing reate below which we ignore that unit
	threshold = 5;
	binsize = 0.05;
	samplerate = 1/binsize;
	offset = 0;
	kl = kernellength;

	sigma_fr = 0.25;
	sz = 30;
   	sigma_fr = sigma_fr*samplerate;
	sz = sigma_fr*3*2;
    x = linspace(-sz/2, sz/2, sz);
    gaussFilter = exp(-x.^2/(2*sigma_fr^2));
   	gaussFilter = gaussFilter/sum(gaussFilter);

	%Preprocess spike and torque data
	[bs rates torque dtorque ddtorque unitnames] = preprocess_shoham_nev(nevfile, fn_out, binsize, threshold, offset);
	nU = length(unitnames);

    %%%%%%%%%%%%
    %Fit models%
    %%%%%%%%%%%%
    %Produce shifted torque data to fit model to
    laggedtorque = [];
    sumR2 = 0;
    for i=1:2
		for j=1:kl
	    	laggedtorque = [laggedtorque, torque(j:end-kl+j,i)];
    	end
    end
    rates = rates(1:end-kl+1,:);
    szRates = size(rates,1);
    fittedrates = zeros(nU, szRates);
    LNfittedrates = zeros(nU, szRates);
    %Do fitting for a range of kernel lengths
    models = {};
    lambda_0 = zeros(nU);
    kernelRU = zeros(nU, kl);
    kernelFE = zeros(nU, kl);
    r2 = zeros(nU);
	display(['Fitting model with kernel length: ' num2str(i)])
	for j=1:nU
		%Estimate parameters as a linear model
		X = [laggedtorque(:,1:kl), laggedtorque(:,(kl+1):(2*kl))];
		Y = rates(:,j);
		model = LinearModel.fit(X,Y);
		models{j} = model;
		r2(j) = model.Rsquared.Ordinary;
		sumR2 = sumR2 + r2(j);
		coeffs = model.Coefficients.Estimate;
		lambda_0(j) = coeffs(1);
		kernelRU(j,1:kl) = coeffs(2:kl+1);
		kernelFE(j,1:kl) = coeffs(kl+2:end);
		%Do the same thing but do the least squares fit ourselves...
		%%Add a constant intercept term
		X = [ones(szRates, 1), X];
		%%Solve normal equations
		beta = (X'*X)\(X'*Y);
    	%%...these give the same result, but it's easy to get the predicted Y values from here
    	fittedrates(j,:) = X*beta;
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
    	subplot(1,2,1);
		hold on
		plot((1:kl)*binsize, squeeze(kernelRU(i,1:kl)));
		keys1 = {keys1{:}, ['R^2=' num2str(r2(i))]};
		keys2 = {keys2{:}, ['lambda_0=' num2str(lambda_0(i))]};
		title(['Unit: ' unitnames(i)])
		xlabel('time (s)')
		ylabel('k [Radial-Ulnar]')
		legend(keys1, 'Location', 'NorthOutside')
		subplot(1,2,2)
		hold on
		plot((1:kl)*binsize, squeeze(kernelFE(i,1:kl)));
		legend(keys2, 'Location', 'NorthOutside')
		xlabel('time (s)')
		ylabel('k [Flexion-Extension]')
		saveplot(gcf, [fn_out '_kernel_unit_' unitnames{i} '_R2_' num2str(r2(i)) '.eps'])
    end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Fit non-linear part of model%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%For each neuron
	for i = 1:nU
		%Put linear predictions into bins
		nbins = 50;
		umin = min(fittedrates(i,:))-0.001;
		umax = max(fittedrates(i,:));
		bins = cell(1,nbins);
		for j=1:szRates
			u = fittedrates(i,j);
			count = rates(j,i);
			bin = ceil((u-umin)/(umax-umin)*nbins);
			bins{bin} = [bins{bin}; count];
		end
		%Compute mean for each bin, and variance
		meanspikes = zeros(nbins,1);
		varspikes = zeros(nbins,1);
		for j=1:nbins
			meanspikes(j)=mean(bins{j});
			stdspikes(j)=std(bins{j});
		end
		%Fit polynomial to mean firing rate
		upts = linspace(umin, umax, nbins);
		upts = upts(~isnan(meanspikes));
		pp = polyfit(upts', meanspikes(~isnan(meanspikes)), polydeg);
		%Plot
		figure;
		hold on
		plot(upts, polyval(pp, upts), upts, meanspikes(~isnan(meanspikes)), '.');
		%errorbar(upts, meanspikes(~isnan(meanspikes)), stdspikes(~isnan(meanspikes)), '.');
		xlabel('Linear prediction')
		ylabel('Mean spikes/s. Error bars = std')
		saveplot(gcf, [fn_out '_expectedfiringrate_unit_ ' unitnames{i} '.eps']);
		%Compute predicted firing rate (lambda) using our LN model 
		LNfittedrates(i,:) = polyval(pp, fittedrates(i,:));
		%Put observed spikes/s vs LN predicted rate into bins

		%Compute variance per bin, plot, show not Poisson (or are they?)

	end 

	%%%%%%%%%%%%%%%%%%%%
	%Plot filtered data%
	%%%%%%%%%%%%%%%%%%%%
	for j=1:nU
		subplot(1,4,1)
		%Plot histogram of filtered data
		hist(squeeze(fittedrates(j,:)));
		title(['Unit: ' unitnames{j} ' kernel length: ' num2str(kl)])
		xlabel('predicted lambda')
		%Plot scatter plot of filtered data vs firing rate
		subplot(1,4,2)
		%smoothhist2D([reshape(squeeze(fittedrates(i,j,:)),szRates,1), -rates(:,j)], 5, [100 100], 0.05);
		smoothhist2D([reshape(squeeze(fittedrates(j,:)),szRates,1), rates(:,j)], 5, [100 100], 0.05);
		axis xy;
		xlabel('linear prediction (lambda)')
		ylabel('estimated instantaneous rate (hat lambda)')
		subplot(1,4,3)
	   	t = 30;
    	times = (1:(t*samplerate))*binsize;
		plot(times, rates(1:(t*samplerate), j), times, squeeze(fittedrates(j,1:(t*samplerate))))
		title('Linear model');
		xlabel('time (s)')
		ylabel('spikes/s')
		legend('Actual', 'Estimation')
		subplot(1,4,4)
    	times = (1:(t*samplerate))*binsize;
		plot(times, conv(rates(1:(t*samplerate), j), gaussFilter, 'same'), times, conv(squeeze(LNfittedrates(j,1:(t*samplerate))), gaussFilter, 'same'));
		title('LN model')
		xlabel('time (s)')
		ylabel('spikes/s')
		legend('Actual', 'Estimation')
		saveplot(gcf, [fn_out '_LNmodel_' unitnames{j} '_klength_' num2str(kl) '.eps'], 'eps', [12 3]);
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Compute likelihood of LN model given data%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



end
