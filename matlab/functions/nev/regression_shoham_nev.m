function [r2 LNr2 loglikelihoods, lp] = regression_shoham_nev(nevfile, fn_out, kernellength, polydeg)
	%regression_shoham_nev	Function to fit the following model to spike and torque data:
	%			lambda_i(t) = \lambda_0 + \sum_j^N k_j^1 . x_i^1(t+\tau+jh) + \sum_j^N k_j^2 . x_i^2(t+\tau+jh)
	%		That is, it fits a linear filter to the torque data. In the above formula, time-step size h
	%		is given by parameter binsize; N is given by kernellength.
	%		
	%		Then fits a non-linear model to the filtered data
	%
	%		Matches methods used by Shoham et al 2005 in fitting LN model. 50ms time bins.
	%
	%		Usage:
	%			[r2 LNr2 loglikelihoods lp] = regression_shoham_nev(nevfile, fn_out, kernlelength, polydeg)
	%
	% 		Input:
	%			nevfile = file to process
	%			fn_out = base name for output plots
	%			kernellength = number of points in linear filter to fit (recommend 4)
	%			polydeg = degree of polynomial for N part of model (recommend 4)
	%		
	%		Output:
	%			r2 = R^2 values for linear model with data
	%			LNr2 = R^2 values for LN model with data
	%			loglikelihood = log likelihood of data given fitted model, assuming poisson neurons
	%			lp = penalized log likelihood. Subtract #params/2*log(#datapoints) (Bayesian information criterion)
	%
	%		Test code:
	%			nevfile = './testdata/20130117SpankyUtah001.nev';
	%			kernellength = 4;
	%			polydeg = 4;
	%			fn_out = './worksheets/LNmodel/plots/20130117SpankyUtah001';
	%			[r2 LNr2 ll lp] = regression_shoham_nev(nevfile, fn_out, kernellength, polydeg);
	
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
    smthLNfittedrates = zeros(nU, szRates);
    smthfittedrates = zeros(nU, szRates);
    smthrates = zeros(nU, szRates);
    %Do fitting for a range of kernel lengths
    models = {};
    %lambda_0 = zeros(nU);
    kernelRU = zeros(nU, kl);
    kernelFE = zeros(nU, kl);
    r2 = zeros(nU,1);
    LNr2 = zeros(nU,1);
	for j=1:nU
		%Estimate parameters as a linear model
		X = [laggedtorque(:,1:kl), laggedtorque(:,(kl+1):(2*kl))];
		Y = rates(:,j);
		%model = LinearModel.fit(X,Y);
		model = LinearModel.fit(X,Y, 'Intercept', false);
		models{j} = model;
		r2(j) = model.Rsquared.Ordinary;
		sumR2 = sumR2 + r2(j);
		coeffs = model.Coefficients.Estimate;
		%lambda_0(j) = coeffs(1);
		%kernelRU(j,1:kl) = coeffs(2:kl+1);
		%kernelFE(j,1:kl) = coeffs(kl+2:end);

		kernelRU(j,1:kl) = coeffs(1:kl);
		kernelFE(j,1:kl) = coeffs(kl+1:end);
		%Do the same thing but do the least squares fit ourselves...
		%%Add a constant intercept term (or not)
		%X = [ones(szRates, 1), X];
		%%Solve normal equations
		beta = (X'*X)\(X'*Y);
    	%%...these give the same result, but it's easy to get the predicted Y values from here
    	fittedrates(j,:) = X*beta;
    end
    maxR2 = max(max(r2));
    display(['Done fitting linear model. Max R^2 value: ' num2str(maxR2) ' Mean R^2 value: '  num2str(sumR2/nU)])

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Fit non-linear part of model, plot kernels for each unit, plot filtered data%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for i=1:nU
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
			varspikes(j)=std(bins{j});
		end
		%Fit polynomial to mean firing rate
		upts = linspace(umin, umax, nbins);
		upts = upts(~isnan(meanspikes));
		pp = polyfit(upts', meanspikes(~isnan(meanspikes)), polydeg);
		%Plot non-linear fit to mean firing rate
		figure;
		subplot(2,4,5)
		hold on
		plot(upts, polyval(pp, upts), upts, meanspikes(~isnan(meanspikes)), '.');
		%errorbar(upts, meanspikes(~isnan(meanspikes)), stdspikes(~isnan(meanspikes)), '.');
		xlabel('linear prediction');
		ylabel('mean spikes/s');
		%Compute predicted firing rate (lambda) using our LN model 
		LNfittedrates(i,:) = polyval(pp, fittedrates(i,:));
		%Make sure always positive
		%LNfittedrates(i,:) = LNfittedrates(i,:).*(LNfittedrates(i,:)>=0);
		%Smooth for plotting
		smthLNfittedrates(i,:) = conv(LNfittedrates(i,:), gaussFilter, 'same');
		smthfittedrates(i,:) = conv(fittedrates(i,:), gaussFilter, 'same');
		smthrates(i,:) = conv(squeeze(rates(:,i)), gaussFilter, 'same');
		subplot(2,4,8);
		scatter(smthrates(i,:),smthLNfittedrates(i,:));
		xlabel('LN mean firing rate');
		ylabel('spikes/s');
		LNr2(i) = corr(smthrates(i,:)',smthLNfittedrates(i,:)');
		title(['R^2= ' num2str(LNr2(i))]);
		
		%Plot linear filters fitted
    	keys1 = {};
    	keys2 = {};
    	subplot(2,4,3);
		hold on
		plot((0:kl-1)*binsize, squeeze(kernelRU(i,1:kl)));
		title(['R^2=' num2str(r2(i))]);
		xlabel('time (s)')
		ylabel('k [RU]')
		subplot(2,4,4)
		hold on
		plot((0:kl-1)*binsize, squeeze(kernelFE(i,1:kl)));
		xlabel('time (s)')
		ylabel('k [FE]')
		%Plot histogram of filtered data
		subplot(2,4,1)
		hist(squeeze(fittedrates(i,:)));
		title(['Unit: ' unitnames{i}]);
		xlabel('linear prediction');
		%Plot scatter plot of filtered data vs firing rate
		subplot(2,4,2);
		scatter(reshape(squeeze(fittedrates(i,:)),szRates,1), rates(:,i));
		xlabel('linear prediction');
		ylabel('spikes/s');
		subplot(2,4,6);
	   	t = 100;
	   	dt = 30;
    	times = ((t*samplerate):((t+dt)*samplerate))*binsize;
		%plot(times, rates((t*samplerate):((t+dt)*samplerate), i), times, smthfittedrates(i,(t*samplerate):((t+dt)*samplerate)));
		plot(times, smthrates(i,(t*samplerate):((t+dt)*samplerate)), times, smthfittedrates(i,(t*samplerate):((t+dt)*samplerate)));
		title('linear model');
		xlabel('time (s)');
		ylabel('spikes/s');
		legend('Actual', 'Model', 'Location', 'NorthOutside');
		subplot(2,4,7);
		%plot(times, rates((t*samplerate):((t+dt)*samplerate),i), times, smthLNfittedrates(i,(t*samplerate):((t+dt)*samplerate)));
		plot(times, smthrates(i,(t*samplerate):((t+dt)*samplerate)), times, smthLNfittedrates(i,(t*samplerate):((t+dt)*samplerate)));
		title('LN model');
		xlabel('time (s)');
		ylabel('spikes/s');
		legend('Actual', 'Model', 'Location', 'NorthOutside');
		saveplot(gcf, [fn_out '_LNmodel_' unitnames{i} '_klength_' num2str(kl) '.eps'], 'eps', [12 6]);
	end
    display(['Done fitting LN model. Max R^2 value: ' num2str(max(max(LNr2))) ' Mean R^2 value: '  num2str(sum(LNr2)/nU)])
	figure
	subplot(1,2,1)
    hist(r2');
    xlabel('R^2 value (linear model)');
	subplot(1,2,2)
	hist(LNr2'.^2);
	xlabel('R^2 value (LN model)');
    saveplot(gcf, [fn_out '_r2vals.eps']);

    maxspks = max(max(bs));
    pn_given_lambda = zeros(nU, nbins, max(4,maxspks+1));

    loglikelihoods = zeros(nU);

	for i=1:nU
    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    	%How Poisson is each neuron?%
    	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   		%Put observed spikes/s vs LN predicted rate into bins
		umin = min(LNfittedrates(i,:))-0.001;
		umax = max(LNfittedrates(i,:));
		bins = cell(1,nbins);
		for j=1:szRates
			u = LNfittedrates(i,j);
			count = rates(j,i);
			bin = ceil((u-umin)/(umax-umin)*nbins);
			bins{bin} = [bins{bin}; count];
			pn_given_lambda(i,bin,bs(j,i)+1) = pn_given_lambda(i,bin,bs(j,i)+1)+1;
		end
		%Compute mean for each bin, and variance
		meanspikes = zeros(nbins,1);
		varspikes = zeros(nbins,1);
		for j=1:nbins
			meanspikes(j)=mean(bins{j});
			varspikes(j)=std(bins{j}).^2/samplerate;
		end

		figure
		subplot(2,2,1)
		hold on;
		bins = 1:nbins;
		plot(bins/nbins*(umax-umin)+umin, varspikes, '.');
		plot(0:umax,0:umax, 'r');
		xlabel('Expected rate');
		ylabel('Conditional variance');
		title(['Unit: ' unitnames{i}])
		subplot(2,2,3)
		hold on;
%		plot(bins/nbins*(umax-umin)+umin, meanspikes);
		plot(meanspikes, varspikes, '.');
		plot(0:umax,0:umax, 'r');
		xlabel('Conditional expected rate');
		ylabel('Conditional variance');
		title(['Unit: ' unitnames{i}])

		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%Compute likelihood of LN model given data%
		%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%Estimate parameters of NG distribution for each neuron from cursor data:
		%Compute observed P(N|lambda) distributions for each neuron
		for j=1:nbins
			pn_given_lambda(i,j,:) = pn_given_lambda(i,j,:)/sum(squeeze(pn_given_lambda(i,j,:)));
		end
		%Fit sigma to these distributions

		%Plot empirical distribution
		%Test how well they fit
		subplot(2,2,2);
		plot(bins, squeeze(pn_given_lambda(i,:,1)), bins, squeeze(pn_given_lambda(i,:,2)), bins, squeeze(pn_given_lambda(i,:,3))...
			, bins, squeeze(pn_given_lambda(i,:,4)));
		xlabel('\lambda')
		ylabel('Pr(N|\lambda)')
		legend('P(N=0|\lambda)', 'P(N=1|\lambda)', 'P(N=2|\lambda)', 'P(N=3|\lambda)', 'Location', 'NorthOutside')
		saveplot(gcf, [fn_out '_pnlambda_unit_' unitnames{i} '.eps'], 'eps', [6 6]);
		%use those estimates to compute likelihood of data, assuming each time point is independent
		%
		%Compute poisson likelihood

		for j=1:szRates
			lambda = LNfittedrates(i,j);
			pj = poisson(lambda, binsize, bs(j,i));
			loglikelihoods(i)=loglikelihoods(i)+real(log(pj));
		end
	end
	loglikelihood = sum(sum(loglikelihoods));
	display(['Log likelihood of model: '  num2str(loglikelihood)]);
	nparams = polydeg+2*kernellength;
	ndatapts = szRates;
	lp = loglikelihood - nparams*log(ndatapts)/2;
	display(['Penalized likelihood: ' num2str(lp)]);
end

function pn = pn_lambda(N, lambda, sigma)
	if (N == 0)
		pn = alambda(lambda, sigma).*exp(-lambda.^2./(sigma.^2)/2);
	else
		pn = blambda(lambda, sigma).*exp(-(lambda-N).^2./(sigma.^2)/2);
	end
end

function a=alambda(lambda, sigma)
	Ns=1:1000;
	b = blambda(lambda, sigma);
	a = (1-b*sum(exp(-(lambda-Ns).^2./(sigma.^2)/2)))*exp(lambda^2/2/sigma^2);

end

function b=blambda(lambda, sigma)
	Ns=1:1000;
	b=lambda/sum(Ns.*exp(-(lambda-Ns).^2/(sigma^2)/2));
end

function p=poisson(lambda, dt, n)
	if n < 0
		display('n less than zero')
	end
	p = ((lambda*dt)^n)*exp(-lambda*dt)/factorial(n);
end
