function [y, tspks, rho, dev] = glmsim_network_rescale(processed, model, data)
	%Compute deviance of a set of data points given a set of fitted coefficients. The deviance is given by:
	%
	%	D(y,mu) = 2\sum y_i ln (y_i / \mu_i) - y_i + \mu_i
	%
	%Usage:
	%	[y, tspks, rho, dev] = glmsim_network_rescale(processed, model, data)
	%     
	%Input:
	%	processed = structure output from preprocess() function
	%	model = a structure of fit coefficients from MLE_glmfit
	%	data = a structure of stimulus and spike history data from ./models
	%   
	%Output:
	%	y = a matrix of simulated spike trains given cursor data for each unit
	%	tspks = 
	%  
	%Test code:
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_short24.mat');
	%	%Make 3x3 network for testing
	%	p = pre.processed; p.binnedspikes = p.binnedspikes(:,1:3); p.unitnames = p.unitnames(1:3);
	%	const = 'on';
	%	nK_sp = 75; 
	%	nK_pos = 6;
	%	dt_sp = p.binsize;
	%	dt_pos = 1/50;
	%	data = filters_sp_pos_network(p, nK_sp, nK_pos, dt_sp, dt_pos);
	%	model = MLE_glmfit_network(data, const);
	%	trains = glmsim_network_rescale(p, model, data);

	bs = processed.binsize;
	nU = size(data.k,2);
	N = size(data.X,1);
	nK_sp = length(data.k{1,2});
	y = zeros(N, nU);
	rho = zeros(N, nU);
	dev = zeros(1,nU);

	%indices of spike history filter
	sp_indices = (1:(nU*nK_sp))+1;
	nspi = length(sp_indices);
	%we're simulating new spikes, so we create a new spike history filter
	sp_hist = zeros(nU,nK_sp);
	%time of spikes
	elecs = cell(1,nU);
	spikemuas = struct('times', elecs);
	for idx=1:nU
		spikemuas(idx).times = [0];    
	end
	%for all units
	%extract all the unit's filters
	b_hat = zeros(size(model.b_hat));
	k_sp = zeros(nU, nspi);
	mu = zeros(nU, N);
	sp_hist = zeros(nspi, 1);
	for i = 1:nU
		b_hat(i,:) = model.b_hat(i,:);
		%and the spike history filter in particular
		k_sp(i,:) = b_hat(i,sp_indices);
		%then set the spike history filter coefficients to zero, since we're generating new spikes and not using data.X's spike history
		b_hat(i,sp_indices) = 0;
		%compute the component of mu = e^eta that comes from the remaining filters used
		mu(i,:) = glmval(b_hat(i,:)', squeeze(data.X(:,:)), 'log');
	end

	%start clocks and draw initial spike times
	clocks = zeros(nU,1);
	tau = exprnd(ones(nU,1));
	%for each time point
	for j = 1:N
		%compute mu that incorporates the current spike history
		mu_sp = mu(:,j).*exp(k_sp*sp_hist);
		%advance the clocks
		clocks = clocks + mu_sp;
		%mark which cells spiked
		yij = clocks >= tau;
		y(j,:) = yij;
		rho(j,:) = mu_sp;
		%for each unit
		for i = 1:nU
			%If there's a spike, add the time to spikemuas
			if (yij(i) > 0)
				sptime = bs*j;%+0.005*(1:yij(i))');
				display(['t = ' num2str(sptime) ' spike! unit = ' num2str(i)])
				spikemuas(i).times = [spikemuas(i).times; sptime];
				%reset this unit's clock
				clocks(i) = 0;
				%draw another spike time for this unit
				tau(i) = exprnd(1);
			end
			dev(i) = dev(i)-2*yij(i)*log(mu_sp(i))-2*yij(i)+2*mu_sp(i);
			%then update the spike history filter
			sph = reshape(sp_hist, nK_sp, nU);
			sph = [sph(2:end,:); yij'];
			sp_hist = reshape(sph, [], 1);
			%sp_hist = [sp_hist(2:end), yij];
		end
	end
	tspks = spikemuas;
