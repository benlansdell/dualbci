function model = MLE_glmfit_network_groupL1(data, lambda, const)
	%Fit GLM to spike data from blackrock recording file for each unit above a specified threshold
	%
	%Input:
	%	data = covariate data output structure from any function in ./models
	%	lambda = weight for prior term
	%	const = (optional, default = 'on') whether to fit a constant term to the model or not, 'on' or 'off'
	%
	%Output:
	%	model is a structure containing the following fields:
	%		b_hat = [nU x (nK + 1)] array with spikes from all channels binned according to binsize. nB = no. bins, nU = no. units.
	%			Note: if a constant term is not fit, a column of zeros is appended to b_hat to make dimensions consistent
	%		dev = [nU x 1] cell array listing deviance of each unit's fit
	%		stats = [nU x 1] cell array listing fitting statistics output from glmfit
	%		converged = [nU x 1] array listing 1 if the IRLS converged within iteration limit, 0 if not
	%		conditioned = [nU x 1] array listing 1 if the IRLS did not issue an ill-conditioned warning, 0 if it did
	%
	%Test code:
	%	const = 'on';
	%	nK_sp = 30; 
	%	nK_pos = 2;
	%	lambda = 30;
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_60hz_short24.mat');
	%	data = filters_sp_pos_network(pre.processed, nK_sp, nK_pos);
	%	model = MLE_glmfit_network_groupL1(data, lambda, const);

	if (nargin < 3) const = 'on'; end
	nU = size(data.y,1);
	nK = size(data.X,2);
	if strcmp(const, 'on')
		model.b_hat = zeros(nU, nK+1);
		model.mask = zeros(nU, nK+1);
	else
		model.b_hat = zeros(nU, nK);
		model.mask = zeros(nU, nK);
	end
	model.dev = cell(nU,1);
	model.stats = cell(nU,1);
	model.converged = ones(nU,1);
	model.conditioned = ones(nU,1);
	cutoff = 1e-10;
	%For each unit, fit a GLM to the torque data
	display(['Fitting GLM by MLE with group L1 penalty. Fitting ' num2str(nU) ' units.'])
	for idx=1:nU 
		display(['Fitting unit ' num2str(idx)])
		%Mask columns that don't vary... they cannot be estimated.
		%mask = (std(data.X) > cutoff);
		%if strcmp(const, 'off')
		%	m = mask;
		%else
		%	m = [1==1 mask];
		%end
		%model.mask(idx,:) = m;
		%[b, dev] = groupL1fit(data.X(:,mask), data.y(idx,:), data, lambda, const, idx, nU);

		%Initial guess for parameters
		b0 = glmfit(data.X, data.y(idx,:), 'poisson', 'constant', const);
		[b, ll, dev] = groupL1fit(data.X, data.y(idx,:), data, lambda, const, idx, nU, b0);
		%Catch if a warning was raised about badly conditioned matrix
		[warn, warnid] = lastwarn;
		if ~strcmp(warn, '')
	   		switch warnid
        	case 'stats:glmfit:IterationLimit'
        		model.converged(idx) = 0;
        	case 'stats:glmfit:BadScaling'
        		model.conditioned(idx) = 0;
       		end
	    end
	    lastwarn('')
	    %Extract filters fitted...
		%model.b_hat(idx,m) = b;	
		model.b_hat(idx,:) = b;	
		model.dev{idx} = dev;
		model.stats{idx} = {};

	end
	model.logli = ll_network(model, data, 'poisson');
	if ~strcmp(const, 'on')
		model.b_hat = [zeros(nU, 1), model.b_hat]
	end
	display('Done')
end

function [b, ll, dev] = groupL1fit(X, y, data, lambda, const, iU, nU, prs)

	method = 'spg';
	groups = zeros(1,size(X,2)+1)';
	%Change number of groups
	nG = size(data.k,1);
	nP = size(data.k{1,2},2);
	
	%Change lambda accordingly
	lambda = 30;
	gOptions.maxIter = 500;
	gOptions.verbose = 2; % Set to 0 to turn off output
	gOptions.corrections = 10; % Number of corrections to store for L-BFGS methods
	gOptions.norm = 2; % Set to inf to use infinity norm
	options = gOptions;
	for idx = 1:nG
		indices = data.k{idx,2}+1;
		if idx ~= iU
			groups(indices) = idx;
		end
		data.k{idx,2};
	end

	%Initial guess
	%prs = zeros(size(groups));
	lambdaVect = lambda*ones(nG, 1);	
	logli = @(prs) ll_network_L1(prs, data, iU);
	options.method = method;
	b = L1GeneralGroup_Auxiliary(logli,prs,lambdaVect,groups,options);
	%dev = dev_network_L1(prs, data, 'poisson');
	ll = logli(b);

	mu = glmval(b, data.X, 'log')';
	nz = y > 0;
	%Return deviance
	dev = 2*sum(y(nz).*log(y(nz)./mu(nz)))-2*sum(y-mu);

end