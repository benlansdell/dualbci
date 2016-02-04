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

	mu = glmval(b, data.X, 'log');
	nz = y > 0;
	%Return deviance
	dev = 2*sum(y(nz).*log(y(nz)./mu(nz)))-2*sum(y-mu);

end