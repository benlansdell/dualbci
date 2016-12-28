function model = MLE_glmfit_L1(data, lambda, const)
	%Fit GLM to spike data from blackrock recording file for each unit above a specified threshold
	%     
	%Input:
	%	data = covariate data output structure from any function in ./models
	%	lambda = weight for L1 prior term
	%	const = (optional, default = 'on') whether to fit a constant term to the model or not, 'on' or 'off'
	%   
	%Output:
	%	model is a structure containing the following fields:
	%		b_hat = [nU x (nK + 1)] array with spikes from all channels binned according to binsize. nB = no. bins, nU = no. units.
	%			Note: if a constant term is not fit, a column of zeros is appended to b_hat to make dimensions consistent
	%		dev = [nU x 1] cell array listing deviance of each unit's fit
	%		stats = [nU x 1] cell array listing fitting statistics output from glmfit
	%
	%Test code:
	%	const = 'on';
	%	nK_sp = 6; 
	%	nK_pos = 6;
	%	%Load test preprocessed data
	%	%pre = load('./testdata/test_preprocess_spline_short.mat');
	%	pre = load('./testdata/test_preprocess_bs_50ms_short.mat');
	%	data = filters_sp_pos(pre.processed, nK_sp, nK_pos);
	%	model = MLE_glmfit_L1(data, const);

	if (nargin < 3) const = 'on'; end

	nU = size(data.y,1);
	nK = size(data.X,3);
	model.b_hat = zeros(nU, nK+1);
	model.dev = cell(nU,1);
	model.stats = cell(nU,1);

	method = 'spg';
	groups = zeros(1,size(nK));
	
	%Change number of groups
	nG = nU;
	nP = nK;	
	gOptions.maxIter = 500;
	gOptions.verbose = 2; % Set to 0 to turn off output
	gOptions.corrections = 10; % Number of corrections to store for L-BFGS methods
	gOptions.norm = 2; % Set to inf to use infinity norm
	options = gOptions;
	istart = size(gg.kt,1)*size(gg.kt,2)+1+size(gg.ih,1);
	for idx = 1:nG
		indices = istart + (((idx-1)*nP+1):(idx*nP));
		groups(indices) = idx; 
	end
	lambdaVect = lambda*ones(nG, 1);	
	options.method = method;
	w = L1GeneralGroup_Auxiliary(@Loss_GLM_logli,prs,lambdaVect,groups,options);

	%For each unit, fit a GLM to the torque data
	for idx=1:nU
		X = squeeze(data.X(idx,:,:));
		%If a constant term is to be added, then add a column of ones to data
		if isequal(const, 'on')
			X = [ones(size(X, 1), 1), X];
		end
		[b, fitinfo] = lassoglm(X,data.y(idx,:),'poisson');
		%Extract filters fitted...
		model.b_hat(idx,:) = b;	
		model.dev{idx} = fitinfo.Deviance;
		model.stats{idx} = fitinfo;
	end
	if ~isequal(const, 'on')
		model.b_hat = [zeros(nU, 1), model.b_hat]
	end