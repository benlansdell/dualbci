function model = MLE_glmfit_lasso(data, const)
	%Fit GLM to spike data from blackrock recording file for each unit above a specified threshold
	%     
	%Input:
	%	data = covariate data output structure from any function in ./models
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
	%	model = MLE_glmfit_lasso(data, const);

	if (nargin < 2) const = 'on'; end

	nU = size(data.y,1);
	nK = size(data.X,3);
	model.b_hat = zeros(nU, nK+1);
	model.dev = cell(nU,1);
	model.stats = cell(nU,1);

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