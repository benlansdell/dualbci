function [model, intermediates] = MLE_SD(data, const, fn_out)
	%Fit GLM to spike data from blackrock recording file for each unit above a specified threshold using steepest descent
	%     
	%Input:
	%	data = covariate data output structure from any function in ./models
	%	const = (optional, default = 'on') whether to fit a constant term to the model or not, 'on' or 'off'
	%	fn_out = (optional, default = none) if provided then will plot the (2D) likelihood surface as a function of two parameters
	%   
	%Output:
	%	model is a structure containing the following fields:
	%		b_hat = [nU x (nK + 1)] array with spikes from all channels binned according to binsize. nB = no. bins, nU = no. units.
	%			Note: if a constant term is not fit, a column of zeros is appended to b_hat to make dimensions consistent
	%		dev = [nU x 1] cell array listing deviance of each unit's fit
	%		stats = [nU x 1] cell array listing fitting statistics output from glmfit
	%	intermediates is an array of intermediate estimate values of b_hat at each iteration
	%
	%Test code:
	%	const = 'off';
	%	nK_sp = 0; 
	%	nK_pos = 1;
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_short.mat');
	%	data = filters_sp_pos(pre.processed, nK_sp, nK_pos);
	%	[model, intermediates] = MLE_SD(data, const);

	if (nargin < 2) const = 'on'; end
	if (nargin < 3) fn_out = ''; end

	nU = size(data.y,1);
	nK = size(data.X,3);
	if strcmp(const, 'on')
		model.b_hat = zeros(nU, nK+1);
	else
		model.b_hat = zeros(nU, nK);
	end
	model.dev = cell(nU,1);
	model.stats = cell(nU,1);

	b0 = -1*ones(size(model.b_hat))';
	%b0= [-3.2023; 0.5575; -0.5688];
	%b0= [0.5575; -0.5688];
	%b0 = zeros(2,1);

	%For each unit, fit a GLM to the torque data
	for idx=1:nU 
		[b, dev] = glmfit_SD(squeeze(data.X(idx,:,:)),data.y(idx,:), const, b0);
		%Extract filters fitted...
		model.b_hat(idx,:) = b;	
		model.dev{idx} = dev;
		model.stats{idx} = 0;
	end
	if ~strcmp(const, 'on')
		model.b_hat = [zeros(nU, 1), model.b_hat];
	end
	intermediates = 0;

	%model2 = MLE_glmfit(data, const);
	%Plot the difference between them
	%b_hat1 = model.b_hat(1,:)
	%b_hat2 = model2.b_hat(1,:)
	%if length(fn_out) > 0
	%	%find also the glmfit estimates...
	%	validateFit(squeeze(data.X(1,:,:)),data.y(1,:),b0, b_hat1, b_hat2)
	%end
end

function [b, dev] = glmfit_SD(X, y, const, b0)
	%X is a matrix of data for that unit
	%y is the set of observations for that unit
	if strcmp(const, 'on')
		X = [ones(size(X,1),1), X];		
	end

	db = 0.1;
	N = size(X, 1);
	nIter = 5e7;
	nB = size(X, 2);
	tolfun = 1e-9;
	opts = optimset('Display', 'iter', 'MaxIter', 4000, 'GradObj', 'on', 'Hessian', 'on', 'TolFun', tolfun);
	[b, fval, exitflag, output] = fminunc(@(b) logli(X,y,b), b0, opts);

	%Dbtol = 0.001;
	%Dbmax = 1;
	%Iterate until change in b is very small, or max iterates is reached
	%for idx = 1:nIter
	%	Db = db*dl(X, y, b);
	%	Db = min(Db, Dbmax);
	%	b = b + Db
	%	ll = l(X,y,b)
	%	if norm(Db) < Dbtol
	%		break
	%	end
	%	pause
	%end
	%idx

	%Plot

	%Compute deviance
	dev = 0;
end

function [ll, d, h] = logli(X,y,b_hat)
	ll = -l(X,y,b_hat);
	d = -dl(X,y,b_hat);
	h = -hessian(X,y,b_hat);
end

function ll = l(X,y,b_hat)
	%Log likelihood
	ll = y*X*b_hat-sum(exp(X*b_hat)+log(y+0.00001)');
end

function d = dl(X,y,b_hat)
	%Derivative of log likelihood with respect to each covariate b
	%d = zeros(size(b_hat));
	%for j = 1:length(d)
	%	d(j) = y*X(:,j)-X(:,j)'*exp(X*b_hat);
	%end
	b_hat;
	d = (y*X)' - X'*exp(X*b_hat);
end

function h = hessian(X,y,b_hat)
	%h = -X'*X'*exp(X*b_hat);
	nB = length(b_hat);
	h = zeros(nB);
	for i = 1:nB
		for j = 1:nB
			h(i,j) = -sum(X(:,i).*X(:,j).*exp(X*b_hat));
		end
	end
end

function validateFit(X,y,b0, b_hat1, b_hat2)
	nX = 100;
	brange1 = linspace(-25,25,nX);
	brange2 = linspace(-3,5,nX);
	ll = zeros(nX, nX);
	%Compute log-likelihood as a function of 
	for i = 1:nX
		i;
		for j = 1:nX
			b = [brange1(i), brange2(j)]';
			%ll(i,j) = log(max(-l(X,y,b), 0.001));
			ll(i,j) = -l(X,y,b);
		end
	end
	%Plot likelihood function and the estimates for its maximum
	clf
	pcolor(brange1, brange2, ll);
	colorbar
	hold on
	plot(b0(1), b0(2), 'k.')
	plot(b_hat1(2), b_hat1(3), 'r.')
	plot(b_hat2(2), b_hat2(3), 'g.')
	legend('', 'b0', 'fminunc', 'IRLS')

end
