function [F, Q, mu] = fit_AR(data, order)
	%Fit an auto-regressive model to torque data of a given order using least squares. That is, fits a model of the form:
	%
	%	[x_2 y_2] ~ [mu_1, mu_2] + [x_{1}, x_{0}, y_{1}, y_{0},...]x[F_11 F_21] + [\epsilon_i^1, \epsilon_i^2]
	%	[x_3 y_3]				   [x_{2}, x_{1}, y_{2}, y_{1},...] [F_12 F_22]
	%	[...    ]				   [...                           ] [F_13 F_23]
	%											  	                [...      ]
	%
	%Usage:
	%	[F, Q, mu] = fit_AR(data, order)
	%     
	%Input:
	%	data = [N x 2] matrix where N is the number of data points. Could be cursor position data, 
	%			or cursor (direction, velocity) data...
	%	order = (optional, default = 1) Order of AR model to fit
	%   
	%Output:
	%	F = fitted coefficient matrix
	%	Q = covariance matrix giving errors
	%	mu = mean data
	%  
	%
	%Test code:
	%	%Load test preprocessed data
	%	pre = load('./testdata/test_preprocess_spline_short.mat');
	%	order = 1;
	%	data = pre.processed.torque;
	%	[F, Q, mu] = fit_AR(data, order);

	if (nargin < 2) order = 1; end
	F = zeros(order*2, 2);
	mu = zeros(1, 2);
	residuals = zeros(size(data,1)-order, size(data,2));

	for idx = 1:2
		Y = data(order+1:end,idx);
		%Set up matrices
		c = data(order:end-1, 1);
		r = data(order:-1:1,1);
		X_RU = toeplitz(c,r);
		c = data(order:end-1,2);
		r = data(order:-1:1,2);
		X_FE = toeplitz(c,r);
		%Fit a constant term
		X = [ones(size(X_RU,1),1), X_RU, X_FE];
		%Do the least squares fit:
		beta_hat = (transpose(X)*X)\(transpose(X)*Y);
		mu(idx) = beta_hat(1);
		F(1:2*order,idx) = beta_hat(2:end);
		residuals(:,idx) = Y - X*beta_hat;
	end
	%Compute the variance/covariance of residuals, Q.
	Q = cov(residuals(:,1), residuals(:,2));