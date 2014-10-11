function output = glm_decode(processed, data, model, F, Q, mu, fn_out)
	%Predict cursor position using recursive Bayesian filter for GLM given a model 
	%	compares to cursor data.
	%
	%Usage:
	%	output = glm_decode(processed, data, model, F, Q, mu, fn_out)
	%     
	%Input:
	%	processed = a structure output from a preprocess function
	%	data = a structure of stimulus and spike history data from ./models
	%	model = a structure of fit coefficients from MLE_glmfit
	%	F = cursor dynamics from AR model
	%	Q = covariance matrix from AR model
	%	mu = average cursor data from AR model
	%	fn_out = base file name to write plots to
	%   
	%Output:
	%	output = Predicted decoded cursor from spike trains
	%  
	%Test code:
	%	%Load test preprocessed data
	%	fn_out = './worksheets/09_23_2014/testdecode.eps';
	%	%pre = load('./testdata/test_preprocess_spline_short24.mat');
	%	pre = load('./testdata/test_preprocess_spline.mat');
	%	const = 'on';
	%	nK_sp = 50; 
	%	nK_pos = 1;
	%	dt_sp = 0.002;
	%	dt_pos = 0.05;
	%	order = 1;
	%	data = filters_sp_pos(pre.processed, nK_sp, nK_pos, dt_sp, dt_pos);
	%	model = MLE_glmfit(data, const);
	%	[F, Q, mu] = fit_AR(data.torque, order);
	%	%load('./testdata/testglm.mat')
	%	decoded_torque = glm_decode(pre.processed, data, model, F, Q, mu, fn_out);

	binsize = processed.binsize;
	output = zeros(size(data.torque));
	nU = size(data.X,1);
	N = size(data.X,2);
	N_sub = 2000;
	%Only decode a short sample...
	N = 20000;


	%Make a Gaussian filter to smooth estimates
	sigma = 0.25;
	%sigma = 0.001;
	sigma = sigma/binsize;
	sz = sigma*3*3;
	x = linspace(-sz/2, sz/2, sz);
	gaussFilter = exp(-x.^2/(2*sigma^2));
	gaussFilter = gaussFilter/sum(gaussFilter);

	%Raise error if filter lengths are greater than 1... not implemented yet :(
	if length(data.k{2,2}) > 1 | length(data.k{3,2}) > 1
		error('NotImplementedError: cursor filters must currently be of length 1')
	end

	%GLM params
	%Cursor filters
	kx = model.b_hat(:,data.k{2,2}+1);
	ky = model.b_hat(:,data.k{3,2}+1);

	%Our initial estimate
	%xk = mu';
	%The actual point
	xk = data.torque(1,:)';
	Wk = Q;%eye(2);
	gradloglambda = [kx, ky];
	p = 0;

	%At each time step
	for idx = 1:N
		if floor(idx*10/N) ~= p
			display([num2str(floor(idx*100/N)) '% decoding completed'])
			p = p + 1;
		end
		f = squeeze(F(1,:,:));
		%Do prediction step 
		xk_p = mu'+f*xk;
		Wk_p = f*Wk*f'+Q;
		%Wk_p = f*Wk*f+Q;
		%Do update of covariance
		updateW = 0;
		updateX = 0;
		for j = 1:nU
			%j=9;
			%Evaluate GLM for each unit, at given time step
			%replace in data.X the actual cursor data for what's been predicted
			x = squeeze(data.X(j,idx,:));
			x(end-1:end) = xk_p;
			b_hat = model.b_hat(j,:);
			lambda = glmval(b_hat', x', 'log');
			%Sum update contributions for each unit
			%updateW = updateW + lambda*gradloglambda(j,:)'*gradloglambda(j,:)*binsize;
			%updateX = updateX + gradloglambda(j,:)'*(data.y(j,idx)-lambda*binsize);
			updateW = updateW + lambda*gradloglambda(j,:)'*gradloglambda(j,:);
			updateX = updateX + gradloglambda(j,:)'*(data.y(j,idx)-lambda);
		end
		%Don't let updateW get too big... or Wkinv will get too small...
		if updateW > 1e10
			display('Very large W')
		end
		updateW = min(updateW, 1e10);

		Wkinv = inv(Wk_p)+updateW;
		Wk = inv(Wkinv);
		%Update position
		xk = xk_p+Wk*updateX;
		%Save result
		output(idx,:) = xk;
		%display(['idx=' num2str(idx) ', time=' num2str(idx*binsize) ', lambda=' num2str(lambda)])
		%xk
		%Wk 
		%pause
	end

	%Plot ten seconds worth of data
	t_i = 0;
	%t_f = 40;
	t_f = 40; 
	ii = 1:N;
	tt = ii*binsize;
	ii = ii(tt > t_i & tt < t_f);
	tt = tt(ii);

	%Plot cursor data
	subplot(2,1,1);
	plot(tt, output(ii,1), 'r--', tt, data.torque(ii,1), 'r')
	xlim([t_i, t_f])
	ylabel('Pred. RU')
	%Take a random sample... so it doesn't take forever...
	ii_sub = datasample(ii, N_sub, 'Replace', false);
	corrRU = corr(output(ii_sub,1), data.torque(ii_sub,1))
	corrFE = corr(output(ii_sub,2), data.torque(ii_sub,2))
	title(['Correlation RU:' num2str(corrRU) ' Correlation FE: ' num2str(corrFE)])
	subplot(2,1,2);
	plot(tt, output(ii,2), 'r--', tt, data.torque(ii,2), 'r')
	xlim([t_i, t_f])
	ylabel('Pred. FE')
	xlabel('time (s)')		
	
	saveplot(gcf, fn_out, 'eps', [9 6]);
	saveas(gcf, [fn_out])
end