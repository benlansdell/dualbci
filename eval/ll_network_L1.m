function [ll, g] = ll_network_L1(prs, data, iU)
	nU = size(data.y,1);
	N = size(data.X,1);
	y = data.y(iU,:)';
	b_hat = prs;
	mu = glmval(b_hat, data.X, 'log');
	ll = sum(y.*log(mu))-sum(mu);
	g = grad(y, data.X, mu)';
end

function g = grad(y, X, mu)
	nP = size(X,2);
	g = zeros(1,nP+1);
	for j = 1:nP
		g(j+1) = sum(X(:,j).*(y - mu));
	end
	g(1) = sum(y - mu);
end 