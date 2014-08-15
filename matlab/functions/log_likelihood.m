function l = log_likelihood(rho_hat, y)
	L = exp(-rho_hat).*(rho_hat).^y./factorial(y);
	l = sum(log(L));
end