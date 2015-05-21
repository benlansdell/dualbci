
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sensitivity to small changes in data for GLM%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
beta1 = 0.5; beta2 = 1;

%First, without correlations between x's
N = 100;
eps = 0.5;

x1 = randn(N,1);
x2 = randn(N,1);
X = [x1, x2];
mu = exp(beta1*x1+beta2*x2);
y = poissrnd(mu);
%Fit for beta1_hat and beta2_hat through IRLS
bhat1 = glmfit(X, y, 'poisson', 'const', 'off');

%Plot in 3D
clf
subplot(1,2,1)
tt = -3:0.1:3;
[X, Y] = meshgrid(tt,tt);
Z = exp(bhat1(1).*X+bhat1(2).*Y);
S = surf(X,Y,Z, 'facecol', [0.75, 0.75, 1], 'edgecol', 'no');%,...
alpha(0.5)
hold on
scatter3(x1, x2, y, '.b')

%Plot for another set
x1 = randn(N,1);
x2 = randn(N,1);
X = [x1, x2];
mu = exp(beta1*x1+beta2*x2);
y = poissrnd(mu);
%Fit for beta1_hat and beta2_hat through IRLS
bhat1 = glmfit(X, y, 'poisson', 'const', 'off')

%Plot in 3D
[X, Y] = meshgrid(tt,tt);
Z = exp(bhat1(1).*X+bhat1(2).*Y);
S = surf(X,Y,Z, 'facecol', [1, 0.75, 0.75], 'edgecol', 'no');%,...
alpha(0.5)
scatter3(x1, x2, y, '.r')

subplot(1,2,2)
%Plot 'energy function'
tt = -5:0.1:5;
[X, Y] = meshgrid(tt,tt);
energy = @(b1, b2) log(sum((y+0.25).*log((y+0.25)./exp(b1*x1+b2*x2))-((y+0.25)-exp(b1*x1+b2*x2))));
Z = zeros(size(X));
for i = 1:size(X,1)
	for j = 1:size(X,2)
		Z(i,j) = energy(tt(i), tt(j));
	end
end
pcolor(X, Y, Z);
hold on
xlabel('\beta_1')
ylabel('\beta_2')
plot(beta2, beta1, 'r.', 'linewidth', 20)
plot2svg('./worksheets/04_30_2015/correlatedfitsplot_GLM_indep.svg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Then do the same for correlated X data...%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
beta1 = 0.5; beta2 = 1;

%First, without correlations between x's
N = 100;
eps = 0.5;
eps2 = 0.2;

x1 = randn(N,1);
x2 = x1 + eps2*randn(N,1);
X = [x1, x2];
mu = exp(beta1*x1+beta2*x2);
y = poissrnd(mu);
%Fit for beta1_hat and beta2_hat through IRLS
bhat1 = glmfit(X, y, 'poisson', 'const', 'off')

%Plot in 3D
clf
subplot(1,2,1)
tt = -3:0.1:3;
[X, Y] = meshgrid(tt,tt);
Z = exp(bhat1(1).*X+bhat1(2).*Y);
S = surf(X,Y,Z, 'facecol', [0.75, 0.75, 1], 'edgecol', 'no');%,...
alpha(0.5)
hold on
scatter3(x1, x2, y, '.b')

x1 = randn(N,1);
x2 = x1 + eps2*randn(N,1);
X = [x1, x2];
mu = exp(beta1*x1+beta2*x2);
y = poissrnd(mu);
%Fit for beta1_hat and beta2_hat through IRLS
bhat1 = glmfit(X, y, 'poisson', 'const', 'off')

%Plot in 3D
[X, Y] = meshgrid(tt,tt);
Z = exp(bhat1(1).*X+bhat1(2).*Y);
S = surf(X,Y,Z, 'facecol', [1, 0.75, 0.75], 'edgecol', 'no');%,...
alpha(0.5)
scatter3(x1, x2, y, '.r')

subplot(1,2,2)
%Plot 'energy function'
tt = -5:0.1:5;
[X, Y] = meshgrid(tt,tt);
energy = @(b1, b2) log(sum((y+0.25).*log((y+0.25)./exp(b1*x1+b2*x2))-((y+0.25)-exp(b1*x1+b2*x2))));
Z = zeros(size(X));
for i = 1:size(X,1)
	for j = 1:size(X,2)
		Z(i,j) = energy(tt(i), tt(j));
	end
end
pcolor(X, Y, Z);
hold on
xlabel('\beta_1')
ylabel('\beta_2')
plot(beta2, beta1, 'r.', 'linewidth', 20)
plot2svg('./worksheets/04_30_2015/correlatedfitsplot_GLM_corr.svg')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Then do the same for intermediately correlated X data...%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
beta1 = 0.5; beta2 = 1;

%First, without correlations between x's
N = 100;
eps = 0.5;
eps2 = 0.5;

x1 = randn(N,1);
x2 = x1 + eps2*randn(N,1);
X = [x1, x2];
mu = exp(beta1*x1+beta2*x2);
y = poissrnd(mu);
%Fit for beta1_hat and beta2_hat through IRLS
bhat1 = glmfit(X, y, 'poisson')

%Plot in 3D
clf
subplot(1,2,1)
tt = -3:0.1:3;
[X, Y] = meshgrid(tt,tt);
Z = exp(bhat1(1).*X+bhat1(2).*Y);
S = surf(X,Y,Z, 'facecol', [0.75, 0.75, 1], 'edgecol', 'no');%,...
alpha(0.5)
hold on
scatter3(x1, x2, y, '.b')

x1 = randn(N,1);
x2 = x1 + eps2*randn(N,1);
X = [x1, x2];
mu = exp(beta1*x1+beta2*x2);
y = poissrnd(mu);
%Fit for beta1_hat and beta2_hat through IRLS
bhat1 = glmfit(X, y, 'poisson')

%Plot in 3D
[X, Y] = meshgrid(tt,tt);
Z = exp(bhat1(1).*X+bhat1(2).*Y);
S = surf(X,Y,Z, 'facecol', [1, 0.75, 0.75], 'edgecol', 'no');%,...
alpha(0.5)
scatter3(x1, x2, y, '.r')

subplot(1,2,2)
%Plot 'energy function'
tt = -5:0.1:5;
[X, Y] = meshgrid(tt,tt);
energy = @(b1, b2) log(sum((y+0.25).*log((y+0.25)./exp(b1*x1+b2*x2))-((y+0.25)-exp(b1*x1+b2*x2))));
Z = zeros(size(X));
for i = 1:size(X,1)
	for j = 1:size(X,2)
		Z(i,j) = energy(tt(i), tt(j));
	end
end
pcolor(X, Y, Z);
hold on
xlabel('\beta_1')
ylabel('\beta_2')
plot(beta2, beta1, 'r.', 'linewidth', 20)
plot2svg('./worksheets/04_30_2015/correlatedfitsplot_GLM_intermedcorr.svg')