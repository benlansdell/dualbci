
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Sensitivity to small changes in data%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
beta1 = 0.5; beta2 = 1;

%First, without correlations between x's
N = 100;
eps = 0.3;

x1 = randn(N,1);
x2 = randn(N,1);
X = [x1, x2];
y = beta1*x1+beta2*x2 + eps*randn(N,1);
%Fit for beta1_hat and beta2_hat
bhat1 = (X'*X)\X'*y;

%Plot in 3D
clf
subplot(1,2,1)
tt = -3:0.1:3;
[X, Y] = meshgrid(tt,tt);
Z = bhat1(1).*X+bhat1(2).*Y;
S = surf(X,Y,Z, 'facecol', [0.75, 0.75, 1], 'edgecol', 'no');%,...
alpha(0.5)
hold on
scatter3(x1, x2, y, '.b')

%Then repeat with another set of estimates
x1 = randn(N,1);
x2 = randn(N,1);
X = [x1, x2];
y = beta1*x1+beta2*x2 + eps*randn(N,1);
%Fit for beta1_hat and beta2_hat
bhat1 = (X'*X)\X'*y;

%Plot in 3D
[X, Y] = meshgrid(tt,tt);
Z = bhat1(1).*X+bhat1(2).*Y;
S = surf(X,Y,Z, 'facecol', [1, 0.75, 0.75], 'edgecol', 'no');%,...
alpha(0.5)
scatter3(x1, x2, y, '.r')
xlabel('x_1'); ylabel('x_2'); zlabel('y')

subplot(1,2,2)
%Plot 'energy function'
energy = @(b1, b2) (y-b1*x1-b2*x2)'*(y-b1*x1-b2*x2);
Z = zeros(size(X));
for i = 1:size(X,1)
	for j = 1:size(X,2)
		Z(i,j) = energy(tt(i), tt(j));
	end
end
pcolor(X, Y, Z);
hold on
xlabel('\beta_2')
ylabel('\beta_1')
plot(beta1, beta2, 'r.', 'linewidth', 20)
%saveplot(gcf, './worksheets/04_30_2015/correlatedfitsplot_indep.eps')
plot2svg('./worksheets/04_30_2015/correlatedfitsplot_indep.svg')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Then do the same for correlated X data...%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
beta1 = 0.5; beta2 = 1;

%First, without correlations between x's
N = 100;
eps = 0.3;
eps2 = 0.1;

x1 = randn(N,1);
x2 = x1+eps2*randn(N,1);
X = [x1, x2];
y = beta1*x1+beta2*x2 + eps*randn(N,1);
%Fit for beta1_hat and beta2_hat
bhat1 = (X'*X)\X'*y;

%Plot in 3D
clf
subplot(1,2,1)
tt = -3:0.1:3;
[X, Y] = meshgrid(tt,tt);
Z = bhat1(1).*X+bhat1(2).*Y;
S = surf(X,Y,Z, 'facecol', [0.75, 0.75, 1], 'edgecol', 'no');%,...
alpha(0.5)
hold on
scatter3(x1, x2, y, '.b')

%Then repeat with another set of estimates
x1 = randn(N,1);
x2 = x1+eps2*randn(N,1);
X = [x1, x2];
y = beta1*x1+beta2*x2 + eps*randn(N,1);
%Fit for beta1_hat and beta2_hat
bhat1 = (X'*X)\X'*y;

%Plot in 3D
[X, Y] = meshgrid(tt,tt);
Z = bhat1(1).*X+bhat1(2).*Y;
S = surf(X,Y,Z, 'facecol', [1, 0.75, 0.75], 'edgecol', 'no');%,...
alpha(0.5)
scatter3(x1, x2, y, '.r')
xlabel('x_1'); ylabel('x_2'); zlabel('y')

subplot(1,2,2)
%Plot 'energy function'
energy = @(b1, b2) (y-b1*x1-b2*x2)'*(y-b1*x1-b2*x2);
Z = zeros(size(X));
for i = 1:size(X,1)
	for j = 1:size(X,2)
		Z(i,j) = energy(tt(i), tt(j));
	end
end
pcolor(X, Y, Z);
hold on
xlabel('\beta_2')
ylabel('\beta_1')
plot(beta1, beta2, 'r.', 'linewidth', 20)
%saveplot(gcf, './worksheets/04_30_2015/correlatedfitsplot_corr.eps')
plot2svg('./worksheets/04_30_2015/correlatedfitsplot_corr.svg')
