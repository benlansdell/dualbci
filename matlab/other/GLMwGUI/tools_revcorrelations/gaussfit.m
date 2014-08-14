function [sigma,mu] = gaussfit(x,y)

% Maximum number of iterations
Nmax = 50;

n = length(x);
x = reshape(x,n,1);
y = reshape(y,n,1);

%sort according to x
X = [x,y];
X = sortrows(X);
x = X(:,1);
y = X(:,2);

%Checking if the data is normalized
dx = diff(x);
dy = 0.5*(y(1:length(y)-1) + y(2:length(y)));
s = sum( dx .* dy );
if( s > 1.5 || s < 0.5 )
    y = y./s;
end

X = zeros(n,3);
X(:,1) = 1;
X(:,2) = x;
X(:,3) = (x.*x);

% try to estimate mean mu from the location of the maximum
[ymax,index]=max(y);
mu = x(index);

% estimate sigma
sigma = 1/(sqrt(2*pi)*ymax);

for i=1:Nmax

    dfdsigma = -1/(sqrt(2*pi)*sigma^2)*exp(-((x-mu).^2)/(2*sigma^2));
    dfdsigma = dfdsigma+1/(sqrt(2*pi)*sigma).*exp(-((x-mu).^2)/...
        (2*sigma^2)).*((x-mu).^2/sigma^3);

    dfdmu = 1/(sqrt(2*pi)*sigma)*exp(-((x-mu).^2)/(2*sigma^2)).*(x-mu)/...
        (sigma^2);

    F = [dfdsigma dfdmu];
    a0 = [sigma;mu];
    f0 = 1/(sqrt(2*pi)*sigma).*exp(-(x-mu).^2/(2*sigma^2));
    a = (F'*F)^(-1)*F'*(y-f0) + a0;
    sigma = a(1);
    mu = a(2);
end