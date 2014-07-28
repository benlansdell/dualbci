close all
clear
clc

tic

load('generated_data_1.mat')

[x_len n_inputs ]=size(x);
[k_len k_co]=size(k);
if n_inputs~= k_co
    error('x and k must have the same number of rows (i.e. inputs)');
end
x=double(x); % data was in integer format, so this makes sure it is double as required by lassoglm
y=double(y);

x_mat=[]; % builds the toeplitz matrix with x.  
          % horizontally concatinates toeplitz matrices for each input channel
          % each row is a time step
for i_input=1:n_inputs
    x_mat=[x_mat toeplitz(x(:,i_input),[x(1,i_input) zeros(1,k_len-1)])];
end

x_mat_pad=[ones(x_len,1) x_mat]; % adds a column of ones for the constant term.  This is only used for evaluating the model later.

lever_hist_mat=toeplitz([0; y(1:(end-1))],[0 zeros(1,length(h)-1)]);

figure
subplot(1,5,1:4)
plot(x)
hold on
plot(rho,'r')
plot(y,'g')
hold off
subplot(1,5,5)
plot(k)
set(gcf,'Position',[74 535 1319 263])

% [bb,dev,stats] = glmfit([x_mat lever_hist_mat],y,'binomial','link','logit','constant','on'); % glmfit could be used instead.  glmfit outputs different stuff that lassoglm, so simply commenting this one in won't work right away
% [B,FitInfo] = lassoglm([x_mat lever_hist_mat],y,'binomial','Alpha',.25,'CV',10,'NumLambda',50,'LambdaRatio',10^-4); % if you have more time, more CV groups and more lambda
[B,FitInfo] = lassoglm([x_mat lever_hist_mat],y,'binomial','Alpha',.25,'CV',3,'NumLambda',20,'LambdaRatio',10^-4);

lassoPlot(B); % makes a plot of the coefficients over the L1 norm (a lasso plot)
lassoPlot(B,FitInfo,'plottype','CV') % plots deviance vs. lambda
%%

selected_lambda=FitInfo.IndexMinDeviance % index of best model based on min deviance
% selected_lambda=1 % alternatively just pick a model by hand

bb=B(:,selected_lambda); % bb is [k1; k2; h]
bbb=[FitInfo.Intercept(selected_lambda); bb]; % bbb is [b; k1, k2; h] This is the vector with all of fit coefficients

rho_hat = glmval(bbb,[x_mat_pad lever_hist_mat],'logit','constant','off'); % Uses the fit coefficients and the original input data to generate the ouput rho

figure % plot the original rho against the fit rho on the same input data
subplot(2,1,1)
plot(rho,'k')
subplot(2,1,2)
plot(rho_hat,'g')


figure % plot all of the fit coefficients against the orignal coefficents
k_fit=reshape(bb(1:(k_co*k_len)),k_len,k_co);
for i_input=1:n_inputs
    subplot(2,n_inputs,i_input)
    plot(k(:,i_input),'b')
    hold on
    plot(k_fit(:,i_input),'r')
end
hold off
subplot(2,n_inputs,(n_inputs+1):2*n_inputs)
plot(h)
hold on
h_fit=bb(k_co*k_len+1:end);
plot(h_fit,'r')
hold off
title(['generative b = ' num2str(b) ', fit b = ' num2str(FitInfo.Intercept(selected_lambda))])
legend('generative parameters','fit parameters')

toc

