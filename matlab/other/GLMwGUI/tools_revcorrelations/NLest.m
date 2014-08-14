function Thresh = NLest(FR,PriorDist,FiltteredDist)
% This function will calculate the estimation of nonlinear relation between
% input and output data. The nonlinearity can be estimated by using Bayes'
% rule P(spike|s) = P(spike)P(s|spike)/P(s). FR := mean firing rate
%
% 31st July 2013
% Author : SangWook Lee
len = length(PriorDist);
    
thresh = (FR.*FiltteredDist)./PriorDist;

% Cutting off the last 90 percent of whole to prevent high fluctuation of
% the last part(where the prior and spike distribution are near zero).
thresh = thresh(1:round(len*0.8));
thresh = smoothing(thresh,9);

Thresh = interp1(1:length(thresh),thresh,...
    linspace(1,length(thresh),len));
Thresh = Thresh(:);