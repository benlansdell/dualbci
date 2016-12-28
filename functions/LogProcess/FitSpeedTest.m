% Test of execution speed for fit I've been using vs. doing the log mapping
% then a linear fit. And of course reproducing same results.

% Conclusion: shaves a few 0.01 seconds off PeakFit3, still less than 10%
% improvement. Not worth refactoring code to accomodate this.

% assumes metaData is loaded

model = metaData.Models{141};

logDW = log2(model.D+model.W/2);
logW = log2(model.W);

start = [1 1 1];
Nmin = 10;
dbFlag = false;

tic
[fitObj, gof, output, X, peaks, weights] = ...
                    PeakFit3([logDW,logW], model.MT, 'poly11', start, Nmin,false,1,dbFlag);
toc  
                
                % create Shannon-Welford fit model
                SWfit = fittype(@(a,b,k,d,w) min(100,max(-100,a + b .* log2( (d + w/2 ) ./ w.^k))),...
                    'independent',{'d','w'});
                
                SWopts = fitoptions(SWfit);
                set(SWopts,'Lower',[-100, -10, 0],'Upper',[100, 10, 10],'Robust','On');
                start = [1 1 1];    % starting parameter values
                
                tic
                [fitObj, gof, output, X, peaks, weights] = ...
                    PeakFit3([model.D,model.W], model.MT, SWfit, start, Nmin,false,1,dbFlag);
                toc