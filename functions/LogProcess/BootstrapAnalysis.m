% 
% Bootstrapping analysis function for use on a single model instance
% (tests Fitts and Shannon-Welford)

function bootstrap = BootstrapAnalysis(model,Nreps,Nmin, plotFlag)

if (nargin < 4)
    plotFlag = true;
end

failCount = 0;

bootstrap = [];
paramSet = zeros(Nreps,6);      % each row is [a b k R^2 AdjR^2]

% create Shannon-Welford fit model
SWfit = fittype(@(a,b,k,d,w) min(100,max(-100,a + b .* log2( (d + w/2 ) ./ w.^k))),...
    'independent',{'d','w'});
% Fitts's law model for 2 indep variables
FittsFit = fittype(@(a,b,d,w) a + b*log2((d-w/2)./w + 1),...
    'independent',{'d','w'});
SWopts = fitoptions(SWfit);
set(SWopts,'Lower',[-100, -10, 0],'Upper',[100, 10, 10],'Robust','On');
start = [1 1 1];    % starting parameter values

for i = 1:Nreps
    
    % Reselect with replacement from the trial set
    Ntrials = model.TrialCount;
    SelectIdx = ceil(rand(1,model.TrialCount)*(model.TrialCount-1));
    
    W = model.W(SelectIdx);
    D = model.D(SelectIdx);
    MT = model.MT(SelectIdx);
    
    % See what our selection of trials looks like against the model
    
    % Find the representative trial times from the data set
    [X, peaks, weights] = PeakFind([D, W],MT,Nmin,false,false);
    
    
    % %%%%% OLD SLOW CODE TO GET RID OF %%%%%%%%%%%%
    %     [fitObj, gof2, output2, X2, peaks2, weights2] = ...
    %      PeakFit3([D,W],...
    %      MT, FittsFit, [1 1], Nmin,false,true,false);
    %   clear gof1
    %
    %  [fitObj, gof1, output1, X1, peaks, weights1] = ...
    %      PeakFit3([D,W], MT, SWfit, start, Nmin,false,true,false);
    
    %%%%%% NEW SHIT THAT MIGHT BE WRONG %%%%%%%%%%%%%
    if (~isempty(peaks))
        [fitObj2, gof2, output2] = fit(X,peaks,FittsFit,'Weight',weights,'StartPoint',[1 1]);
        [fitObj, gof1, output1] = fit(X,peaks,SWfit,'Weight',weights,'StartPoint',[1 1 1]);
        
        % run the fit by pre-processing input terms and using poly11
        %%% Don't know why this stuff is not worknig but it is not!
        % Shannon-Welford
        %     [fo, gofp, output1] = fit([log2(X(:,1)+X(:,2)/2) log2(X(:,2))],peaks,'poly11','Weight',weights);
        %     k = -fo.p01/fo.p10;
        %     b = fo.p10;
        %     a = fo.p00;
        %     % Fitts
        %     [fo, gofp, output2] = fit([log2( (X(:,1)./X(:,2) + 0.5) )], peaks, 'poly1','Weight',weights);
        %
        
        if (isempty(gof1) || isempty(gof2))
            paramSet(i,:) = nan(1,6);
            failCount = failCount + 1;
            if (failCount > Nreps/2)
                fprintf('Too many data failures, bailing on bootstrap...\n');
                %paramSet = nan(1,5);
                bootstrap = [];
                return
            end
        else
            % sanity check
            
            %assert(X1 == X2 && weights2 == weights1,'Data points used were different!')
            
            % F-test sanity check
            sseFitts = sum((output2.residuals.^2));% .* weights2);
            sseSW = sum((output1.residuals.^2)); %.* weights1);
            
            %      if(sseSW > sseFitts)
            %          fprintf('F-test abysmal failure\n');
            %      end
            %
            N = size(unique(model.SWX,'rows'),1);   % number of points fit
            Fcrit = finv(0.05,1,N-3);                           % critical value of F-stat
            F = (sseFitts - sseSW)/(3-2)/sseSW/(N-3);   % value of F-stat
            paramSet(i,:) = [fitObj.a fitObj.b fitObj.k gof1.rsquare gof2.rsquare F > Fcrit];
            
            %disp([fitObj.a fitObj.b fitObj.k])
        end
    end
    
end

%%
paramSet = paramSet(~isnan(paramSet));
paramSet = reshape(paramSet,numel(paramSet)/6,6);
if (isempty(paramSet))
    return;
end
fprintf('\ta\tb\tk\tR^2\tFittsR^2\tF-test\n');
fprintf('Means:\n')
disp(mean(paramSet))
fprintf('Std Dev:\n')
disp(std(paramSet))
if (plotFlag)
    figure(10), clf
    
    labels = {'a','b','k','R^2','FR^2','F-test'};
    for s = 1:size(paramSet,2)
        subplot(size(paramSet,2),1,s)
        hist(paramSet(:,s),round(Nreps/3))
        title(labels{s})
    end
    
%     title('Fit Parameters')
%     subplot(1,2,2),boxplot(paramSet(:,4:5),'labels',{'R^2','Adj R^2'})
%     title('Goodness-of-fit')
end

bootstrap.k = mean(paramSet(:,3));
bootstrap.kCI = prctile(paramSet(:,3),[2.5,97.5]);
bootstrap.a = mean(paramSet(:,1));
bootstrap.b = mean(paramSet(:,2));
bootstrap.kstd = std(paramSet(:,3));
bootstrap.rsquare = mean(paramSet(:,4));
bootstrap.rsquarestd = std(paramSet(:,4));
bootstrap.results = paramSet;
bootstrap.Ftests = mean(paramSet(:,6));




bootstrap.rsquareFitts = mean(paramSet(:,5));
bootstrap.rsquarestdFitts = std(paramSet(:,5));

end
