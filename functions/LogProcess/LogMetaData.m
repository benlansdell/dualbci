% LogMetaData.m
% by C. Matlack 2013-06-12
%
% USAGE:    LogMetaData(name, matpath, <dbFlag>)
%
% This script compiles trial-by-trial data for *all* days of data matching
% an animal name.
%
% DEPENDS:      ARTANOVA
%               GetFilesByName
%               LogMetaParser

function LogMetaData(name, dbFlag, nofit)
global matpath

if (nargin < 2)
    dbFlag = 0;
end

if (nargin < 3)
    nofit = false;
end

global f1 f2 a1 a2  % figure handles in PeakFit

%%%%%%%%% IMPORTANT SUBJECTIVE PARAMETERS %%%%%%%%%%%%

Nmin = 10;      % Minimum number of trials at a DxW condition for peakfit

tol = 0.01;     % rounding tolerance for Indicies of Difficulty and Time (sec)

testDxW = [0.4 0.3];

%
% Cell array columns:
%
% Hand-entered:
%  {Neurons = [A.a, B.b, C.c, D.d]}
%  {Important = {-1: ignore; 0: default; 1: exceptional}
%
% {Animal = <Lamar, Spanky>}
% {Date = <YYYY-MM-DD-HHMM>}

%% Get all filenames
fprintf('\nLooking for .mat files in:\n%s\n',matpath);
fileList = GetFilesByName(name,matpath)';

special = find(ismember(fileList,'Spanky_2013-02-08-1251'));

if (dbFlag)
    nfiles = length(fileList);
    fileList = fileList([special max(end-5,1):end]);
    fprintf('Cutting file list from %d to most recent %d for testing...\n',...
        nfiles,length(fileList));
end

metaData.name = name;
metaData.fileList = fileList;
metaData.numFiles = length(metaData.fileList);
fprintf('Processing %d data files...\n',metaData.numFiles);

% create index mapping from existing database entries to files
%idxMap = zeros(metaData.numFiles,1);
%for i = 1:metaData.numFiles
%    if 

%% Extract trial data from all

metaData.trials = cell(metaData.numFiles,1);

fits = 0;   % number of segments fittable w/ Fitts' Law
sigCount = 0;   % count of number of segments w/ significant effect of width/distance
keep = false(metaData.numFiles,1);
i = 0;

for fi = 1:metaData.numFiles
    fprintf('\nFile #%d, %s ...\n',fi, metaData.fileList{fi});
    
    
    load(strcat(matpath,filesep,metaData.fileList{fi},'.mat'));  % creates 'data' object
    trials = LogMetaParser(data, data.messages.time(1),data.messages.time(end),dbFlag);
    
    if (~isempty(trials))
        i = i + 1;
        keep(fi) = true;
        
        %metaData.FittsAnalysis{i} = [];
        
        
        % Check sets of trials by two criteria:
        % 1) continuous segments w/ same mapping
        % subtract one because orphan trials are labeled with 0
        nSegs = length(unique(trials.segment)) - 1;
        % 2) just same mapping
        if (nSegs > 1)
            nMaps = nSegs + length(setdiff(unique(trials.Map),{''}));
        else
            nMaps = nSegs;
        end
        FittsAnalysis.Maps = cell(nMaps,1);
        Maps = setdiff(unique(trials.Map),{''});
        
        
        %nMaps = length(FittsAnalysis.Maps);
        Day = metaData.fileList{fi};
        %         FittsAnalysis.A = zeros(nMaps,1);
        %         FittsAnalysis.Fscores = cell(nMaps,1);
        %         FittsAnalysis.Rsquared = zeros(nMaps,1);
        %         FittsAnalysis.Model = cell(nMaps,1);
        %         FittsAnalysis.Fits = {};
        
        for j = 1:nMaps
            if (j <= nSegs)
                %FittsAnalysis.Maps{j} = ...
                try
                Map = ...
                    trials.Map{find(trials.segment == j,1,'first')};
                trialIdx = find(trials.segment == j & trials.valid & trials.success == 1);
                seg = j;    % indicates which part of expt this data is from
                catch me
                    keyboard
                end
                %fprintf('(skipped)\n')
                %continue % skip these for now, just look at condition-based sets
                
            else    
                %FittsAnalysis.Maps{j} =
                Map = Maps{j-nSegs};
                
                if (length(unique(trials.segment(ismember(trials.Map, Map)))) > 1)
                    
                trialIdx = find(ismember(trials.Map, Map)' & trials.valid & trials.success == 1)';
                seg = 0;    % indicates that this is a condition-based data set
                
                segIdx = unique(trials.segment(trialIdx));   % which segments in this one
                
                else
                    fprintf('Skipping redundant mapping group for %s\n',Map);
                    continue
                end
            end
            fprintf('%d\t%25s\t %d trials\n',seg,Map,length(trialIdx));
            
            if (isempty(trialIdx))
                continue
            end
            MT = trials.MT(trialIdx);
            D = trials.D(trialIdx);
            W = trials.W(trialIdx);
            
            % round off trial parameters for good repeated-measures set
            D = round(D/tol)*tol;
            W = round(W/tol)*tol;
            MT = round(MT/tol)*tol;
            
            % generate enumeration of unique DxW combinations
            DxW = unique([D W], 'rows');
            DxWi = 0*D;
            for k = 1:length(D)
                DxWi(k) = find( DxW(:,1) == D(k) & DxW(:,2) == W(k));
            end
            
            if (~nofit)
                
                % Test #1: ART + ANOVA test of significant effects for either
                % distance and width combination, or their ratio
                
                [F,P] = ARTANOVA( [D W], MT, {'D','W','MT'});
                [f,p] = ARTANOVA( log2((D+W/2)./W + 1), MT,{'ID','MT'});                
                F = [F f];
                P = [P p];
                
               
                
                clear Models;
                
                Models.D = D;
                Models.W = W;
                Models.DxW = DxW;
                Models.DxWi = DxWi;
                Models.MT = MT;
                Models.idxSet = trialIdx;    % indices of this segment's trials
                Models.trialsIdx = i;        % index of this day
                Models.P = P;
                Models.F = F;
                Models.IDcount = length(unique(D./W));
                Models.TrialCount = length(trialIdx);
                Models.seg = seg;
                Models.Map = Map;
                %Fitts.Fscores = FP;
                %Fitts.A = A;
                Models.Day = Day;
                
                if (Models.DxW > Models.IDcount)
                    %  keyboard
                end
                
                sigCount = sigCount + 1;
                
                % create Shannon-Welford fit model
                SWfit = fittype(@(a,b,k,d,w) min(100,max(-100,a + b .* log2( (d + w/2 ) ./ w.^k))),...
                    'independent',{'d','w'});
                
                SWopts = fitoptions(SWfit);
                set(SWopts,'Lower',[-100, -10, 0],'Upper',[100, 10, 10],'Robust','On');
                start = [1 1 1];    % starting parameter values
                
                % Fit the Shannon-Welford model
                %if (length(unique(W)) > 1)
                
                [X, peaks, weights] = PeakFind([D, W],MT,Nmin,false,false);
                
                if (~isempty(peaks))
                    clear gof
                    
                    [fitObj, gof, output] = fit(X,peaks,SWfit,'Weight',weights,'StartPoint',[1 1 1]);
                    
                    %                 [fitObj, gof, output, X, peaks, weights] = ...
                    %                     PeakFit3([D,W], MT, SWfit, start, Nmin,false,1,dbFlag);
                    
                    
                    if ( ~isempty(gof)) % we have a significant data set
                        %fprintf('Got adj R^2 of %1.2f, recording...\n',gof.adjrsquare);
                        
                        Models.SWfitObj = fitObj;
                        Models.SWgof = gof;
                        Models.SWoutput = output;
                        Models.SWk = fitObj.k;
                        Models.SWrsquare = gof.rsquare;
                        Models.SWX = X;
                        Models.SWpeaks = peaks;
                        Models.SW = true;
                        %                     if (gof.rsquare > 0.50)       DO THIS SEPARATELY
                        %                         paramSet = BootstrapAnalysis(Models,100,false);
                        %                         clear bootstrap
                        %                         try
                        %                         bootstrap.paramSet = mean(paramSet);
                        %                         bootstrap.stdParams = std(paramSet);
                        %                         bootstrap.count = sum(~isnan(paramSet(:,1)));
                        %                         Models.bootstrap = bootstrap;
                        %                         Models.SWk = paramSet(3);
                        %                         Models.SWrsquare = paramSet(4);
                        %                         catch me
                        %                             keyboard
                        %                         end
                        %                     end
                    else
                        Models.SW = false;
                    end
                else
                    Models.SW = false;
                end
                
                % define Fitts Law model
                % FittsFit = fittype( @(a,b,x) a + b*x);
                % ID = log2((D-W/2)./W + 1);
                
                % Fitts's law model for 2 indep variables
                FittsFit = fittype(@(a,b,d,w) a + b*log2((d-w/2)./w + 1),...
                    'independent',{'d','w'});
                start = [1 1]; % starting parameter values
                
                if (~isempty(peaks))
                    clear gof
                    [fitObj, gof, output] = fit(X,peaks,FittsFit,'Weight',weights,'StartPoint',[1 1]);
                    
                    %                 [fitObj, gof, output, X, peaks, weights] = ...
                    %                     PeakFit3([D,W],...
                    %                     MT, FittsFit, start, Nmin,false,1,dbFlag);
                    %gof = [];
                    
                    if (~isempty(gof))
                        Models.fitObj = fitObj;
                        Models.gof = gof;
                        Models.output = output;
                        Models.X = X;
                        Models.peaks = peaks;
                        Models.weights = weights;
                        Models.Fitts = true;
                    else
                        Models.Fitts = false;
                    end
                else
                    Models.Fitts = false;
                end
                
                % Do a targets-per-minute analysis for comparison
                    [test,mi] = ismember(testDxW,DxW,'rows');
                    if (test)
                        testMTs = MT(DxWi == mi);
                        Models.TPMTrials = testMTs;
                        Models.TargetsPerMin = 60.0 / mean(testMTs);
                    else
                        Models.TPMTrials = [];
                        Models.TargetsPerMin = NaN;
                    end
                
                if (Models.Fitts || Models.SW)
                    fits = fits + 1;
                    try
                    metaData.Models{fits} = Models;
                    catch me 
                        keyboard
                    end
                    fprintf('...saved model fit data.\n');
                end
                %end
                
                %[A,FP] = FittsModelTest(MT,D,W,Nmin);
               
            end
        end
    end
    metaData.trials{fi} = trials;
    metaData.sigCount = sigCount;       % number of days with significant trend
    
    % every time we save, the *entire* mat file gets overwritten, gets very
    % slow near the end if there are many files.
    if dbFlag
        save(strcat(matpath,'/',name,'_LogMetaData-test.mat'),'metaData')
    else
        %save(strcat(matpath,'/',name,'_LogMetaData.mat'),'metaData')
    end
end
try
    
metaData.trials = metaData.trials(keep);  
metaData.fileList = metaData.fileList(keep);

catch me
    keyboard
end
% by days with insufficient data
    if dbFlag
        save(strcat(matpath,'/',name,'_LogMetaData-test.mat'),'metaData')
    else
        save(strcat(matpath,'/',name,'_LogMetaData.mat'),'metaData')
    end
