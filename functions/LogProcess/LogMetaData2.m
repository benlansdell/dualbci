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
%
% VERSION 2:    Do log transform of trial times before doing a straight-up
%               Fitts's law or Shannon-Welford analysis


function LogMetaData2(name, dbFlag, nofit)
global matpath
try 
if (nargin < 2)
    dbFlag = 0;
end

if (nargin < 3)
    nofit = false;
end

if (nargin == 0)
    name = 'Spanky'
    dbFlag = 1
    clc
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
        Maps = setdiff(unique(trials.Map),{''});
        nMaps = length(Maps);
        nConds = nSegs + nMaps;

        Day = metaData.fileList{fi};

        
        for j = 1:nMaps % only look at all trials w/ same map (change to nConds to check every segment)
            
            if (j > nMaps)      % per-segment (trial block) analysis
                try
                Map = ...
                    trials.Map{find(trials.segment == j-nMaps,1,'first')};
                if (length(unique(trials.segment(ismember(trials.Map, Map)))) > 1)

                trialIdx = find(trials.segment == j-nMaps & trials.valid & trials.success == 1);
                seg = j-nMaps;    % indicates which part of expt this data is from
              
                else
                    fprintf('Skipping redundant segment for %s\n',Map);
                    continue
                end
                catch me
                    keyboard
                end
                
            else                % per-mapping condition analysis    
                Map = Maps{j};
                trialIdx = find(ismember(trials.Map, Map)' & trials.valid & trials.success == 1)';
                seg = 0;    % indicates that this is a condition-based data set                
                %segIdx = unique(trials.segment(trialIdx));   % which segments in this one
                
            end
            fprintf('%d\t%25s\t %d trials\t',seg,Map,length(trialIdx));
            
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
            
            % kill stupid values (negative difficulty, 0 movement time ?!)
            tkeep = (MT > 0) & (D-W/2)./W + 1 > 0;
            D = D(tkeep);
            W = W(tkeep);
            MT = MT(tkeep);
            
            % generate enumeration of unique DxW combinations
            DxW = unique([D W], 'rows');
            DxWi = 0*D;
            for k = 1:length(D)
                DxWi(k) = find( DxW(:,1) == D(k) & DxW(:,2) == W(k));
            end
 
            clear Models;
            
            %Models.D = D;  % can be reconstructed from [D W] = DxW(DxWi,:)
            %Models.W = W;
            Models.DxW = DxW;
            Models.DxWi = DxWi;
            Models.logMT = log(MT);
            Models.idxSet = trialIdx;    % indices of this segment's trials
            Models.trialsIdx = i;        % index of this day
            Models.IDcount = length(unique(D./W));
            Models.TrialCount = length(trialIdx);
            Models.seg = seg;
            Models.Map = Map;
            Models.Day = Day;
            
            try
            
            if (Models.DxW > Models.IDcount)
                %  keyboard
            end
            X = [D, W];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Per-dataset curve fitting and stats
            if (~nofit && length(unique(D)) > 1 && size(DxW,1) > 3)

                sigCount = sigCount + 1;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Null model fit (for 'absolute' F-tests)
                
                [B,BINT,R] = regress(Models.logMT, ones(size(Models.logMT)));
                Models.NULL.B = B;
                Models.NULL.BINT = BINT;
                Models.NULL.R = R;
                SSnull = sum(R.^2);
                DFnull = length(Models.logMT) - 1;
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Shannon-Welford fit model
                SWfit = fittype(@(a,b,k,d,w) log(max(eps,a + b .* log2( (d + w/2 ) ./ w.^k))),...
                    'independent',{'d','w'});
                
                SWopts = fitoptions(SWfit);
                set(SWopts,'Lower',[-100, -10, 0],'Upper',[100, 10, 10],'Robust','On');
                start = [1 1 1];    % starting parameter values
                
                clear fo gof out
                [fo, gof, out] = fit(X,log(MT),SWfit,'StartPoint',start);
               
                Models.SW.fo = fo;
                Models.SW.gof = gof;
                Models.SW.res = out.residuals;
                
                % F-testing
                SS_SW = sum(out.residuals.^2);
                DF_SW = length(Models.logMT) - 3;
                Models.SW.F = (SSnull - SS_SW)/(DFnull - DF_SW) /...
                                (SS_SW / DF_SW);
                Models.SW.P = 1-fcdf(Models.SW.F, DFnull, DF_SW);
                
                
                %%%% Optionally plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if (dbFlag || gof.rsquare > 0.5)
                    figure(1),subplot(2,1,2)
                    plot(D,log(MT),'k.')
                    if (gof.rsquare > 0.1)
                        plot(fo,X,log(MT),'Style','Residuals','YLim',[min(log(MT)) max(log(MT))])
                        
                    end
                end
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Fitts's law model (Shannon form)
                FittsFit = fittype(@(a,b,d,w) log(a + b*log2((d-w/2)./w + 1)),...
                    'independent',{'d','w'});
                FittsOpts = fitoptions(FittsFit);
                set(FittsOpts,'Lower',[0, 0],'Upper',[5, 5],'Robust','On');
                start = [0 1]; % starting parameter values

                clear fo gof out
                [fo, gof, out] = fit(X,log(MT),FittsFit,'StartPoint',start);
              
                Models.Fitts.fo = fo;
                Models.Fitts.gof = gof;
                Models.Fitts.res = out.residuals;
                
                % F-testing
                SS_Fitts = sum(out.residuals.^2);
                DF_Fitts = length(Models.logMT) - 2;
                Models.Fitts.F = (SSnull - SS_Fitts) / (DFnull - DF_Fitts) / ...
                                    (SS_Fitts / DF_Fitts);
                Models.Fitts.P = 1-fcdf(Models.Fitts.F, DFnull, DF_Fitts);
               
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % F-test of Shannon-Welford vs. Fitts's law
                Models.Fversus = (SS_Fitts - SS_SW) / (DF_Fitts - DF_SW) /...
                                    (SS_SW / DF_SW);
                Models.Pversus = 1-fcdf(Models.Fversus, DF_Fitts, DF_SW);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if (~isempty(Models.Fitts.gof) && ...
                        ~isempty(Models.SW.gof))
                    Models.valid = true;
                else
                    Models.valid = false;
                    fprintf('fitting failure!\n')
                end
                    fits = fits + 1;
                    metaData.Models{fits} = Models;
                
                    fprintf('...saved model fit data.\n');
            else
                fprintf('\n');
            end
            
            catch me
                Models.valid = false;
                Models.error = me.getReport;
                fits = fits + 1;
                metaData.Models{fits} = Models;
                if (dbFlag)
                    keyboard
                end
            end
            
        end
    end
    metaData.trials{fi} = trials;
    
    % every time we save, the *entire* mat file gets overwritten, gets very
    % slow near the end if there are many files.
    if dbFlag
        save(strcat(matpath,'/',name,'_LogMetaData-test2.mat'),'metaData')
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
        save(strcat(matpath,'/',name,'_LogMetaData-test2.mat'),'metaData')
    else
        save(strcat(matpath,'/',name,'_LogMetaData2.mat'),'metaData')
    end
    
catch me
    me.getReport
    keyboard
end
