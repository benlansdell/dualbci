% VerifyLogs.m
%
% USAGE:    [M,L,C] = VerifyLogs(name, logpath, matpath, debug)
%
%                     VerifyLogs(name, logpath, matpath)    prints output
%
%                     VerifyLogs(name)              will prompt for paths
%
%                     VerifyLogs()     will check *all* animals' logfiles
%
% This script will ask for you to specify two folders. It will search one
% for .mat files and the other for .csv/.log files, and report back on
% where any inconsistencies exist between the two. Useful for confirming
% .csv/.log backups of existing .mats, and whether metada files should be
% updated.
%
% M = list of orphan MAT files
% L = ...           .log files
% C = ...           .csv files

function [M,L,C] = VerifyLogs(name, logpath, matpath, dbFlag)

if (nargin < 4)
    dbFlag = 1;     % print output to console
end

if (nargin < 3)
    sprintf('Select Log File (.csv/.log) Directory\n\n');
    logpath = uigetdir(pwd,'Select Log File (.csv/.log) Directory')

    sprintf('Select Output (.mat) Directory\n\n');
    matpath = uigetdir(pwd,'Select Output (.mat) Directory')
    
    if (nargin == 0)    % this will run the check on all animals at once
        name = '\w+';
    end
end

% Quick sanity check before we get to business.
assert(isdir(logpath) && isdir(matpath),'Valid .csv/.log/.mat directories not provided, aborting.\n')

% get directory listing of MAT and log file directories
matdir = dir(matpath);
logdir = dir(logpath);

% search for files matching time/date format in name + correct extension
mats = unique(regexp({matdir(:).name},strcat(name,'_[-0-9]{15}.mat'),'match','once'));
mats = regexp(mats,'[^.]+','match','once'); % strip extension

csvs = unique(regexp({logdir(:).name},strcat(name,'_[-0-9]{15}.csv'),'match','once'));
csvs = regexp(csvs,'[^.]+','match','once'); % strip extension
logs = unique(regexp({logdir(:).name},strcat(name,'_[-0-9]{15}.log'),'match','once'));
logs = regexp(logs,'[^.]+','match','once'); % strip extension


% compare lists to see if any orphans exist
M = union(setdiff(mats,csvs),setdiff(mats,logs));
L = setdiff(logs,mats);
C = setdiff(csvs,mats);

if (dbFlag)
fprintf('=====MAT files without both corresponding logfiles:\n');
fprintf('%s\n',M{:})
fprintf('\n=====Log files (.log) without corresponding mat file:\n');
fprintf('%s\n',L{:})
fprintf('\n=====Log files (.csv) without corresponding mat file:\n');
fprintf('%s\n',C{:})
end

