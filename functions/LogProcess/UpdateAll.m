% UpdateAll.m
% by C. Matlack 2013-10-15
%
% USAGE: 
% 
% UpdateAll('AnimalName','LogPath','OutPath')
%
% UpdateAll('AnimalName')  will pop up selection dialogs for the paths
%
% Parses .csv and .log files created by labview, creates .mat file for each
% valid pair of .csv and .log files, and creates metadata .mat and .csv
% files. Only creates .mat files as needed, based on modify date of
% .csv/.log compared to .mat file (and existence of .mat file).
% 
% ARGUMENTS:    'AnimalName'    e.g. one of {'Lamar','Optimus','Spanky'}
%               'LogPath'       full path to folder with log
%               'OutPath'       full path for output .mat etc files


function UpdateAll(name, dbFlag)
global logpath matpath

if (nargin == 0)
    sprintf('Select Log File (.csv/.log) Directory\n\n');
    logpath = uigetdir(pwd,'Select Log File (.csv/.log) Directory')

    sprintf('Select Output (.mat) Directory\n\n');
    matpath = uigetdir(pwd,'Select Output (.mat) Directory')
end

if (nargin < 2)
    dbFlag = false;
end

% Quick sanity check before we get to business.
assert(isdir(logpath) && isdir(matpath),'Valid .csv/.log/.mat directories not provided, aborting.\n')

% Check to see if we have any mat files missing their log files
[M,~,~] = VerifyLogs(name,logpath,matpath,dbFlag);
if (~isempty(M))
    fprintf('WARNING: Found .mat files with no corresponding logfiles.\n')
    fprintf('%s\n',M{:})
    fprintf('Press any key to continue.\n')
    pause
end

% todo: popup message telling user how many days' logfiles will be
% processed

% First update all mat files
UpdateMatFiles(name,logpath, matpath, dbFlag)

fprintf('\n========== Done processing daily logs, now generating meta data ===========\n')
% Then update LogMetaData file. Call with the nofit flag, so it just pulls
% the trial data and segment (trial block) data, latter needed for NEV sync
LogMetaData(name,false,true)

% will be needed for NEV processing
load([matpath filesep name '_LogMetaData.mat'])

% Then update summary files
LogMetaData_Summary(name,matpath)

UpdateLogSummary(name,matpath)

% Finally update mapping from NEV file times and channels to Labview timing
% and channels. Requires metaData segment information; takes several hours
NEV_UpdateAllMAT(name)

end
    