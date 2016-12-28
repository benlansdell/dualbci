function UpdateMatFiles(fNameBase,logpath, matpath, dbFlag)
% USAGE: UpdateMatFiles(name, logpath, matpath, dbFlag)
%
% This function checks all .csv files matching fNameBase in the current
% directory, and runs LogRead() on them. LogRead() checks for existing .mat
% and modify dates to determine whether it needs to parse the data or not.
% In this way .mat files for every logfile can automatically be
% kept up to date as new experimental data is recorded, or the log parsing
% code is updated. You're welcome.
%
% INPUTS:   fNameBase       e.g. name of animal
%           updateSummary   (flag) if true, re-generate summary text file
%                           of all manually-entered log messages for this
%                           base name
%

% So far, just one monkey and no rats.
if (nargin == 1)
    logpath = uigetdir(pwd,'Select Log File (.csv/.log) Directory');
    matpath = uigetdir(pwd,'Select Output (.mat) Directory');
end

if (nargin < 4)
    dbFlag = false;
end

fileList = GetFilesByName(fNameBase,logpath);
%fileList = fileList(end:-1:1);      % Update most recent files first

for idx = 1:length(fileList)
        fname = fileList{idx};

        fprintf('Processing %s ...\n',fname )
        % parse the .log and .csv files into a matlab data structure, which
        % is written to a .mat file by this script (and left in workspace)
        LogRead(fname,logpath,matpath,dbFlag);  
        
        clear data ans
     
end

end



    
    

    