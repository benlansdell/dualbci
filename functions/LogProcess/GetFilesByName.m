function fileList = GetFilesByName( name , workingDir)
%
% USAGE: GetFilesByName(name)
%
% This pulls a list of *unique* file base names in the current folder 
% matching a base string (e.g. animal name) with a .csv/.mat extension, 
% to facilitate running batch analyses on them.
%
% Note that it doesn't return the full filename, but everything before the
% extension, assuming you'll want to feed the resultant strings to LogRead
% which automatically appends .csv and .log

dirData = dir(workingDir);

fileList = {};
% yes, this loop is inefficient because it doesn't preallocate the array

for i = 1:length(dirData)
    
    % filenames must have form <name>_YYYY-MM-DD-HHMM.{mat,csv}
    match = regexp(dirData(i).name, strcat(name,'_[-0-9]{15}\.(mat|csv)'),'match');
    
    if (~isempty(match))
        fileList = [fileList; match{1}(1:(end-4))];  % append, stripping extension
    end

end

fileList = unique(fileList);

% now fileList is all the files you want, and none that you don't,
% hopefully