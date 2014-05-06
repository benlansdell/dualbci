function UpdateLogSummary(name, matpath)

% UpdateLogSummary

% Run after running UpdateMatFiles if there are new .csv or .log files
% Unlike UpdateMatFiles, this script has to trash all previous results and
% re-traverse every .mat file matching the animal name.

fSummaryName = strcat(matpath,filesep,name,'_UserLogSummary.txt');


% Get all unique files with fnamebase in the name, stripping off extensions
fileList = GetFilesByName(name,matpath);



%% open messages summary file, clearing existing contents
fid = fopen(fSummaryName,'w+');

for idx = 1:length(fileList)
        base = fileList{idx};
        load(strcat(matpath,filesep,base,'.mat'));
        fprintf(strcat(base,'\n')); % status info
        
        % Pull out textual log messages and add to summary text file
        fprintf(fid,strcat('\n-------------------------',base,'--------------------\n'));
        for i = 1:length(data.messages.text)
            fprintf(fid, strcat(data.messages.text{i},'\n'));
        end
        
        clear data
end

fclose(fid);