% Adds NEV catalog to every MAT file in path matching animal name

% USAGE:    NEV_UpdateAllMAT(name, matpath, nevpath)
%
% INPUTS:   name    animal name
%
%


function ME = NEV_UpdateAllMAT(name)

global metaData matpath nevpath

% Quick sanity check before we get to business.
assert(isdir(matpath) && isdir(nevpath),'Valid .nev/.mat directories not provided, aborting.')
assert(~isempty(metaData),'You didn''t load a metaData file...')
cd(matpath)
fileList = GetFilesByName(name,matpath);
%fileList = fileList(end:-1:1);      % Update most recent files first
ME = {};

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % first pass is only unvisited files
% for idx = 1:length(fileList)
%         fname = fileList{idx};
% 
%         fprintf('Processing %s ...\n',fname )
%         % parse the .log and .csv files into a matlab data structure, which
%         % is written to a .mat file by this script (and left in workspace)
%         load(fname)
%         try
%             if (data.nevUpdateNum > 0)
%                 continue
%             end
%         catch me    % this section is fallback conditional for unvisited 
%             
%             try
%                 NEV_FindTimesForMAT1b(data,matpath,nevpath,metaData,0,0);
%             catch me
%                 me.getReport()
%                 ME = [ME, me];
%             end
%             clear data    
%         end
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % second pass is files not updated since most recent script update

for idx = 90:length(fileList)
        fname = fileList{idx};

        fprintf('Processing %s ...\n',fname )
        % parse the .log and .csv files into a matlab data structure, which
        % is written to a .mat file by this script (and left in workspace)
        load(fname)

            try
                NEV_FindTimesForMAT1b(data,matpath,nevpath,metaData,0,0);
            catch me
                me.getReport()
                ME = [ME, me];
            end
            clear data
end
