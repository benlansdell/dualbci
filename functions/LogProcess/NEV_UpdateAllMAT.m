% Adds NEV catalog to every MAT file in path matching animal name

% USAGE:    NEV_UpdateAllMAT(name, matpath, nevpath)
%
% INPUTS:   name    animal name

function ME = NEV_UpdateAllMAT(name, matpath, nevpath)

global matpath nevpath metaData

if (nargin < 1) name = 'Spanky'; end
if (nargin < 2) matpath = '/home/lansdell/projects/bci/matlab/labview/'; end
if (nargin < 3) nevpath = '/home/lansdell/projects/bci/matlab/blackrock/'; end

fn_in = '/home/lansdell/projects/bci/matlab/nevmappings.csv_needsprocessing';
fn_out_success = '/home/lansdell/projects/bci/matlab/nevmappings.csv_success';
fn_out_fail = '/home/lansdell/projects/bci/matlab/nevmappings.csv_fail';

% Quick sanity check before we get to business.
assert(isdir(matpath) && isdir(nevpath),'Valid .nev/.mat directories not provided, aborting.')
assert(~isempty(metaData),'You didn''t load a metaData file...')
cd(matpath)
%fileList = GetFilesByName(name,matpath)
fileList = GetUnprocessedFiles(fn_in);
%fileList = fileList(end:-1:1);      % Update most recent files first
ME = {};

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % first pass is only unvisited files
 %for idx = 1:length(fileList)
 %        fname = fileList{idx};
 %
  %       fprintf('Processing %s ...\n',fname )
   %      % parse the .log and .csv files into a matlab data structure, which
    %     % is written to a .mat file by this script (and left in workspace)
     %    load(fname)
      %   try
 %            if (data.nevUpdateNum > 0)
  %               continue
   %          end
    %     catch me    % this section is fallback conditional for unvisited 
 %            
  %           try
   %              NEV_FindTimesForMAT1b(data,matpath,nevpath,metaData,0,0);
    %         catch me
 %                me.getReport()
  %               ME = [ME, me];
   %          end
    %         clear data    
 %        end
 %end
 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % second pass is files not updated since most recent script update

fh_success = fopen(fn_out_success, 'w');
fh_fail = fopen(fn_out_fail, 'w');
for idx = 1:length(fileList)
        fname = fileList{idx};
        fname = [matpath '/' fname];

        fprintf('Processing file %d of %d: %s ...\n',idx, length(fileList), fname )
        % parse the .log and .csv files into a matlab data structure, which
        % is written to a .mat file by this script (and left in workspace)
        load(fname)
        try
            [data, success] = NEV_FindTimesForMAT1b(data,matpath,nevpath,metaData,0,0);
            if success == 1
              fprintf('Counting as success.\n');
              fprintf(fh_success, '%s succeeded.\n', fname);
            else
              fprintf('Counting as failure.\n');
              fprintf(fh_fail, '%s failed.\n', fname);
            end
        catch me
            fprintf('Counting as failure.\n');
            fprintf(fh_fail, '%s failed.\n', fname);
            me.getReport()
            ME = [ME, me];
        end
        clear data
end
fclose(fh_fail);
fclose(fh_success);