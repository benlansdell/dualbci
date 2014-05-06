% This function traverses a log meta data file and dumps out a csv
% spreadsheet with the number of trials done with each mapping type, on
% each day.

function LogMetaData_TrialTypeSummary(name,validOnly)
global matpath

if (nargin < 2)
    validOnly = true;
end

if (nargin == 0)
    name = 'Spanky'
end

fileName = strcat(matpath,filesep,name,'_LogMetaData.mat');
load(fileName);

headers = {'Day','Neurons','1D Horiz Manual Position','2D Manual Position','2D Manual Goofy','1D Horiz Manual Velocity','1D Vertical Manual Velocity','2D Manual Velocity','1D Horiz Brain  Velocity','2D Brain        ','Dual Control','Unknown'};

% strips off any filename extension and replaces with .csv
fid = fopen(strcat(matpath,'/',name,'_TrialSummary.csv'),'w+');            

% write headers to file
for i = 1:length(headers)-1
    fprintf(fid,'%s, ',headers{i});
end
fprintf(fid,'%s\n',headers{end});

maps = headers(3:end);

% iterate over each day, counting valid trials under each mapping
for i = 1:length(metaData.trials)
    try
    if (~isempty(metaData.trials{i}) && ...
            sum(metaData.trials{i}.valid))
        sprintf(metaData.fileList{i})
        fprintf(fid,'%s, ',metaData.fileList{i});
        
        neurons = length(setdiff(unique(metaData.trials{i}.Neurons),[-1,0])); % -1 for inevitable 0
        fprintf(fid,'%d, ',neurons);
        
        for j = 1:length(maps)
            count = strcmp(metaData.trials{i}.Map, maps{j});
            if (validOnly)
                count = count & metaData.trials{i}.valid';
            end
            count = sum(count);
            if (j < length(maps))
                fprintf(fid,'%d, ',count);
            else
                fprintf(fid,'%n',count);
            end
        end
        fprintf(fid,'\n');
    end
    catch me
        keyboard
    end
end

fclose(fid);

end