% midx is index of metaData entry. 
% This function will write [index, distance, width, movement time, map] to
% a csv file named the same as the .mat file, and put it in the .mat
% folder.

function MetaDataToCSV(midx)

if (nargin == 0)
    midx = 116;
end

try
global metaData matpath

trials = metaData.trials{midx};
fname = metaData.fileList{midx};

headers = {'Index','Distance','Width','Movement Time','Condition'};

% strips off any filename extension and replaces with .csv
fid = fopen(strcat(matpath,filesep,fname,'_TrialSummary.csv'),'w+');    

% write headers to file
for i = 1:length(headers)-1
    fprintf(fid,'%s, ',headers{i});
end
fprintf(fid,'%s\n',headers{end});


for i=1:length(trials.index)
    
    fprintf(fid,'%d, %1.2f, %1.2f, %2.2f, %s\n',trials.index(i),trials.D(i),trials.W(i),trials.MT(i),trials.Map{i});
    
end

fclose(fid);
catch me
    me.getReport
    keyboard
end
end