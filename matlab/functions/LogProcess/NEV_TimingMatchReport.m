% This script reports on problematic match failures of channels & time
% offsets between NEV files and mat files.

% it assumes that metaData is loaded
function NEV_TimingMatchReport(fid)
global metaData matpath

% results written to file by default, set fid = 1 for console output
if (nargin == 0)
    fid = fopen([matpath filesep '_NEV_TimingMatchReport.txt'],'w');
end


flagIdx = [];

for fi = 1:length(metaData.fileList)
    f = metaData.fileList{fi};
    load([matpath filesep f])
    
    try
        for ni = 1:(length(data.nev)-1)             % ignore first, last block
            if data.nev(ni).DurationSec >= 60*5    % skip waitfiles (~10-15sec)
                
                if isempty(data.nev(ni).segInfo)
                    fprintf(fid,'No trial block match for %s NEV #%d %.0fs\n',...
                        f, ni,data.nev(ni).DurationSec);
                    if ~ismember(fi,flagIdx)
                        flagIdx = [flagIdx, fi];
                    end
                end
                
                if diff(data.nev(ni).Toffset([2 3])) > 1    % 0.1s @ 60Hz
                    fprintf(fid,'Timing mismatch for %s NEV #%d %.0fs\n',...
                        f, ni, data.nev(ni).DurationSec);
                    if ~ismember(fi,flagIdx)
                        flagIdx = [flagIdx, fi];
                    end
                end
                if sum(data.nev(ni).chans == -1) > 0 && ...
                        ~isempty(regexp(data.nev(ni).map,'.*Brain.*')) && ...
                        sum(ismember(find(data.nev(ni).chans == -1),find(data.nev(ni).segInfo.enabled)))
                    fprintf(fid,'Failed channel IDs for %s NEV #%d %.0fs ch: %s\n',...
                        f, ni, data.nev(ni).DurationSec,...
                        num2str(find(data.nev(ni).chans == -1)));
                    if sum([data.nev(ni-1).chans data.nev(ni+1).chans] ~= -1)
                        fprintf(fid,'**Recovery possible\n');
                        if ~ismember(fi,flagIdx)
                            flagIdx = [flagIdx, fi];
                        end
                    end
                end
            end
            
        end
    catch me
        %fprintf('Failed to read NEV data from %s\n',f)
        %me.getReport()
    end
    fprintf('.\n')
end

if fid ~= 1
    fclose(fid)
end
    