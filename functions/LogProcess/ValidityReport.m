% This script goes through a loaded LogMetaData file and reports on counts
% of valid trials on each experimental day

for i=1:length(metaData.trials)
    if (~isempty(metaData.trials{i}))
        valid = sum(metaData.trials{i}.valid);
        count = length(metaData.trials{i}.valid);
        if (valid/count < 0.95)
        fprintf('%s: %d valid of %d total\n',...
            metaData.fileList{i},...
            valid,count)
        end
    end
end