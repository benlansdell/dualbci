
%%

function BootstrapScript(metaData,matpath,saveas)

% if (nargin == 0)
%     fprintf('Loading metadata, since you didn''t give it to me...\n')
%     load Spanky_LogMetaData
% end

Nmin = 10;

conds = {'1D.+Brain','1D.+Manual.+Velocity','2D.+Brain.+Velocity'};


try
    metaData.Models{1}.BS
    fprintf('Bootstrap analysis already done for %d of %d segments...\n',...
            sum(metaData.Models.BS), length(metaData.BS));
    
catch me        % if above throws error, assume BS doesn't exist
    for i = 1:length(metaData.Models)
    metaData.Models{i}.BS = false;
    end
end
clc

for cond = conds
    cond = cond{1}; % make string, not cell array
    
fprintf('\n\n****Bootstrapping analysis for condition: %s\n',cond)

stats = DailyStatsForCondition(metaData,cond);  % sorted by Shannon-Welford

%idxSet = [stats([stats.SWrsquare] > 0.90).idx];
idxSet = row([stats( ([stats.SWrsquare] > 0.75 | [stats.Frsquare] > 0.75) & ([stats.Wcount] > 1)).idx ]); 

for idx = idxSet
    model = metaData.Models{idx};
    if (model.BS)
        continue    % already did bootstrap analysis on this segment
    end
    fprintf('Condition: %s, %d of %d\n',cond, find(idxSet == idx), length(idxSet));
    fprintf('%s a=%2.1f\tb=%2.1f\tk=%2.1f\tR^2=%0.3f\n',model.Day,model.SWfitObj.a,model.SWfitObj.b,model.SWfitObj.k,model.SWgof.rsquare)
    
    metaData.Models{idx}.bootstrap = BootstrapAnalysis(model,500,Nmin,true);
    if (~isempty(metaData.Models{idx}.bootstrap))
        metaData.Models{idx}.BS = true;
        metaData.Models{idx}.SW
    end
    
end
fprintf('Saving metaData...\n')
save([matpath filesep saveas],'metaData')
fprintf('\n\n***** SAVE POINT ****\n\n')
end

%%  % day with R^2 = 0.98 for SW model
% load Spanky_2013-03-20-1330    


%% Manual control day with R^2 = 0.95
% model = metaData.Models{341}

