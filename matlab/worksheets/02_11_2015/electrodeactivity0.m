duration = 360; %at least six minutes
fn_out = './worksheets/01_17_2015/weeklynevs.mat';
conds = {'2D Manual Position', 'Dual Control'};
load(fn_out);
binsize = 1/60;
threshold = 0;
offset = 0;
const = 'on';
pval = 0.001;
nK_sp = 6;
nK_pos = 6;
nE = 128;
nC = length(conds);

activity = {};

%For each condition and for each file in each condition:
for i = 1:nC
    condition = conds{i};
    nevs = weeklynevs(condition);
    nR = length(nevs);
    activity{i} = zeros(nE, nR);
    for j = 1:nR
        recinfo = nevs(j);
        %Get files
        nevfile = ['./blackrock/' recinfo.nevfile];
        matfile = ['./labview/' recinfo.matfile];
        currdate = strrep(recinfo.date, '/', '-');
        if (exist(nevfile, 'file') == 2) & (exist(matfile, 'file') == 2)
            %Preprocess data
            processed = preprocess_spline_lv(nevfile, matfile, binsize, threshold, offset);     
            %Combine data to same electrode
            processed_mua = combine_mua(processed);
            %Look at firing rate of each electrode
            totalspks = sum(processed_mua.binnedspikes,1);
            dur = size(processed_mua.binnedspikes,1)*processed_mua.binsize;
            fr = totalspks/dur;
            %Add to activity matrix for this condition
            electrodes = cellfun(@str2num, processed_mua.unitnames);
            activity{i}(electrodes, j) = fr;
        else
            display(['Missing files: ' nevfile ' or ' matfile '. Continuing to next file.'])
        end
    end
    %Plot
    clf;
    fn_out = ['./worksheets/02_11_2015/plots/' condition '_avefiring.eps'];
    colormap(hot);
    imagesc(activity{i});
    title('Averate firing rate per recording (spikes/s)');
    ylabel('Electrode');
    xlabel('Week');
    %set(gca,'XTick',1:nU);
    %set(gca,'YTick',1:nU);
    %set(gca,'XTickLabel',processed.unitnames);
    %set(gca,'YTickLabel',processed.unitnames);
    %rotateXLabels(gca, 90);
    colorbar;
    saveplot(gcf, fn_out);
end

save(['./worksheets/02_11_2015/electrodeactivities.mat'], 'activity');