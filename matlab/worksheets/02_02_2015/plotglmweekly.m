%Plot results from glmweekly for files that match a given pattern
duration = 360; %at least six minutes
fn_out = './worksheets/01_17_2015/weeklynevs.mat';
conds = {'2D Manual Position', 'Dual Control'};
load(fn_out);
binsize = 0.002;
dt_sp = binsize;
dt_pos = 0.2;
threshold = 5;
offset = 0;
const = 'on';
nK_sp = 100;
nK_pos = 6;
modeltype = 'sp_pos_lv';

devs = {};
topten = {};
nE = 128;

%Load each file and extract the deviances for each unit
for i = 1:length(conds)
	condition = conds{i};
	nevs = weeklynevs(condition);
	nW = length(nevs);
	devs{i} = zeros(nW, nE);
	topten{i} = zeros(nW, 1);
	for j = 1:nW
		recinfo = nevs(j);
		%Get files
		nevfile = ['./blackrock/' recinfo.nevfile];
		matfile = ['./labview/' recinfo.matfile];
		currdate = strrep(recinfo.date, '/', '-');
		curr_fn = ['./worksheets/02_02_2015/glmresults/' modeltype '_' condition '_' currdate '.mat'];
		curr_fn = strrep(curr_fn, ' ', '_')
		if exist(curr_fn, 'file')
			%Add data to matrix for plotting
			display([curr_fn ' exists, reading in deviances.'])
			load(curr_fn);
			nU = size(model.b_hat,1);
			%Extract deviances
			curr_devs = [];
			for k = 1:nU
				curr_devs(k) = model.dev{k};
			end
			%Find units used
			units = cellfun(@str2num, processed_mua.unitnames);
			%Add these to devs matrix
			devs{i}(j,units) = curr_devs;
			tt = sort(curr_devs, 'ascend');
			topten{i}(j) = sum(tt(1:8));
		end
	end
	%Once loaded all files, plot the results as a heatmap for all electrodes, 
    clf;
    fn_out = ['./worksheets/02_02_2015/plots/' modeltype '_' condition '.eps'];
    colormap(hot);
    imagesc(devs{i}');
    title(['Deviance under ' modeltype ' for ' condition ' recordings']);
    ylabel('Electrode');
    xlabel('Week');
    %set(gca,'XTick',1:nU);
    %set(gca,'YTick',1:nU);
    %set(gca,'XTickLabel',processed.unitnames);
    %set(gca,'YTickLabel',processed.unitnames);
    %rotateXLabels(gca,90);
    colorbar;
    saveplot(gcf, fn_out);
    %Look at sum of top 10 performing units for each recording
    clf
    plot(topten{i}, 'o')
    xlabel('Week')
    ylabel('Sum of smallest 8 deviances')
    fn_out = ['./worksheets/02_02_2015/plots/' modeltype '_' condition '_toptendev.eps'];
    saveplot(gcf, fn_out)
    %Look at boxplot of all units for each week
    clf
    [X,Y] = meshgrid(1:nW, 1:nE);
    grps = reshape(X, nW*nE,1);
    rshpdevs = reshape(devs{i}, nW*nE,1);
    boxplot(rshpdevs(rshpdevs>0), grps(rshpdevs>0));
    xlabel('Week')
    ylabel('Deviance')
    fn_out = ['./worksheets/02_02_2015/plots/' modeltype '_' condition '_boxplot.eps'];
    saveplot(gcf, fn_out)    
	%Save data for comparing to other conditions run separately
end