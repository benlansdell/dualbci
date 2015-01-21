%Granger analysis on weekly dataset
duration = 360; %at least six minutes
fn_out = './worksheets/01_17_2015/weeklynevs.mat';
conds = {'2D Manual Position',
    'Dual Control'}
load(fn_out);
binsize = 1/60;
threshold = 5;
offset = 0;
const = 'on';
pval = 0.001;
nK_sp = 6;
nK_pos = 6;

%For each condition and for each file in each condition:
for i = 1:length(conds)
	condition = conds{i};
	nevs = weeklynevs(condition);
	for j = 1:length(nevs)
		recinfo = nevs(j);
		%Get files
		nevfile = ['./blackrock/' recinfo.nevfile];
		matfile = ['./labview/' recinfo.matfile];
		currdate = strrep(recinfo.date, '/', '-');
		%Preprocess data
		processed = preprocess_spline_lv(nevfile, matfile, binsize, threshold, offset);		
		%Combine data to same electrode
		processed_mua = combine_mua(processed);
		%Run with position filters_sp_pos_network
		fn_out = ['./worksheets/01_17_2015/plots/granger_pos_' condition '_' currdate '.eps'];
		data = filters_sp_pos_network_lv(processed_mua, nK_sp, nK_pos);
		[GCdev, GCpval, GCsig] = granger(processed, data, fn_out, pval);		
		%Save the GC matrices for producing GC plots in cytoscape
		fn_out = ['./worksheets/01_17_2015/grangerresults/' condition '_' currdate '.mat'];
		save(fn_out, 'GCdev', 'GCpval', 'GCsig', 'processed_mua');
	end
end