%Go through and add to GC matrices missing entries for electrodes that were not including in the GC analysis
%This is so that when plotted using cytoscape the plots from different weeks will have all electrodes, not just
%the analysed subset. The plots will look the same this way, and can be more easiliy compared.

%Open all .mat files in ./worksheets/01_17_2015/grangerresults

fn_out = './worksheets/01_17_2015/weeklynevs.mat';
conds = {'2D Manual Position', 'Dual Control'};
load(fn_out);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add bci units, if relevant%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(conds)
	condition = conds{i};
	nevs = weeklynevs(condition);
	%Only if under some form of brain control (not manual)
	if isempty(strfind(condition, 'Manual'))
		display(['Adding BCI units used in ' condition ' recordings.'])
		for j = 1:length(nevs)
			recinfo = nevs(j);
			%Get files
			nevfile = ['./blackrock/' recinfo.nevfile];
			matfile = ['./labview/' recinfo.matfile];
			currdate = strrep(recinfo.date, '/', '-');
			curr_fn = ['./worksheets/01_17_2015/grangerresults/' condition '_' currdate '.mat'];
			curr_fn = strrep(curr_fn, ' ', '_')
			if exist(curr_fn, 'file')
				display(['Adding BCI units to ' condition '_' currdate '.mat'])
				load(curr_fn)
				%Load the matfile
				load(matfile);
				%Find the nev file 
				for k = 1:length(data.nev)
					if strcmp(data.nev(k).nevfile , recinfo.nevfile)
						%Found the entry, grab the BCI units
						chans = data.nev(k).chans;
						onchans = floor(chans(chans>0));
						bciunits = arrayfun(@num2str, onchans, 'unif', 0);
					end
				end
				save(curr_fn, 'GCdev', 'GCpval', 'GCsig', 'processed_mua', 'bciunits', 'nevfile', 'matfile');
			end
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Extend to all electrodes, not just the ones above 5Hz%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

directory = './worksheets/01_17_2015/grangerresults/';
fns = dir([directory '*.mat']);
nE = 128;

for idx = 1:length(fns)
	fn = fns(idx).name;
	fullpath = [directory fn];
	%Load
	load(fullpath);
	%Manipulate stuff
	electrodes = cellfun(@str2num, processed_mua.unitnames);
	GCdevext = zeros(nE, nE);
	GCpvalext = zeros(nE, nE);
	GGsigext = zeros(nE, nE);
	GCdevext(electrodes,electrodes) = GCdev;
	GCpvalext(electrodes,electrodes) = GCpval;
	GCsigext(electrodes,electrodes) = GCsig;
	unitnames = processed_mua.unitnames;
	unitnamesext = arrayfun(@num2str, 1:nE, 'unif', 0);
	display(fullpath)
	%Save mat files for visualization in cytoscape
	%Need: GCdev, GCpval, GCsig, unitnames, brainunits (cell array of (string) BCI unit names)
	save(fullpath, 'GCdev', 'GCpval', 'GCsig', 'GCdevext', 'GCpvalext', 'bciunits', 'GCsigext', 'unitnames', 'unitnamesext', 'processed_mua');
end