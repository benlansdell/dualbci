function nevmappings(fn_out, matpath)
	%Produces csv file listing all processed nev files with date, time,
	%duration, and mapping used
	%
	%Usage:
	%	nevmappings(fn_out, matpath)
	%     
	%Input:
	%	fn_out = (optional, default = './nevmappings.csv') filename to
	%		write data to in csv
	%	matpath = (optional, default = './labview/') path to processed
	%		labview .mat files
	%  
	%Test code:
	%	matpath = './testdata/';
	%	fn_out = './nevmappings.csv';
	%	nevmappings(fn_out, matpath);

	if (nargin < 1) fn_out = './nevmappings.csv'; end
	if (nargin < 2) matpath = './labview/'; end

	fh = fopen(fn_out, 'w');
	fh_process = fopen([fn_out '_needsprocessing'], 'w');
	fprintf(fh, 'nev file, mat file, date-time, duration, mapping\n');
	%Find all .mat files in matpath directory
	matfiles = dir([matpath '/*.mat']);
	for i = 1:length(matfiles)
		%For each mat file, load it, see if it has an nev structure, 
		%extract the recording info from it, write file, clear structure
		fn = matfiles(i).name
		clear data;
		load([matpath '/' fn]);
		needsprocessing = 1;
		if exist('data', 'var') 
			if isfield(data, 'nev')
			if isfield(data.nev, 'map')
			needsprocessing = 0;
			%Find nev file name, mat file name, date-time, duration, mapping
			for j = 1:length(data.nev)
				nevfile = data.nev(j).nevfile
				matfile = fn;
				nevdur = data.nev(j).DurationSec;
				nevmap = data.nev(j).map;
				if strcmp(nevmap, 'Error') | strcmp(nevmap, '') | strcmp(nevfile, '')
					needsprocessing = 1;
				end
				nevdate = strrep(data.fname, 'Spanky_', '');
				fprintf(fh, '%s,%s,%s,%f,%s\n', nevfile, matfile, nevdate, nevdur, nevmap);
			end
			end
			end
		end
		if needsprocessing == 1
			fprintf(fh_process, '%s\n', fn);
		end
	end
	fclose(fh);
	fclose(fh_process);
end 