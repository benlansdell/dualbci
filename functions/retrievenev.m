function nevs = retrievenev(csv, year, month, day, condition)
	%Produces csv file listing all processed nev files from a given date and condition
	%
	%Usage:
	%	nevs = retrievenev(csv, year, month, day, condition)
	%     
	%Input:
	%	csv = filename to csv file output from nevmappings containing all nev and mapping info
	%	year = integer or list of integers specifying year(s) to retrieve nev files
	%	month = integer or list of integers specifying month(s) to retrieve nev files
	%	day = integer or list of integers specifying day(s) of the month to retrieve nev files
	%  
	%Output:
	%	nevs = cell array containing a list of file names under requested conditions and duration
	%
	%Test code:
	%	csv = './nevmappings2.csv';
	%	year = 2013;
	%	month = 1;
	%	%First week's worth of recordings
	%	days = 1:8;
	%	condition = '2D Manual Position';
	%	nevs = retrievenev(csv, year, month, days, condition);

	%Spanky_2013-01-08-1334.mat
	fh = fopen(csv, 'r');
	tline = fgetl(fh);
	mappings = {};
	c = 0;
	counts = containers.Map;
	totaldur = containers.Map;
	totaldurhrs = containers.Map;
	cnt = 1;

	while ischar(tline)
		%20130108SpankyUtah012.nev,Spanky_2013-01-08-1334.mat,2013-01-08-1334,16.323467,1D Horiz Brain  Velocity
		%Split 
		words = strsplit(tline,',');
		datetime = words{3};
		dur = str2num(words{4});
		mapping = words{5};
		nevfile = words{1};
		matfile = words{2};
		curr_year = str2num(datetime(1:4));
		curr_month = str2num(datetime(6:7));
		curr_day = str2num(datetime(9:10));
		if ~ismember(mapping, mappings)
			counts(mapping) = 1;
			totaldur(mapping) = dur;
			totaldurhrs(mapping) = dur/60/60;
		else
			counts(mapping) = counts(mapping) + 1;
			totaldur(mapping) = totaldur(mapping) + dur;
			totaldurhrs(mapping) = totaldurhrs(mapping) + dur/60/60;
		end
		if ismember(curr_year,year) & ismember(curr_day,day) & ismember(curr_month,month) & strcmp(mapping, condition)
			nevs(cnt).nevfile = nevfile;
			nevs(cnt).dur = dur;
			nevs(cnt).matfile = matfile;
			nevs(cnt).date = [num2str(curr_year) '/' num2str(curr_month) '/' num2str(curr_day)];
			cnt = cnt + 1;
		end
		mappings = {mappings{:}, mapping};
		tline = fgetl(fh);
		c = c + 1;
%		if c > 10
%			break
%		end
	end
	%display('Mappings')
	%keys(counts)'
	%display('nev file counts')
	%values(counts)'
	%display('total recorded duration (s)')
	%values(totaldur)'
	%display('duration (hrs)')
	%values(totaldurhrs)'
	if ~exist('nevs', 'var')
		nevs = [];
	end
	fclose(fh);
end 