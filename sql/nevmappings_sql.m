function np = nevmappings_sql(matpath)
	%Produces csv file listing all processed nev files with date, time,
	%duration, and mapping used
	%
	%Usage:
	%	nevmappings_sql(matpath)
	%     
	%Input:
	%	matpath = (optional, default = './labview/') path to processed
	%		labview .mat files
	%  
	%Test code:
	%	matpath = './labview/';
	%	nevmappings(matpath);

	if (nargin < 1) matpath = './labview/'; end
	%Set up database connection
	conn = database('',databaseuser,databasepwd,'com.mysql.jdbc.Driver', ...
	databaseurl);
	tablename = 'Recordings2';
	colnames = {'`labview file`', '`labview date`', 'sessionID', '`nev file`', '`nev date`',...
	 'starttime', 'endtime', 'duration', 'tasktype', 'axis', 'map', 'Toffset', 'ch1', ...
				'ch2', 'ch3', 'ch4', '`labview task type`', '`confirmedmap`'};
	%Find all .mat files in matpath directory
	matfiles = dir([matpath '/*.mat']);
	nrows = 0;
	nnp = 0;
	newrows = {};
	np = {};
	for i = 1:length(matfiles)
		fn = matfiles(i).name
		matfile = fn;
		matdate = [matfile(8:17) ' ' matfile(19:20) ':' matfile(21:22) ':00' ];

		%Only process dates after Dec 5 2013
		mattime = datenum(matdate);
		dec = datenum('2013-12-08');
		if mattime < dec
			continue
		end

		%Don't process Log metadata files
		if ~isempty(findstr(fn, 'Log'))
			continue
		end

		%Check if mat file exists already in database, if so then skip the file
		alreadyadded = fetch(exec(conn, ['SELECT * FROM ' tablename ' WHERE `labview file` = "' fn  '"']));
		alreadyadded = alreadyadded.Data;
		if ~strcmp(alreadyadded, 'No Data')
			display([fn ' already added to table'])
			continue
		end
		%Then, for each mat file, load it 
		%extract the recording info from it, write to table, clear structure
		clear data;
		load([matpath '/' fn]);
		needsprocessing = 1;
		if exist('data', 'var') 
			if isfield(data, 'nev')
			if isfield(data.nev, 'map')
			needsprocessing = 0;
			%Find nev file name, mat file name, date-time, duration, mapping
			for j = 1:length(data.nev)
				Toffset = data.nev(j).Toffset(1);
				nevfile = data.nev(j).nevfile;
				nevdate = datestr(addtodate(datenum(matdate), Toffset/60, 'second'), 'yyyy-mm-dd hh:MM:ss');
				nevdur = data.nev(j).DurationSec;
				starttime = Toffset/60;
				endtime = starttime + nevdur;
				labviewtasktype = data.nev(j).map;
				tasktype = 'NUL';
				if ~isempty(findstr(labviewtasktype, 'Manual'))
					tasktype = 'manual';
				elseif ~isempty(findstr(labviewtasktype, 'Brain'))
					tasktype = 'brain';
				elseif ~isempty(findstr(labviewtasktype, 'Dual'))
					tasktype = 'dual';
				elseif ~isempty(findstr(labviewtasktype, 'Multi'))
					tasktype = 'multi';
				elseif ~isempty(findstr(labviewtasktype, 'wait'))
					tasktype = 'wait';
				end
				map = 'NUL';
				if ~isempty(findstr(labviewtasktype, 'Vel'))
					map = 'vel';
				elseif ~isempty(findstr(labviewtasktype, 'Pos'))
					map = 'pos';
				end
				axis = 'NUL';
				if ~isempty(findstr(labviewtasktype, 'Horiz'))
					axis = 'horiz';
				elseif ~isempty(findstr(labviewtasktype, 'Vert'))
					axis = 'vert';
				elseif ~isempty(findstr(labviewtasktype, '2D'))
					axis = '2D';
				end
				ch1 = data.nev(j).chans(1);
				ch2 = data.nev(j).chans(2);
				ch3 = data.nev(j).chans(3);
				ch4 = data.nev(j).chans(4);
				confirmed = 0;
				if strcmp(labviewtasktype, 'Error') | strcmp(labviewtasktype, '') | strcmp(nevfile, '') | strcmp(labviewtasktype, 'Need to redo this using within-data map info')
					needsprocessing = 1;
				end
				%Add to database
				%fprintf(fh, '%s,%s,%s,%f,%s\n', nevfile, matfile, nevdate, nevdur, nevmap);
				nrows = nrows + 1;
				%colnames = {'`labview file`', '`labview date`', 'sessionID', '`nev file`', '`nev date`',...
	 			%'starttime', 'endtime', 'duration', 'tasktype', 'axis', 'map', 'Toffset', 'ch1', ...
				%'ch2', 'ch3', 'ch4', '`labview task type`', '`confirmed map`'};
				if length(nevfile) > 0
					newrows = vertcat(newrows, {matfile, matdate, j, nevfile, nevdate, starttime, endtime,...
						nevdur, tasktype, axis, map, Toffset, ch1, ch2, ch3, ch4, labviewtasktype,...
						confirmed});
				end
			end
			end
			end
		end
		if needsprocessing == 1
			%fprintf(fh_process, '%s\n', fn);
			display([fn ' needs to be processed by Charlies nev map code']);
			nnp = nnp + 1;
			np{nnp} = fn;
			system(['cp ' matpath '/' fn ' ~/projects/bci/unprocessedmat']);
		end
	end
	fclose(fh_process);
	display('Adding all rows to database');
	datainsert(conn,tablename,colnames,newrows)
end 

%LogMetaData2('spanky_db', 0, true)
%global matpath
%global nevpath
%global metaData
%matpath = '../unprocessedmat'; nevpath = './blackrock';
%load([matpath '/Spanky_LogMetaData2.mat']);
%anp = dir([matpath '/*2014*.mat']);
%
%for idx = 1:length(anp)
%	idx
%	fn = anp(idx).name
%	%system(['cp ' matpath '/' fn ' ~/projects/bci/unprocessedmat']);
%	load([matpath '/' fn]);
%	NEV_FindTimesForMAT4(data)
%end