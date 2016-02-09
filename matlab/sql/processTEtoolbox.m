function processTEtoolbox(conn, modelID, blackrock, labviewpath, MCnevfile1, BCnevfile1, MCnevfile2, DCnevfile, expt_id, paramcode)
	%Keep the units of the network the same between recordings
	%Load parameters
	eval(paramcode);

	%Figure out units to use... above 5hz in all recordings
	units = fetch(exec(conn, ['SELECT u1.unit FROM '...
	'`experiment_tuning` et1 '...
	'INNER JOIN `units` u1 '...
	'ON u1.`nev file` = et1.`manualrecording` '...
	'INNER JOIN `units` u2 '...
	'ON u2.`nev file` = et1.`1DBCrecording` '...
 	'INNER JOIN `units` u4 '...
	'ON u4.`nev file` = et1.`dualrecording` '...
	'WHERE u1.unit = u2.unit AND u1.unit = u4.unit AND '...
	'u1.firingrate > ' num2str(threshold) ' AND u2.firingrate > ' num2str(threshold)  ' AND '...
	'u4.firingrate > ' num2str(threshold) ' AND et1.`manualrecording` = "' MCnevfile1 '"']));
	otherunits = units.Data;
 
	%Make sure BC units from both dual and brain control are in there
	bciunits = exec(conn, ['SELECT `unit` FROM bci_units WHERE `ID` = "' BCnevfile1 '"']);
	bciunits = fetch(bciunits);
	bciunits = bciunits.Data;
	dualunits = exec(conn, ['SELECT `unit` FROM bci_units WHERE `ID` = "' DCnevfile '"']);
	dualunits = fetch(dualunits);
	dualunits = dualunits.Data;
	allunits = unique([otherunits; bciunits; dualunits]);

	%Tag with computer run on, date, last git commit
	host = hostname();
	stamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
	comm = currCommit();

	previous = fetch(exec(conn, ['SELECT id FROM analyses WHERE `experiment_id` = "' num2str(expt_id) '" AND modelID = ' num2str(modelID)]));
	if ~strcmp(previous.Data{1}, 'No Data')
		display(['modelID ' num2str(modelID) ' and experiment_id ' num2str(expt_id) ' already analysed. Skipping'])
		return
	end

	%Setup an analysis
	%Insert into analyses
	tablename = 'analyses';
	fitcols = {'modelID', '`experiment_id`', 'unit', 'unitnum', 'ncoeff', 'computer', '`analysis date`', 'commit'};
	sqldata = { modelID, expt_id, 'NULL', 1, 1, host, stamp, comm};
	datainsert(conn,tablename,fitcols,sqldata);
	%Get the analysis_id used
	analysis_id = fetch(exec(conn, 'SELECT LAST_INSERT_ID()'));
	analysis_id = analysis_id.Data{1};

	runTE(conn, analysis_id, modelID, blackrock, labviewpath, MCnevfile1, allunits, expt_id, paramcode);
%	runTE(conn, analysis_id, modelID, blackrock, labviewpath, MCnevfile2, allunits, expt_id, paramcode);
	runTE(conn, analysis_id, modelID, blackrock, labviewpath, BCnevfile1, allunits, expt_id, paramcode);
	runTE(conn, analysis_id, modelID, blackrock, labviewpath, DCnevfile, allunits, expt_id, paramcode);
end 

function runTE(conn, analysis_id, modelID, blackrock, labviewpath, nevfile, units, expt_id, paramcode)
	nevpath = [blackrock nevfile];
	%Load parameters
	eval(paramcode);
	nU = length(units);

	%Preprocess data(nevfile, binsize, threshold, offset, fn_out, verbose, units)
	processed = preprocess_spline(nevpath, binsize, threshold, offset, [], [], units);
	%Truncate to specified duration
	processed = truncate_recording(processed, dur);
	%Run with position filters
    
	%Estimates the Transfer Entropy using the TE toolbox from Shinya Ito, Indiana University
    
    % * transent.c  must be compiled prior to this
    asdf = SparseToASDF(processed.binnedspikes', 1);
    %check if processed or processed' needs to be used
    [te_estimate, CI, TEdelays_big]=ASDFTE(asdf,j_delay, i_order, j_order);
    
   
	%Tag with computer run on, date, last git commit
	host = hostname();
	comm = currCommit();
	stamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
	%For each unit, save the results 
	for idx = 1:nU
		%Extract and save regression fiticients
		unit = processed.unitnames{idx};
	
		%If already in database, skip
		%unit = '21.3'; modelID = 2; nevfile = '20140610SpankyUtah002.nev';
		previous = fetch(exec(conn, ['SELECT id FROM fits WHERE `analyses_id` = ' num2str(analysis_id) ' AND `nev file` = "' nevfile '" AND modelID = ' num2str(modelID) ' AND unit = "' unit '"']));
		if ~strcmp(previous.Data{1}, 'No Data')
			display(['Model ' num2str(modelID) ' nevfile ' nevfile ' and unit ' unit ' already analysed. Skipping'])
			continue
		end

		tablename = 'fits';
		fitcols = {'modelID', '`analyses_id`', '`nev file`', 'unit', 'unitnum', 'computer', '`analysis date`', 'commit'};
		sqldata = { modelID, analysis_id, nevfile, unit, idx, host, stamp, comm};
		%sqldata = { 1, '20130920SpankyUtah001.nev', 999, 1, 3, 3, '3', '2013-12-09 12:12:12', '12'};
		datainsert(conn,tablename,fitcols,sqldata);
		%Get the fit id used
		fitid = fetch(exec(conn, 'SELECT LAST_INSERT_ID()'));
		fitid = fitid.Data{1};

		%Insert into fits_te
		tablename = 'fits_te';
		fitcols = {'id', 'j_delay', 'i_order', 'j_order'};
		sqldata = { fitid, j_delay', i_order, j_order};
		datainsert(conn,tablename,fitcols,sqldata);

		%Insert into NetworkEstimates
		tablename = 'estimates_te';
		fitcols = {'id', 'fromnum', 'fromunit', 'score'};
		for j = 1:nU
			if j ~= idx
				unitj = processed.unitnames{j};
				te = te_estimate(idx, j);		
				sqldata = { fitid, j, unitj, te};
				datainsert(conn,tablename,fitcols,sqldata);
			end
		end
	end	
end 
