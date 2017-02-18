function processGPFAPaired(conn, modelID, blackrock, labviewpath, labviewfile, MCnevfile1, BCnevfile1, DCnevfile, MCnevfile2, expt_id, paramcode, basename)
	%Keep the units of the network the same between recordings
	%Load parameters
	eval(paramcode);
	%threshold = 5;

	%Figure out units to use... above 5hz in all recordings
	units = fetch(exec(conn, ['SELECT u1.unit FROM '...
	'`experiment_tuning` et1 '...
	'INNER JOIN `units` u1 '...
	'ON u1.`nev file` = et1.`manualrecording` '...
	'INNER JOIN `units` u2 '...
	'ON u2.`nev file` = et1.`1DBCrecording` '...
	'INNER JOIN `units` u3 '...
	'ON u3.`nev file` = et1.`manualrecordingafter` '...
	'INNER JOIN `units` u4 '...
	'ON u4.`nev file` = et1.`dualrecording` '...
	'WHERE u1.unit = u2.unit AND u1.unit = u4.unit AND u1.unit = u3.unit AND '...
	'u3.firingrate > ' num2str(threshold) ' AND ' ...
	'u1.firingrate > ' num2str(threshold) ' AND u2.firingrate > ' num2str(threshold) ' AND '...
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

	runGPFAPaired(conn, analysis_id, 1, modelID, blackrock, labviewpath, labviewfile, MCnevfile1, allunits, expt_id, paramcode, basename);
	runGPFAPaired(conn, analysis_id, 2, modelID, blackrock, labviewpath, labviewfile, BCnevfile1, allunits, expt_id, paramcode, basename);
	runGPFAPaired(conn, analysis_id, 3, modelID, blackrock, labviewpath, labviewfile, DCnevfile, allunits, expt_id, paramcode, basename);
	runGPFAPaired(conn, analysis_id, 4, modelID, blackrock, labviewpath, labviewfile, MCnevfile2, allunits, expt_id, paramcode, basename);
end 

function runGPFAPaired(conn, analysis_id, id, modelID, blackrock, labviewpath, labviewfile, nevfile, units, expt_id, paramcode, basename)
	nevpath = [blackrock nevfile];
	lbpath = [labviewpath labviewfile];
	%Load parameters
	eval(paramcode);
	nU = length(units);

	processed = preprocess_spline_target(nevpath, lbpath, binsize, threshold, offset, [], [], units);
	processed = truncate_recording(processed, dur);
	%Time bins which occur within a trial
	spikes = processed.binnedspikes;
	trials = sum(abs(processed.target),2)>0;

	%Run all trials together
	[dat, octs, quads] = makedata(trials, spikes, processed.cursor, processed.target);
	method = 'gpfa';
	% Select number of latent dimensions
	kernSD = 30;
	% Extract neural trajectories
	runIdx = analysis_id*10 + id;
	result = neuralTraj(runIdx, dat, 'method', method, 'xDim', ndim, 'kernSDList', kernSD);
	[estParams, seqTrain, seqTest, DD] = postprocess(result, 'kernSD', kernSD);
	%Save results
	fn_out = [basename '_' num2str(analysis_id) '_' num2str(id) '.mat'];
	save(fn_out, 'estParams', 'DD');
end 
