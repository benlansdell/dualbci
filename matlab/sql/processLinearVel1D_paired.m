function processLinearVel1D_Paired(conn, modelID, blackrock, labviewpath, MCnevfile1, BCnevfile1, MCnevfile2, DCnevfile, expt_id, paramcode)
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
	'INNER JOIN `units` u3 '...
	'ON u3.`nev file` = et1.`manualrecordingafter` '...
	'INNER JOIN `units` u4 '...
	'ON u4.`nev file` = et1.`dualrecording` '...
	'WHERE u1.unit = u2.unit AND u2.unit = u3.unit AND u3.unit = u4.unit AND '...
	'u1.firingrate > ' num2str(threshold) ' AND u2.firingrate > ' num2str(threshold) ' AND u3.firingrate > ' num2str(threshold) ' AND '...
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

	%Find BC axis (1 = RU, 2 = FE)
	bciaxis = exec(conn, ['SELECT `direction` FROM bci_units WHERE `ID` = "' BCnevfile1 '"']);
	bciaxis = fetch(bciaxis);
	bciaxis = bciaxis.Data;
	if strcmp(bciaxis{1}, 'east') | strcmp(bciaxis{1}, 'west')
		%We want to fit the model to the _non_ BC axis, so, if BCI = FE = 'east/west' = 2 
		%we'll fit the model to RU = 1
		axe = 1;
	else
		axe = 2;
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

	runLinearPaired(conn, analysis_id, modelID, blackrock, labviewpath, MCnevfile1, allunits, expt_id, paramcode, axe);
	runLinearPaired(conn, analysis_id, modelID, blackrock, labviewpath, MCnevfile2, allunits, expt_id, paramcode, axe);
	runLinearPaired(conn, analysis_id, modelID, blackrock, labviewpath, BCnevfile1, allunits, expt_id, paramcode, axe);
	runLinearPaired(conn, analysis_id, modelID, blackrock, labviewpath, DCnevfile, allunits, expt_id, paramcode, axe);
end 

function runLinearPaired(conn, analysis_id, modelID, blackrock, labviewpath, nevfile,    units,    expt_id, paramcode, axe)
	%nevfile = MCnevfile1;
	%units = allunits;

	nevpath = [blackrock nevfile];
	%Load parameters
	eval(paramcode);
	nU = length(units);

	%Preprocess data(nevfile, binsize, threshold, offset, fn_out, verbose, units)
	processed = preprocess_spline(nevpath, binsize, threshold, offset, [], [], units);
	%Truncate to specified duration
	%processed = truncate_recording(processed, dur);
	[processed, processed_novel] = split_recording(processed, dur, dur+testdur);

	%Fit a linear model to training data
	d = ['processLinearVel1D: Fitting model ' num2str(modelID) ' nevfile ' nevfile];
	display(d);
	%writelog(logfile, d);
	data = filters_sp_vel_1D(processed, nK_sp, nK_pos, axe);
	model = MLE_lmfit(data, const);
	%%Compute MSE on test data
	datanovel = filters_sp_vel_1D(processed_novel, nK_sp, nK_pos, axe);
	nBout = size(datanovel.y,2);
	nBin = size(data.y,2);

	if length(size(data.X)) == 2
		nC = 1
	else
		nC = size(data.X,2);
	end
	if strcmp(const, 'on')
		nC = nC + 1;
	end
	%Tag with computer run on, date, last git commit
	host = hostname();
	comm = currCommit();
	stamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
	%For each unit, save the results 
	for idx = 1:nU
		%Extract and save regression coefficients
		unit = processed.unitnames{idx};
		%Extract deviance
		dev = model.dev{idx, 1};
		%Extract SE
		stats = model.stats{idx};
		%Extract MSE. Need to add cross-validation code to do this... add later
		b_hat = model.b_hat(idx,:);
		rho_hat = glmval(b_hat', squeeze(data.X(idx,:,:)), 'identity');
		rho_hat_novel = glmval(b_hat', squeeze(datanovel.X(idx,:,:)), 'identity');

		mseout = sum((datanovel.y(idx,:)-rho_hat_novel').^2)/nBout;
		meany = mean(data.y(idx,:));
		r2 = 1-sum((data.y(idx,:)-rho_hat').^2)/sum((data.y(idx,:)-meany).^2);
		meany_novel = mean(datanovel.y(idx,:));
		r2_novel = 1-sum((datanovel.y(idx,:)-rho_hat_novel').^2)/sum((datanovel.y(idx,:)-meany_novel).^2);

		%If already in database, skip
		%unit = '21.3'; modelID = 2; nevfile = '20140610SpankyUtah002.nev';
		previous = fetch(exec(conn, ['SELECT id FROM fits WHERE `analyses_id` = ' num2str(analysis_id) ' AND `nev file` = "' nevfile '" AND modelID = ' num2str(modelID) ' AND unit = "' unit '"']));
		if ~strcmp(previous.Data{1}, 'No Data')
			display(['Model ' num2str(modelID) ' nevfile ' nevfile ' and unit ' unit ' already analysed. Skipping'])
			continue
		end

		tablename = 'fits';
		fitcols = {'modelID', '`analyses_id`', '`nev file`', 'unit', 'unitnum', 'ncoeff', 'dev', 'computer', '`analysis date`', 'commit'};
		sqldata = { modelID, analysis_id, nevfile, unit, idx, nC, dev, host, stamp, comm};
		%sqldata = { 1, '20130920SpankyUtah001.nev', 999, 1, 3, 3, '3', '2013-12-09 12:12:12', '12'};
		datainsert(conn,tablename,fitcols,sqldata);
		%Get the fit id used
		fitid = fetch(exec(conn, 'SELECT LAST_INSERT_ID()'));
		fitid = fitid.Data{1};

		%Insert into fits_linear
		tablename = 'fits_linear';
		fitcols = {'id', 'dir', 'size', 'r2', 'r2out'};
		tuningsize = abs(model.b_hat(idx,2));
		if model.b_hat(idx,2) < 0
			
		else
			
		end
		[direction, tuningsize] = unitTheta(regrRU, regrFE);
		sqldata = { fitid, direction, tuningsize, r2, r2_novel};
		datainsert(conn,tablename,fitcols,sqldata);

		%Insert into estimates_parameters
		tablename = 'estimates_parameters';
		fitcols = {'id', 'num', 'label', 'value', 'se', 'mask'};
		for j = 1:nC
			num = j;
			if strcmp(const, 'on')
				labidx = j-1;
				if j == 1
					label = 'const';
				else
					label = findLabel(labidx, data.k);				
				end
			else
				labidx = j;
				label = findLabel(labidx, data.k);
			end
			mask = model.mask(idx,j);			
			sqldata = {fitid, num, label, model.b_hat(idx, j), stats.se(j), mask};
			datainsert(conn,tablename,fitcols,sqldata);
		end
	end	
end 
