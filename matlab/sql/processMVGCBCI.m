function processMVGCBCI(conn, modelID, blackrock, labviewpath, nevfile, paramcode, threshold)
	nevpath = [blackrock nevfile];
	%Load parameters
	eval(paramcode);
	%Get which units are BCI units, get the mat file for this nev file
	bcunits = fetch(exec(conn, ['SELECT `unit` FROM `bci_units` WHERE `ID` = "' nevfile '"']));
	%bcunits = num2str(cell2mat(bcunits.Data));
	bcunits = bcunits.Data;
	if size(bcunits,1) < 2
		display(['Warning: fewer than 2 BCI units are labeled within ' nevfile])
	end
	matfile = fetch(exec(conn, ['SELECT `labview file`,`duration` FROM `recordings` rec WHERE rec.`nev file` = "' nevfile '"']));
	duration = matfile.Data{2};
	matfile = matfile.Data{1};
	matpath = [labviewpath matfile];
	taskaxes = fetch(exec(conn, ['SELECT `axis` FROM `recordings` rec WHERE rec.`nev file` = "' nevfile '"']));
	taskaxes = taskaxes.Data{1};
	%Get how many units are above threshold from units
	%%Still processing...
	%abovethresh = fetch(exec(conn, ['SELECT `unit` FROM `units` WHERE `nev file` = "' nevfile '" AND `firingrate` > ' num2str(threshold)]))
	%abovethresh = abovethresh.Data;
	if duration < dur
		display(['Recording is shorter than specified ' num2str(dur) ' seconds. Skipping'])
		return
	end

	%%Preprocess data
	%processed = preprocess_smooth_lv(nevpath, matfile, binsize, sigma_fr, sigma_trq, threshold, offset, [], [], units);
	processed = preprocess_spline_lv(nevpath, matfile, binsize, threshold, offset);
	%Truncate to specified duration
	processed = truncate_recording(processed, dur);
	nUabove = length(processed.unitnames);
	%Truncate to a random 15 units, plus the BCI units
	bcindices = find(ismember(processed.unitnames, bcunits));
	sample = 1:nUabove;
	sample(ismember(sample, bcindices)) = [];
	indices = [bcindices, randsample(sample, nU)];
	%indices = [bcindices];
	%indices = randsample(sample, nU);

	processed.binnedspikes = processed.binnedspikes(:,indices);
	processed.rates = processed.rates(:,indices);              
	processed.unitnames = processed.unitnames(indices);
	nUtotal = length(processed.unitnames);

	%Make data matrix
	if strcmp(taskaxes, 'horiz')
		idx = [1];
	elseif strcmp(taskaxes, 'vert')
		idx = [2];
	elseif strcmp(taskaxes, '2D')
		idx = [1, 2];
		display([nevfile ' is labeled as 2D; currently unsupported, skipping.'])
		return;
	else
		error([nevfile ' does not specify the task axes.'])
	end

	%MVGC
	X = [processed.cursor(:,idx), processed.rates]';
	%X = [processed.binnedspikes]';
	%Y = resample(X', 1, 2)';
	results = mvgcBCI(X, [], [], pval);

	nC = results.morder;
	causaldensity = sum(results.pwcgc_sig(1,2:end))/(nUtotal);

	%Tag with computer run on, date, last git commit
	host = hostname();
	comm = currCommit();
	stamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');

	%Get the fitID
	%fitid = getFitID(conn);
	%fitid = randi(1e9);
	unit = 'curs';
	unitnum = 1;

	previous = fetch(exec(conn, ['SELECT id FROM fits WHERE `nev file` = "' nevfile '" AND modelID = ' num2str(modelID) ' AND unit = "' unit '"']));
	if ~strcmp(previous.Data{1}, 'No Data')
		display(['Model ' num2str(modelID) ' nevfile ' nevfile ' and unit ' unit ' already analysed. Skipping'])
		continue
	end

	%Insert into fits
	tablename = 'fits';
	fitcols = {'modelID', '`nev file`', 'unit', 'unitnum', 'ncoeff', 'computer', '`analysis date`', 'commit'};
	sqldata = { modelID, nevfile, unit, unitnum, nC, host, stamp, comm};
	datainsert(conn,tablename,fitcols,sqldata);
	%Get the fit id used
	fitid = fetch(exec(conn, 'SELECT LAST_INSERT_ID()'));
	fitid = fitid.Data{1};

	%Insert into fits_mvgc
	tablename = 'fits_mvgc';
	fitcols = {'id', 'alpha', 'units', 'causaldensity'};
	sqldata = { fitid, pval, nUtotal, causaldensity};
	datainsert(conn,tablename,fitcols,sqldata);

	%For each unit, save the results 
	tablename = 'GrangerEstimates';
	fitcols = {'id', 'fromnum', 'fromunit', 'score', 'pval', 'significant'};
	for j = 2:(nUtotal+1)
		%Extract and save regression fiticients
		unit = processed.unitnames{j-1};
		%Extract deviance
		score = results.pwcgc(1,j);
		p = results.pwcgc_pval(1,j);
		sig = results.pwcgc_sig(1,j);
	
		%Insert into fits_mvgc
		sqldata = {fitid, j-1, unit, score, p, sig};
		datainsert(conn,tablename,fitcols,sqldata);
	end	
end
