function processMVGCmanualPaired(conn, modelID, blackrock, labviewpath, nevfile1, BCnevfile, nevfile2, paramcode)
	nevpath1 = [blackrock nevfile1];
	nevpath2 = [blackrock nevfile2];
	BCnevpath = [blackrock BCnevfile];
	nevfiles = {nevfile1, BCnevfile, nevfile2};
	%Load parameters
	eval(paramcode);
	%Get which units are BCI units, get the mat file for this nev file
	bcunits = fetch(exec(conn, ['SELECT `unit` FROM `bci_units` WHERE `ID` = "'...
	 BCnevfile '"']));
	%bcunits = num2str(cell2mat(bcunits.Data));
	bcunits = bcunits.Data;
	if size(bcunits,1) < 2
		display(['Warning: fewer than 2 BCI units are labeled within ' BCnevfile])
	end
	matfile = fetch(exec(conn, ['SELECT `labview file`,`duration` FROM `recordings`'...
		' rec WHERE rec.`nev file` = "' nevfile1 '"']));
	duration = matfile.Data{2};
	matfile = matfile.Data{1};
	matpath = [labviewpath matfile];
	taskaxes = fetch(exec(conn, ['SELECT `axis` FROM `recordings` rec WHERE rec.`nev file` = "' BCnevfile '"']));
	taskaxes = taskaxes.Data{1};
	%Get how many units are above threshold from units
	units = fetch(exec(conn, ['SELECT `unit` FROM `units` WHERE `nev file` = "'...
	 nevfile1 '" AND `firingrate` > ' num2str(threshold)]));
	units = units.Data;
	%If more than 30 units...take a random sample
	if length(units) > nU
		units = randsample(units, nU);
	end
	%Add the BCI units if not already there
	units = unique([units; bcunits]);
	if duration < dur
		display(['Recording is shorter than specified ' num2str(dur) ' seconds. Skipping'])
		return
	end

	%Check if the requisite number of units have already been analysed in this file...
	analysedunits = fetch(exec(conn, ['SELECT `unit` FROM fits WHERE modelID = '...
	 num2str(modelID) ' AND `nev file` = "' nevfile1 '"']));
	analysedunits = analysedunits.Data;
	if all(~strcmp(analysedunits, 'No Data'))
		display('Already analysed this file with this model. Continuing')
		return
	end	

	if strcmp(taskaxes, 'horiz')
		idx = [1];
	elseif strcmp(taskaxes, 'vert')
		idx = [2];
	elseif strcmp(taskaxes, '2D')
		idx = [1, 2];
	else
		display([BCnevfile ' does not specify the task axes.'])
		return;
	end
	nG = length(idx);

	%%%%%%
	%MC I%
	%%%%%%

	processed = preprocess_spline_lv(nevpath1, matfile, binsize, threshold, offset, [], [], units);
	processed = truncate_recording(processed, dur);
	nUtotal = length(processed.unitnames);
	X = [processed.cursor(:,idx), processed.rates]';
	results{1} = mvgcBCI(X, [], [], pval);
	nCs{1} = results{1}.morder;
	causaldensities{1} = sum(results{1}.pwcgc_sig(1,2:end))/(nUtotal);

	%%%%
	%BC%
	%%%%

	processed = preprocess_spline_lv(BCnevpath, matfile, binsize, threshold, offset, [], [], units);
	processed = truncate_recording(processed, dur);
	X = [processed.cursor(:,idx), processed.rates]';
	results{2} = mvgcBCI(X, [], [], pval);
	nCs{2} = results{2}.morder;
	causaldensities{2} = sum(results{2}.pwcgc_sig(1,2:end))/(nUtotal);

	%%%%%%%
	%MC II%
	%%%%%%%

	processed = preprocess_spline_lv(nevpath2, matfile, binsize, threshold, offset, [], [], units);
	processed = truncate_recording(processed, dur);
	X = [processed.cursor(:,idx), processed.rates]';
	results{3} = mvgcBCI(X, [], [], pval);
	nCs{3} = results{3}.morder;
	causaldensities{3} = sum(results{3}.pwcgc_sig(1,2:end))/(nUtotal);

	%%%%%%%%%%%%%%%%%%%
	%Database addition%
	%%%%%%%%%%%%%%%%%%%

	%Tag with computer run on, date, last git commit
	host = hostname()
	stamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
	comm = currCommit();
	unit = 'curs';
	unitnum = 1;

	for i = 1:length(results)
		result = results{i};
		nC = nCs{i};
		causaldensity = causaldensities{i};
		nevfile = nevfiles{i};

		previous = fetch(exec(conn, ['SELECT id FROM fits WHERE `nev file` = "'...
			nevfile '" AND modelID = ' num2str(modelID) ' AND unit = "' unit '"']));
		if ~strcmp(previous.Data{1}, 'No Data')
			display(['Model ' num2str(modelID) ' nevfile ' nevfile ' and unit '... 
				unit ' already analysed. Skipping'])
			continue
		end
	
		%Insert into fits
		tablename = 'fits';
		fitcols = {'modelID', '`nev file`', 'unit', 'unitnum', 'ncoeff', 'computer',...
		 '`analysis date`', 'commit'};
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
		tablename = 'estimates_granger';
		fitcols = {'id', 'fromnum', 'fromunit', 'score', 'pval', 'significant'};
		for j = (nG+1):(nUtotal+nG)
			%Extract and save regression fiticients
			fromunit = processed.unitnames{j-nG};
			%Extract deviance
			score = result.pwcgc(1,j);
			p = result.pwcgc_pval(1,j);
			sig = result.pwcgc_sig(1,j);
			%Insert into fits_mvgc
			sqldata = {fitid, j-nG, fromunit, score, p, sig};
			datainsert(conn,tablename,fitcols,sqldata);
		end	
	end
end
