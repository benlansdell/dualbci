function processGrangerGLM(conn, modelID, blackrock, nevfile, paramcode, threshold)
	nevpath = [blackrock nevfile];
	%Load parameters
	eval(paramcode);
	%Preprocess data
	processed = preprocess_spline(nevpath, binsize, threshold, offset);
	nU = length(processed.unitnames);
	%Truncate to specified duration
	processed = truncate_recording(processed, dur);
	%Truncate to top 15 units by firing rate
	[totalspikes, indices] = sort(sum(processed.binnedspikes,1), 2, 'descend');
	nU = 15; %size(processed.binnedspikes,2);
	nUabove = length(processed.unitnames);
	if nUabove < nU
		display(['Warning: nev file ' nevfile ' contains only ' num2str(nUabove)...
		 ' units above ' num2str(threshold) ' when ' num2str(nU) ' was requested.'])
		nU = nUabove;
	end
	top = indices(1:nU);
	processed.binnedspikes = processed.binnedspikes(:,top);
	processed.rates = processed.rates(:,top);              
	processed.unitnames = processed.unitnames(top);
	%Run with position filters
	gdata = filters_sp_pos_network(processed, nK_sp, nK_pos);
	[GCdev, GCpval, GCsig, fulldevs] = granger(processed, gdata, [], pval);
	causaldensity = sum(sum(GCsig))/nU/(nU-1);

	nC = size(gdata.X,2);
	if strcmp(const, 'on')
		nC = nC + 1;
	end
	%Tag with computer run on, date, last git commit
	host = hostname();
	comm = currCommit();
	stamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
	%For each unit, save the results 
	for idx = 1:nU
		%Extract and save regression fiticients
		unit = processed.unitnames{idx};
		%Extract deviance
		dev = fulldevs(idx);

		%Insert into Fits
		tablename = 'Fits';
		fitcols = {'modelID', '`nev file`', 'unit', 'unitnum', 'ncoeff', 'dev', 'computer', '`analysis date`', 'commit'};
		sqldata = { modelID, nevfile, unit, idx, nC, dev, host, stamp, comm};
		%sqldata = { 1, '20130920SpankyUtah001.nev', 999, 1, 3, 3, '3', '2013-12-09 12:12:12', '12'};
		datainsert(conn,tablename,fitcols,sqldata);
		%Get the fit id used
		fitid = fetch(exec(conn, 'SELECT LAST_INSERT_ID()'));
		fitid = fitid.Data{1};

		%Insert into FitsGranger
		tablename = 'FitsGranger';
		fitcols = {'id', 'alpha', 'units', 'causaldensity'};
		sqldata = { fitid, pval, nU, causaldensity};
		datainsert(conn,tablename,fitcols,sqldata);

		%Insert into NetworkEstimates
		tablename = 'GrangerEstimates';
		fitcols = {'id', 'fromnum', 'fromunit', 'score', 'pval', 'significant'};
		for j = 1:nU
			if j ~= idx
				unitj = processed.unitnames{j};
				gc = GCdev(idx, j);		
				gcpval = GCpval(idx, j);
				gcsig = GCsig(idx, j);
				sqldata = { fitid, j, unitj, gc, gcpval, gcsig};
				datainsert(conn,tablename,fitcols,sqldata);
			end
		end
	end	
end
