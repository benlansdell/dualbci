function processLinear(conn, modelID, blackrock, nevfile, paramcode, threshold, units)
	nevpath = [blackrock nevfile];
	%Load parameters
	eval(paramcode);
	%Preprocess data
	if nargin < 7
		processed = preprocess_smooth(nevpath, binsize, sigma_fr, sigma_trq, threshold, offset);
	else
		processed = preprocess_smooth(nevpath, binsize, sigma_fr, sigma_trq, threshold, offset, 0, 0, units);
	end
	%Truncate to specified duration
	processed = truncate_recording(processed, dur);
	%Truncate to units that haven't been analyzed before using this model
	processed = removeProcessedUnits(processed, conn, modelID);
	nU = length(processed.unitnames);
	%If all units have been analyzed then return
	if nU == 0
		display(['All units in ' nevfile ' have been analyzed by model number ' modelID '. Continuing'])
		return
	end

	%Fit a linear model
	data = filters_sp_pos(processed, nK_sp, nK_pos);
	model = MLE_lmfit(data, const);
	%Tag with computer run on, date, last git commit
	host = hostname();
	comm = currCommit();
	stamp = datestr(now);
	%For each unit, save the results 
	for idx = 1:nU
		%Extract and save regression coefficients
		const = model.b_hat(idx,1);
		regrFE = model.b_hat(idx,2);
		regrRU = model.b_hat(idx,3);
		%Extract deviance
		dev = model.dev{idx, 1};

		%Extract SE

		%Extract MSE

		%Insert into Fits
		tablename = 'Fits';
		fitcols = {'modelID', '`nev file`', 'unit', 'coeffID', 'ncoeff', 'dev', 'mse', 'computer', '`analysis date`', 'commit'};
		sqldata = { modelID, nevfile, unit, coeffid, nC, dev, mse, host, stamp, comm};
		datainsert(conn,tablename,colnames,sqldata);

		%Insert into FitsLinear
		tablename = 'FitsLinear';
		fitcols = {'coeffID', 'dir', 'size'};
		sqldata = { coeffid, direction, tuningsize};
		datainsert(conn,tablename,colnames,sqldata);

		%Insert into ParameterEstimates
		tablename = 'ParameterEstimates';
		fitcols = {'coeffID', 'num', 'label', 'value', 'se'};
		for j = 1:nC

			sqldata = {coeffid, direction, tuningsize};
			datainsert(conn,tablename,colnames,sqldata);
		end
	end
end