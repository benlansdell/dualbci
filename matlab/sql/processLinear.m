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
	%Truncate to units that haven't been analyzed before using this model
	processed = removeProcessedUnits(processed, conn, modelID);
	nU = length(processed.unitnames);
	%Truncate to specified duration
	[processed, processed_novel] = split_recording(processed, dur, dur+testdur);
	%If all units have been analyzed then return
	if nU == 0
		display(['All units in ' nevfile ' have been analyzed by model number ' modelID '. Continuing'])
		return
	end

	%Fit a linear model to training data
	data = filters_sp_pos(processed, nK_sp, nK_pos);
	model = MLE_lmfit(data, const);
	%%Compute MSE on test data
	datanovel = filters_sp_pos(processed_novel, nK_sp, nK_pos);

	nC = size(model.b_hat,2);
	%Tag with computer run on, date, last git commit
	host = hostname();
	comm = currCommit();
	stamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
	%For each unit, save the results 
	for idx = 1:nU
		%Extract and save regression fiticients
		unit = processed.unitnames{idx};
		%Extract deviance
		dev = model.dev{idx, 1};
		%Extract SE
		stats = model.stats{idx};
		%Extract MSE. Need to add cross-validation code to do this... add later
		b_hat = model.b_hat(idx,:);
		rho_hat = glmval(b_hat', squeeze(datanovel.X(idx,:,:)), 'identity');
		mseout = sum((datanovel.y(idx,:)-rho_hat').^2);

		%Get the fitID
		fitid = getFitID(conn);
		%Insert into Fits
		tablename = 'Fits';
		fitcols = {'modelID', '`nev file`', 'unit', 'fitID', 'nCoeff', 'dev', 'mse', 'computer', '`analysis date`', 'commit'};
		sqldata = { modelID, nevfile, unit, fitid, nC, dev, mseout, host, stamp, comm};
		datainsert(conn,tablename,fitcols,sqldata);

		%Insert into FitsLinear
		tablename = 'FitsLinear';
		fitcols = {'fitID', 'dir', 'size'};
		regrRU = model.b_hat(idx,2);
		regrFE = model.b_hat(idx,3);
		[direction, tuningsize] = unitTheta(regrRU, regrFE);
		sqldata = { fitid, direction, tuningsize};
		datainsert(conn,tablename,colnames,sqldata);

		%Insert into ParameterEstimates
		tablename = 'ParameterEstimates';
		fitcols = {'fitID', 'num', 'label', 'value', 'se'};
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
			sqldata = {fitid, num, label, model.b_hat(idx, j), stats.se(j)};
			datainsert(conn,tablename,colnames,sqldata);
		end
	end
end