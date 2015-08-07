function processLinear(conn, modelID, blackrock, labview, nevfile, matfile, paramcode, threshold, units)
	nevpath = [blackrock nevfile];
	matpath = [labview matfile];
	%Load parameters
	eval(paramcode);
	%Preprocess data
	if nargin < 9
		processed = preprocess_spline_lv(nevpath, matpath, binsize, threshold, offset);
	else
		processed = preprocess_spline_lv(nevpath, matpath, binsize, threshold, offset, 0, 0, units);
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
	data = filters_sp_pos_lv(processed, nK_sp, nK_pos);
	model = MLE_lmfit(data, const);
	%%Compute MSE on test data
	datanovel = filters_sp_pos_lv(processed_novel, nK_sp, nK_pos);
	nBout = size(datanovel.y,2);
	nBin = size(data.y,2);

	nC = size(model.b_hat,2);
	%Tag with computer run on, date, last git commit
	host = hostname();
	stamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
	comm = currCommit();
	%For each unit, save the results 
	for idx = 1:nU
		%If already in database, skip
		%unit = '21.3'; modelID = 2; nevfile = '20140610SpankyUtah002.nev';
		%Extract and save regression fiticients
		unit = processed.unitnames{idx};
		previous = fetch(exec(conn, ['SELECT id FROM fits WHERE `nev file` = "' nevfile '" AND modelID = ' num2str(modelID) ' AND unit = "' unit '"']));
		if ~strcmp(previous.Data{1}, 'No Data')
			display(['Model ' num2str(modelID) ' nevfile ' nevfile ' and unit ' unit ' already analysed. Skipping'])
			continue
		end
		%Extract deviance
		dev = model.dev{idx, 1};
		%Extract SE
		stats = model.stats{idx};
		%Extract MSE. Need to add cross-validation code to do this... add later
		b_hat = model.b_hat(idx,:);
		rho_hat = glmval(b_hat', squeeze(datanovel.X(idx,:,:)), 'identity');

		mseout = sum((datanovel.y(idx,:)-rho_hat').^2)/nBout;

		%Get the fitID
		%fitid = getFitID(conn);
		%fitid = randi(1e9);
		%Insert into fits
		tablename = 'fits';
		fitcols = {'modelID', '`nev file`', 'unit', 'ncoeff', 'dev', '`mse out`', 'computer', '`analysis date`', 'commit'};
		sqldata = { modelID, nevfile, unit, nC, dev, mseout, host, stamp, comm};
		datainsert(conn,tablename,fitcols,sqldata);

		%Get the fit id used
		fitid = fetch(exec(conn, 'SELECT LAST_INSERT_ID()'));
		fitid = fitid.Data{1};

		%Insert into fits_linear
		tablename = 'fits_linear';
		fitcols = {'id', 'dir', 'size'};
		regrRU = model.b_hat(idx,2);
		regrFE = model.b_hat(idx,3);
		[direction, tuningsize] = unitTheta(regrRU, regrFE);
		sqldata = { fitid, direction, tuningsize};
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
