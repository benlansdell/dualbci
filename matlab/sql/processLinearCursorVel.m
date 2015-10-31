function processLinearCursorVel(conn, modelID, blackrock, labview, nevfile, matfile, paramcode, threshold, units, rerun)
	logfile = './sqllog.txt';
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
	if nargin < 10
		rerun = 0;
	end
	%Truncate to units that haven't been analyzed before using this model
	if rerun == 0
		processed = removeProcessedUnits(processed, conn, modelID);
	end
	nU = length(processed.unitnames);
	%Truncate to specified duration
	[processed, processed_novel] = split_recording(processed, dur, dur+testdur);
	%If all units have been analyzed then return
	if nU == 0
		display(['All units in ' nevfile ' have been analyzed by model number ' modelID '. Continuing'])
		return
	end

	%Fit a linear model to training data
	d = ['processLinearCursor: Fitting model ' num2str(modelID) ' nevfile ' nevfile];
	display(d);
	writelog(logfile, d);
	data = filters_sp_vel_lv(processed, nK_sp, nK_pos);
	model = MLE_lmfit(data, const);
	%%Compute MSE on test data
	datanovel = filters_sp_vel_lv(processed_novel, nK_sp, nK_pos);
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
		previous = fetch(exec(conn, ['SELECT id FROM fits WHERE `nev file` = "' nevfile '" AND modelID = ' num2str(modelID) ' AND unit = "' unit '"'...
			' AND analyses_id IS NULL']));
		if ~strcmp(previous.Data{1}, 'No Data')
			if ~rerun 
				d = ['processLinearCursor: WARNING Model ' num2str(modelID) ' nevfile ' nevfile ' and unit ' unit ' already analysed. Skipping'];
				display(d);
				writelog(logfile, d);
				continue
			elseif length(previous.Data) > 1
				d = ['processLinearCursor: WARNING Model ' num2str(modelID) ' nevfile ' nevfile ' and unit ' unit ' already analysed. Trying to overwrite '...
					' old fit, but found too many matching fits for this model to unambiguously replace with new data. Skipping'];
				display(d);
				writelog(logfile, d);
				continue				
			else
				d = ['processLinearCursor: rerun = 1: Model ' num2str(modelID) ' nevfile ' nevfile ' and unit ' unit ' already analysed. Deleting and '...
				're-entering.'];
				display(d);
				writelog(logfile, d);
				exec(conn, ['DELETE FROM fits WHERE `nev file` = "' nevfile '" AND modelID = ' num2str(modelID) ' AND unit = "' unit '"'...
					' AND analyses_id IS NULL']);
			end
		end

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
		fitcols = {'id', 'dir', 'size', 'r2', 'r2out'};
		regrRU = model.b_hat(idx,2);
		regrFE = model.b_hat(idx,3);
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
