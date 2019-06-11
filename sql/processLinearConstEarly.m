function processLinearConstEarly(conn, modelID, blackrock, labviewpath, nevfile, paramcode, units)
	nevpath = [blackrock nevfile];
	%Load parameters
	eval(paramcode);
	%If too short then skip this file...
	duration = fetch(exec(conn, ['SELECT `duration` FROM recordings WHERE `nev file` = "' nevfile '"']));
	duration = duration.Data{1};
	if duration < dur+testdur
		display(['Duration of ' nevfile '(' num2str(duration) 's) is less than requested ' num2str(dur+testdur) 's. Continuing'])		
		return
	end
	%Preprocess data
	if nargin < 7
		processed = preprocess_spikes(nevpath, binsize, threshold, offset);
		units = 0;
	else
		processed = preprocess_spikes(nevpath, binsize, threshold, offset, 0, 0, units);
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
	nBout = size(processed_novel.binnedspikes,1);
	nC = 1;
	%Tag with computer run on, date, last git commit
	host = hostname();
	stamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
	comm = currCommit();
	%For each unit, save the results 
	for idx = 1:nU
		unit = processed.unitnames{idx};
		previous = fetch(exec(conn, ['SELECT id FROM fits WHERE `nev file` = "' nevfile '" AND modelID = ' num2str(modelID) ' AND unit = "' unit '"']));
		if ~strcmp(previous.Data{1}, 'No Data')
			display(['Model ' num2str(modelID) ' nevfile ' nevfile ' and unit ' unit ' already analysed. Skipping'])
			continue
		end
		%Extract MSE. Need to add cross-validation code to do this... add later
		b_hat = mean(processed.binnedspikes(:,idx));
		converged = 1;
		conditioned = 1;
		mseout = sum((processed_novel.binnedspikes(:,idx)-b_hat).^2)/nBout;
		%Insert into fits
		tablename = 'fits';
		fitcols = {'modelID', '`nev file`', 'unit', 'ncoeff', '`mse out`', 'computer', '`analysis date`', 'commit', 'conv', 'cond'};
		sqldata = { modelID, nevfile, unit, nC, mseout, host, stamp, comm, converged, conditioned};
		datainsert(conn,tablename,fitcols,sqldata);
	end
end