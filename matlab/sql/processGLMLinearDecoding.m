function processGLMLinearDecoding(conn, modelID, blackrock, nevfile, paramcode, verbose)
	if nargin < 9
		verbose = 0;
	end
	N_sub = 5000;
	nevpath = [blackrock nevfile];
	%Load parameters
	eval(paramcode);
	%%%
	%I%
	%%%

	%Get top 15 units that are not in the BC task
	searchstr = ['SELECT unit FROM `Units` WHERE `nev file` = "' nevfile '"'...
	' ORDER BY `Units`.`firingrate` DESC LIMIT ' num2str(nU)];
	units = fetch(exec(conn, searchstr));
	units = units.Data;

	%Preprocess
	processed = preprocess_spline(nevpath, binsize, threshold, offset, [], [], units);
	%Truncate to specified duration
	[processed, processed_novel] = split_recording(processed, dur, dur+testdur);
	nUtotal = length(processed.unitnames);

	%Prepare data
	data = filters_sprc_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
	nS = min(N_sub, size(data.X, 1));
	%Train
	model = MLE_glmfit(data, const);
	[F, Q, mu] = fit_AR_LS(data.torque, order);
	%Decode
	data = filters_sprc_pos(processed_novel, nK_sp, nK_pos, dt_sp, dt_pos);
	if verbose > 0
		fn_out = ['./worksheets/2015_07_28-LinearVsGLMDecoding/' nevfile '_GLMDecodingMCMC.eps'];
	else
		fn_out = '';
	end
	decoded_cursor_MCMC = glm_decode(processed_novel, data, model, F, Q, mu, fn_out);
	%Take a random sample... so it doesn't take forever...
	ii = 1:size(decoded_cursor_MCMC,1);
	ii_sub = datasample(ii, nS, 'Replace', false);
	corrTorqGLM(1) = corr(decoded_cursor_MCMC(ii_sub,1), data.torque(ii_sub,1));
	corrTorqGLM(2) = corr(decoded_cursor_MCMC(ii_sub,2), data.torque(ii_sub,2));


	%Train linear model
	%Train
	data = filters_sprc_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
	model = MLE_lmfit(data, const);
	%Decode
	data = filters_sprc_pos(processed_novel, nK_sp, nK_pos, dt_sp, dt_pos);
	if verbose > 0
		fn_out = ['./worksheets/2015_07_28-LinearVsGLMDecoding/' nevfile '_LMDecodingMCMC.eps'];
	else
		fn_out = '';
	end
	decoded_cursor_MCMC_lin = lm_decode(processed_novel, data, model, F, Q, mu, fn_out);
	%Take a random sample... so it doesn't take forever...
	ii = 1:size(decoded_cursor_MCMC,1);
	ii_sub = datasample(ii, nS, 'Replace', false);
	corrTorqLM(1) = corr(decoded_cursor_MCMC_lin(ii_sub,1), data.torque(ii_sub,1));
	corrTorqLM(2) = corr(decoded_cursor_MCMC_lin(ii_sub,2), data.torque(ii_sub,2));

	%%%%%%%%%%%%%%%%%%%%%%
	%Insert into database%
	%%%%%%%%%%%%%%%%%%%%%%
	%Tag with computer run on, date, last git commit
	host = hostname();
	stamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
	cursnames = {'cursX', 'cursY'};
	comm = currCommit();
	for idx = 1:2
		unit = cursnames{idx};
		unitnum = idx;
		previous = fetch(exec(conn, ['SELECT id FROM Fits WHERE `nev file` = "' nevfile '" AND modelID = ' num2str(modelID) ' AND unit = "' unit '"']));
		if ~strcmp(previous.Data{1}, 'No Data')
			display(['Model ' num2str(modelID) ' nevfile ' nevfile ' and unit ' unit ' already analysed. Skipping'])
			continue
		end
		%Insert into Fits
		tablename = 'Fits';
		fitcols = {'modelID', '`nev file`', 'unit', 'unitnum', 'ncoeff', 'computer', '`analysis date`', 'commit'};
		sqldata = { modelID, nevfile, unit, unitnum, nU, host, stamp, comm};
		datainsert(conn,tablename,fitcols,sqldata);
		%Get the fit id used
		fitid = fetch(exec(conn, 'SELECT LAST_INSERT_ID()'));
		fitid = fitid.Data{1};
		%Insert into FitsGLMVsLMDecode
		tablename = 'FitsGLMVsLMDecode';
		fitcols = {'id', 'corrTorqGLM', 'corrTorqLM'};
		sqldata = { fitid, corrTorqGLM(idx), corrTorqLM(idx)};
		datainsert(conn,tablename,fitcols,sqldata);	
	end
end
