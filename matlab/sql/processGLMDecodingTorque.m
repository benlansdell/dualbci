function processGLMDecodingTorque(conn, modelID, blackrock, labview, nevfile, BCnevfile, nevfile2, paramcode, verbose)
	if nargin < 9
		verbose = 0;
	end
	N_sub = 5000;
	nevpath = [blackrock nevfile];
	BCnevpath = [blackrock BCnevfile];
	nevpath2 = [blackrock nevfile2];
	%Load parameters
	eval(paramcode);
	%Get which units are BCI units, get the mat file for this nev file
	bcunits = fetch(exec(conn, ['SELECT `unit` FROM `BCIUnits` WHERE `ID` = "' BCnevfile '"']));
	%bcunits = num2str(cell2mat(bcunits.Data));
	bcunits = bcunits.Data;
	if size(bcunits,1) < 2
		display(['Warning: fewer than 2 BCI units are labeled within ' nevfile])
	end
	matfile = fetch(exec(conn, ['SELECT `labview file`,`duration` FROM `Recordings` rec WHERE rec.`nev file` = "' nevfile '"']));
	duration = matfile.Data{2};
	matfile = matfile.Data{1};
	matpath = [labview matfile];
	taskaxes = fetch(exec(conn, ['SELECT `axis` FROM `Recordings` rec WHERE rec.`nev file` = "' BCnevfile '"']));
	taskaxes = taskaxes.Data{1};
	if strcmp(taskaxes, 'horiz')
		mask = [1];
	elseif strcmp(taskaxes, 'vert')
		mask = [2];
	elseif strcmp(taskaxes, '2D')
		mask = [1, 2];
	else
		display([nevfile ' does not specify the task axes.'])
		return;
	end
	%%%
	%I%
	%%%

	%First part: 15 units that are not BC, decode MC, then using the same decoder, decode BC
	%For reference, train those same units on BC, decode on BC
	%For reference, train those same units plus the BC units on BC task, decode on BC

	%Get top 15 units that are not in the BC task
	bcunitsstr = ['unit != ', strjoin(bcunits, ' AND unit != ')];
	searchstr = ['SELECT unit FROM `Units` WHERE `nev file` = "' nevfile '"'...
	' AND ' bcunitsstr ' ORDER BY `Units`.`firingrate` DESC LIMIT ' num2str(nU)];
	units = fetch(exec(conn, searchstr));
	units = units.Data;

	%Preprocess
	processed = preprocess_spline_lv(nevpath, matfile, binsize, threshold, offset, [], [], units);
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
		fn_out = ['./worksheets/2015_07_27/' nevfile '_GLMDecodingMCMC.eps'];
	else
		fn_out = '';
	end
	decoded_cursor_MCMC = glm_decode(processed_novel, data, model, F, Q, mu, fn_out);
	%Take a random sample... so it doesn't take forever...
	ii = 1:size(decoded_cursor_MCMC,1);
	ii_sub = datasample(ii, nS, 'Replace', false);
	corrMCMCcurs(1) = corr(decoded_cursor_MCMC(ii_sub,1), data.cursor(ii_sub,1));
	corrMCMCcurs(2) = corr(decoded_cursor_MCMC(ii_sub,2), data.cursor(ii_sub,2));
	corrMCMCtorq(1) = corr(decoded_cursor_MCMC(ii_sub,1), data.torque(ii_sub,1));
	corrMCMCtorq(2) = corr(decoded_cursor_MCMC(ii_sub,2), data.torque(ii_sub,2));
	%corrMCMC
	%Only use the axis from the BC task in evaluation

	%Decode on BCI data
	%Preprocess
	processed = preprocess_spline_lv(BCnevpath, matfile, binsize, threshold, offset, [], [], units);
	%Truncate to specified duration
	[processed, processed_novel] = split_recording(processed, dur, dur+testdur);
	data = filters_sprc_pos(processed_novel, nK_sp, nK_pos, dt_sp, dt_pos);
	nS = min(N_sub, size(data.X, 1));
	if verbose > 0
		fn_out = ['./worksheets/2015_07_27/' nevfile '_GLMDecodingMCBC.eps'];
	else
		fn_out = '';
	end
	[F, Q, mu] = fit_AR_LS(data.torque(:,mask), order);
	masked = setDiff(mask, [1,2]);
	modelmasked = model;
	modelmasked.mask(:,2+masked) = zeros(nU,1);
	decoded_cursor_MCBC = glm_decode(processed_novel, data, modelmasked, F, Q, mu, fn_out);
	%Eval performance
	%Only use the axis from the BC task in evaluation
	%Take a random sample... so it doesn't take forever...
	ii = 1:size(decoded_cursor_MCBC,1);
	ii_sub = datasample(ii, nS, 'Replace', false);
	corrMCBCcurs = zeros(2,1);
	corrMCBCtorq = zeros(2,1);
	for midx = mask
		corrMCBCcurs(midx) = corr(decoded_cursor_MCBC(ii_sub,midx), data.cursor(ii_sub,midx));
		corrMCBCtorq(midx) = corr(decoded_cursor_MCBC(ii_sub,midx), data.torque(ii_sub,midx));
	end
	%Eval performance

	%Decode on second MC recording
	%Preprocess
	if ~strcmp(nevfile2, 'null')
		processed = preprocess_spline_lv(nevpath2, matfile, binsize, threshold, offset, [], [], units);
		%Truncate to specified duration
		[processed, processed_novel] = split_recording(processed, dur, dur+testdur);
		data = filters_sprc_pos(processed_novel, nK_sp, nK_pos, dt_sp, dt_pos);
		[F, Q, mu] = fit_AR_LS(data.torque, order);
		nS = min(N_sub, size(data.X, 1));
		if verbose > 0
			fn_out = ['./worksheets/2015_07_27/' nevfile '_GLMDecodingMCMC2.eps'];
		else
			fn_out = '';
		end
		decoded_cursor_MCMC2 = glm_decode(processed_novel, data, model, F, Q, mu, fn_out);
		%Eval performance
		%Only use the axis from the BC task in evaluation
		%Take a random sample... so it doesn't take forever...
		ii = 1:size(decoded_cursor_MCBC,1);
		ii_sub = datasample(ii, nS, 'Replace', false);
		corrMCMC2curs(1) = corr(decoded_cursor_MCMC2(ii_sub,1), data.cursor(ii_sub,1));
		corrMCMC2curs(2) = corr(decoded_cursor_MCMC2(ii_sub,2), data.cursor(ii_sub,2));
		corrMCMC2torq(1) = corr(decoded_cursor_MCMC2(ii_sub,1), data.torque(ii_sub,1));
		corrMCMC2torq(2) = corr(decoded_cursor_MCMC2(ii_sub,2), data.torque(ii_sub,2));
	else
		corrMCMC2curs = zeros(2,1);
		corrMCMC2torq = zeros(2,1);
	end

	%Train the same units on BC task
	processed = preprocess_spline_lv(BCnevpath, matfile, binsize, threshold, offset, [], [], units);
	%Truncate to specified duration
	[processed, processed_novel] = split_recording(processed, dur, dur+testdur);
	data = filters_sprc_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
	model = MLE_glmfit(data, const);
	%Only fit AR to cursor components used in task
	[F, Q, mu] = fit_AR_LS(data.torque(:,mask), order);
	%Decode
	data = filters_sprc_pos(processed_novel, nK_sp, nK_pos, dt_sp, dt_pos);
	nS = min(N_sub, size(data.X, 1));	
	if verbose > 0
		fn_out = ['./worksheets/2015_07_27/' nevfile '_GLMDecodingBCBC.eps'];
	else
		fn_out = '';
	end
	decoded_cursor_BCBC = glm_decode(processed_novel, data, model, F, Q, mu, fn_out);
	%Take a random sample... so it doesn't take forever...
	ii = 1:size(decoded_cursor_BCBC,1);
	ii_sub = datasample(ii, nS, 'Replace', false);
	corrBCBCcurs = zeros(2,1);
	corrBCBCtorq = zeros(2,1);
	for midx = mask
		corrBCBCcurs(midx) = corr(decoded_cursor_BCBC(ii_sub,midx), data.cursor(ii_sub,midx));
		corrBCBCtorq(midx) = corr(decoded_cursor_BCBC(ii_sub,midx), data.torque(ii_sub,midx));
	end

	%Train the same units on BC task, but using cursor velocity
	%data = filters_sprc_vel_lv(processed, nK_sp, nK_pos, dt_sp, dt_pos);
	%model = MLE_glmfit(data, const);
	%%Only fit AR to cursor components used in task
	%[F, Q, mu] = fit_AR_LS(data.dcursor(:,mask), order);
	%%Decode
	%data = filters_sprc_vel_lv(processed_novel, nK_sp, nK_pos, dt_sp, dt_pos);
	%nS = min(N_sub, size(data.X, 1));	
	%if verbose > 0
	%	fn_out = './testGLMDecodingBCBCvel.eps';
	%else
	%	fn_out = '';
	%end
	%decoded_cursor_BCBCvel = glm_decode(processed_novel, data, model, F, Q, mu, fn_out);
	%%Take a random sample... so it doesn't take forever...
	%ii = 1:size(decoded_cursor_BCBC,1);
	%ii_sub = datasample(ii, nS, 'Replace', false);
	%corrBCBCvelcurs = zeros(2,1);
	%corrBCBCveltorq = zeros(2,1);
	%for midx = mask
	%	corrBCBCvelcurs(midx) = corr(decoded_cursor_BCBCvel(ii_sub,midx), data.dcursor(ii_sub,midx));
	%	corrBCBCveltorq(midx) = corr(decoded_cursor_BCBCvel(ii_sub,midx), data.dtorque(ii_sub,midx));
	%end

	%Do the same, but include the BC units themselves... as another reference score
	%units = unique([units; bcunits]);
	%processed = preprocess_spline_lv(BCnevpath, matfile, binsize, threshold, offset, [], [], units);
	%%Truncate to specified duration
	%[processed, processed_novel] = split_recording(processed, dur, dur+testdur);
	%data = filters_sprc_vel_lv(processed, nK_sp, nK_pos, dt_sp, dt_pos);
	%model = MLE_glmfit(data, const);
	%data = filters_sprc_pos(processed_novel, nK_sp, nK_pos, dt_sp, dt_pos);
	%[F, Q, mu] = fit_AR_LS(data.dcursor(:,mask), order);
	%decoded_cursor_BCBCBC = glm_decode(processed_novel, data, model, F, Q, mu);
	%%Eval performance
	%%Only use the axis from the BC task in evaluation
	%%Take a random sample... so it doesn't take forever...
	%ii = 1:size(decoded_cursor_BCBCBC,1);
	%ii_sub = datasample(ii, nS, 'Replace', false);
	%corrXBCBCBC = corr(decoded_cursor_BCBCBC(ii_sub,1), data.dcursor(ii_sub,1))
	%corrYBCBCBC = corr(decoded_cursor_BCBCBC(ii_sub,2), data.dcursor(ii_sub,2))	%Eval performance

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
		%Insert into FitsGLMDecode
		tablename = 'FitsGLMDecode';
		fitcols = {'id', 'corrCursMCMC', 'corrCursMCBC', 'corrCursBCBC', 'corrCursMCMC2', 'corrTorqMCMC', 'corrTorqMCBC', 'corrTorqBCBC', 'corrTorqMCMC2'};
		sqldata = { fitid, corrMCMCcurs(idx), corrMCBCcurs(idx), corrBCBCcurs(idx), corrMCMC2curs(idx), corrMCMCtorq(idx), corrMCBCtorq(idx), corrBCBCtorq(idx), corrMCMC2torq(idx)};
		datainsert(conn,tablename,fitcols,sqldata);	
	end
end
