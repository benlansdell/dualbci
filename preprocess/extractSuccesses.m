function processed = extractSuccesses(processed_orig, conn, nevfile, offset)
	%Brain control, etc, data
	processed = processed_orig;
	bs = processed.binsize;
	nB = size(processed.binnedspikes,1);

	%Get all successful trials from `trials`
	trials = fetch(exec(conn, ['SELECT `start`, `end` FROM trials WHERE `nev file` = "' nevfile '" '...
		'AND success = 1']));
	trials = cell2mat(trials.Data);

	%Get offset from `recordings`
	nevoffset = fetch(exec(conn, ['SELECT `Toffset` FROM recordings WHERE `nev file` = "' nevfile '"']));
	nevoffset = nevoffset.Data{1};

	sxidx = zeros(nB,1);
	for idx = 1:size(trials,1)
		tstart = max(1, floor((trials(idx,1)-(offset+nevoffset)/60)/bs));
		tend = min(nB, floor((trials(idx,2)-(offset+nevoffset)/60)/bs));
		sxidx(tstart:tend) = 1;
	end

	%Find success bins
	processed.binnedspikes = processed.binnedspikes(sxidx==1,:);
	processed.rates = processed.rates(sxidx==1,:);
	if isfield(processed, 'torque')
		processed.torque = processed.torque(sxidx==1,:);
		processed.dtorque = processed.dtorque(sxidx==1,:);
		processed.ddtorque = processed.ddtorque(sxidx==1,:);
	end
	if isfield(processed, 'cursor')
		processed.cursor = processed.cursor(sxidx==1,:);
		processed.dcursor = processed.dcursor(sxidx==1,:);
		processed.ddcursor = processed.ddcursor(sxidx==1,:);		
	end
end