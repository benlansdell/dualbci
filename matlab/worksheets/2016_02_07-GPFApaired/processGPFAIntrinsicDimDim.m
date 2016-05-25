function [analysis_id, ll, D] = processGPFAIntrinsicDim(conn, modelID, blackrock, labviewpath, labviewfile, MCnevfile1, BCnevfile1, DCnevfile, expt_id, paramcode)
	%Keep the units of the network the same between recordings
	%Load parameters
	eval(paramcode);
	%threshold = 5;

	%Figure out units to use... above 5hz in all recordings
	units = fetch(exec(conn, ['SELECT u1.unit FROM '...
	'`experiment_tuning` et1 '...
	'INNER JOIN `units` u1 '...
	'ON u1.`nev file` = et1.`manualrecording` '...
	'INNER JOIN `units` u2 '...
	'ON u2.`nev file` = et1.`1DBCrecording` '...
	'INNER JOIN `units` u4 '...
	'ON u4.`nev file` = et1.`dualrecording` '...
	'WHERE u1.unit = u2.unit AND u1.unit = u4.unit AND '...
	'u1.firingrate > ' num2str(threshold) ' AND u2.firingrate > ' num2str(threshold) ' AND '...
	'u4.firingrate > ' num2str(threshold) ' AND et1.`manualrecording` = "' MCnevfile1 '"']));
	otherunits = units.Data;
 
	%Make sure BC units from both dual and brain control are in there
	bciunits = exec(conn, ['SELECT `unit` FROM bci_units WHERE `ID` = "' BCnevfile1 '"']);
	bciunits = fetch(bciunits);
	bciunits = bciunits.Data;
	dualunits = exec(conn, ['SELECT `unit` FROM bci_units WHERE `ID` = "' DCnevfile '"']);
	dualunits = fetch(dualunits);
	dualunits = dualunits.Data;
	allunits = unique([otherunits; bciunits; dualunits]);

	%Tag with computer run on, date, last git commit
	host = hostname();
	stamp = datestr(now, 'yyyy-mm-dd HH:MM:SS');
	comm = currCommit();

	previous = fetch(exec(conn, ['SELECT id FROM analyses WHERE `experiment_id` = "' num2str(expt_id) '" AND modelID = ' num2str(modelID)]));
	if ~strcmp(previous.Data{1}, 'No Data')
		display(['modelID ' num2str(modelID) ' and experiment_id ' num2str(expt_id) ' already analysed.'])
		analysis_id = previous.Data{1};
		[llMC, D] = runGPFAIntPaired(conn, analysis_id, 1, modelID, blackrock, labviewpath, labviewfile, MCnevfile1, allunits, expt_id, paramcode);
		ll = {llMC};
	end
end

function [ll, D] = runGPFAIntPaired(conn, analysis_id, id, modelID, blackrock, labviewpath, labviewfile, nevfile, units, expt_id, paramcode)
	nevpath = [blackrock nevfile];
	lbpath = [labviewpath labviewfile];
	%Load parameters
	eval(paramcode);
	nU = length(units);

	processed = preprocess_spline_target(nevpath, lbpath, binsize, threshold, offset, [], [], units);
	processed = truncate_recording(processed, dur);
	%Time bins which occur within a trial
	spikes = processed.binnedspikes;
	trials = sum(abs(processed.target),2)>0;

	%Run all trials together
	[dat, octs, quads] = makedata(trials, spikes, processed.cursor, processed.target);
	method = 'gpfa';
	% Select number of latent dimensions
	kernSD = 30;
	% Extract neural trajectories
	runIdx = analysis_id*10 + id;
	nd = 2:nU;
	DDs = zeros(nfolds, nU, nU);
	for i = nd
		result = neuralTraj(runIdx, dat, 'method', method, 'xDim', i, 'kernSDList', kernSD, 'numFolds', nfolds);
		for idx = 1:nfolds
			%Load data
			if i > 9
				fn = ['./mat_results/run' num2str(runIdx) '/gpfa_xDim' num2str(i) '_cv0' num2str(idx) '.mat'];
			else
				fn = ['./mat_results/run' num2str(runIdx) '/gpfa_xDim0' num2str(i) '_cv0' num2str(idx) '.mat'];
			end
			result = load(fn);
			% Orthonormalize neural trajectories
			[estParams, seqTrain, seqTest, DD] = postprocess(result, 'kernSD', kernSD);
			if i == nU
				DDs(idx,:,:) = DD;
			end
			%Average log likelihood for test set
			Tlen = 0;
			for j = 1:length(seqTest)
				Tlen = Tlen + seqTest(j).T;
			end
			[s, LL] = exactInferenceWithLL(seqTest, estParams, 'getLL', 1);
			ll(i,idx) = LL/Tlen;
		end
	end
	D = squeeze(mean(DDs, 1));
	%Save results
end 
