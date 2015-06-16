function processed = removeProcessedUnits(processed, conn, modelID)
	unitnames = processed.unitnames;
	%Remove units that have already been processed
	[pathstr, nevfile, ext] = fileparts(processed.nevfile);
	nevfile = [nevfile '.nev'];
	%Find which of these units have been analyzed with this modelID and nevfile
	processedunits = exec(conn, ['SELECT `unit` FROM Fits WHERE modelID = ' num2str(modelID) ' AND `nev file` = "' nevfile '"']);
	processedunits = fetch(processedunits);
	processedunits = processedunits.Data;
	if strcmp(processedunits, 'No Data')
		processedunits = {};
	end
	for idx = 1:length(processedunits)
		processedunits{idx} = num2str(processedunits{idx});
	end
	unprocessed = not(ismember(unitnames, processedunits));
	%Truncate
	processed.binnedspikes = processed.binnedspikes(:,unprocessed);
	processed.rates = processed.rates(:,unprocessed);
	processed.unitnames = processed.unitnames(:,unprocessed);
end