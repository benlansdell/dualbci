%Script to run simple linear regression on manual control between 2013-09-20 and 2014-01-01
modelID = 7;
blackrock = './blackrock/';
labview = './labview/';
threshold = 5;

%Fetch paramcode to load
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
paramcode = exec(conn, ['SELECT `description` FROM models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};

%Fetch each pair of nev files to run
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
toprocess = exec(conn, ['SELECT `manualrecording`, rec.`labview file`, `1DBCrecording`, `manualrecordingafter`, `dualrecording` FROM experiment_tuning et INNER JOIN `recordings` rec ON '...
	'et.`manualrecording` = rec.`nev file`']);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);
	
%Set to 1 with caution...
rerun = 1;

%MC1
for idx = 1:nR
	nevfile = toprocess{idx, 1};
	matfile = toprocess{idx, 2};
	BCnevfile = toprocess{idx, 3};
	bciunits = exec(conn, ['SELECT `unit` FROM bci_units WHERE `ID` = "' BCnevfile '"']);
	bciunits = fetch(bciunits);
	bciunits = bciunits.Data;
	otherunits = exec(conn, ['SELECT `unit` FROM units WHERE `nev file` = "' nevfile '" AND `firingrate` > '...
	 num2str(threshold)]);
	otherunits = fetch(otherunits); 
	otherunits = otherunits.Data;
	allunits = unique([otherunits; bciunits]);
	display(['Processing ' nevfile])
	if exist([blackrock nevfile], 'file')
		processLinearCursor(conn, modelID, blackrock, labview, nevfile, matfile, paramcode, threshold, allunits, rerun);
	else
		display('Cannot find file, continuing')
	end
end

%BC1
for idx = 1:nR
	nevfile = toprocess{idx, 3};
	matfile = toprocess{idx, 2};
	bciunits = exec(conn, ['SELECT `unit` FROM bci_units WHERE `ID` = "' nevfile '"']);
	bciunits = fetch(bciunits);
	bciunits = bciunits.Data;
	otherunits = exec(conn, ['SELECT `unit` FROM units WHERE `nev file` = "' nevfile '" AND `firingrate` > '...
	 num2str(threshold)]);
	otherunits = fetch(otherunits); 
	otherunits = otherunits.Data;
	allunits = unique([otherunits; bciunits]);
	display(['Processing ' nevfile])
	if exist([blackrock nevfile], 'file')
		processLinearCursor(conn, modelID, blackrock, labview, nevfile, matfile, paramcode, threshold, allunits, rerun);
	else
		display('Cannot find file, continuing')
	end
end

%MC2
for idx = 1:nR
	nevfile = toprocess{idx, 4};
	matfile = toprocess{idx, 2};
	BCnevfile = toprocess{idx, 3};
	bciunits = exec(conn, ['SELECT `unit` FROM bci_units WHERE `ID` = "' BCnevfile '"']);
	bciunits = fetch(bciunits);
	bciunits = bciunits.Data;
	otherunits = exec(conn, ['SELECT `unit` FROM units WHERE `nev file` = "' nevfile '" AND `firingrate` > '...
	 num2str(threshold)]);
	otherunits = fetch(otherunits); 
	otherunits = otherunits.Data;
	allunits = unique([otherunits; bciunits]);
	display(['Processing ' nevfile])
	if exist([blackrock nevfile], 'file')
		processLinearCursor(conn, modelID, blackrock, labview, nevfile, matfile, paramcode, threshold, allunits, rerun);
	else
		display('Cannot find file, continuing')
	end
end

%DC
for idx = 1:nR
	nevfile = toprocess{idx, 5};
	matfile = toprocess{idx, 2};
	BCnevfile = toprocess{idx, 3};
	bciunits = exec(conn, ['SELECT `unit` FROM bci_units WHERE `ID` = "' BCnevfile '"']);
	bciunits = fetch(bciunits);
	bciunits = bciunits.Data;
	otherunits = exec(conn, ['SELECT `unit` FROM units WHERE `nev file` = "' nevfile '" AND `firingrate` > '...
	 num2str(threshold)]);
	otherunits = fetch(otherunits); 
	otherunits = otherunits.Data;
	allunits = unique([otherunits; bciunits]);
	display(['Processing ' nevfile])
	if exist([blackrock nevfile], 'file')
		processLinearCursor(conn, modelID, blackrock, labview, nevfile, matfile, paramcode, threshold, allunits, rerun);
	else
		display('Cannot find file, continuing')
	end
end