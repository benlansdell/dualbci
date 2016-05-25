%Script to run simple linear regression on manual control between 2013-09-20 and 2014-01-01
modelID = 40;
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
toprocess = exec(conn, ['SELECT `manualrecording`, rec.`labview file`, `1DBCrecording`, `dualrecording`, `manualrecordingafter` FROM experiment_tuning et INNER JOIN `recordings` rec ON '...
	'et.`manualrecording` = rec.`nev file`']);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);

%Set to 1 with caution...
rerun = 0;

eval(paramcode);

for idx = 1:nR
	nevfile = toprocess{idx, 1};
	matfile = toprocess{idx, 2};
	BCnevfile = toprocess{idx, 3};
	DCnevfile = toprocess{idx, 4};
	MC2nevfile = toprocess{idx, 5};

	%Figure out which units to use:
		%Figure out units to use... above 5hz in all recordings
	units = fetch(exec(conn, ['SELECT u1.unit FROM '...
	'`experiment_tuning` et1 '...
	'INNER JOIN `units` u1 '...
	'ON u1.`nev file` = et1.`manualrecording` '...
	'INNER JOIN `units` u2 '...
	'ON u2.`nev file` = et1.`1DBCrecording` '...
	'INNER JOIN `units` u3 '...
	'ON u3.`nev file` = et1.`manualrecordingafter` '...
	'INNER JOIN `units` u4 '...
	'ON u4.`nev file` = et1.`dualrecording` '...
	'WHERE u1.unit = u2.unit AND u2.unit = u3.unit AND u3.unit = u4.unit AND '...
	'u1.firingrate > ' num2str(threshold) ' AND u2.firingrate > ' num2str(threshold) ' AND u3.firingrate > ' num2str(threshold) ' AND '...
	'u4.firingrate > ' num2str(threshold) ' AND et1.`manualrecording` = "' nevfile '"']));
	otherunits = units.Data;
 
	%Make sure BC units from both dual and brain control are in there
	bciunits = exec(conn, ['SELECT `unit` FROM bci_units WHERE `ID` = "' BCnevfile '"']);
	bciunits = fetch(bciunits);
	bciunits = bciunits.Data;
	dualunits = exec(conn, ['SELECT `unit` FROM bci_units WHERE `ID` = "' DCnevfile '"']);
	dualunits = fetch(dualunits);
	dualunits = dualunits.Data;
	allunits = unique([otherunits; bciunits; dualunits]);

	display(['Processing ' nevfile])
	if exist([blackrock nevfile], 'file') & exist([blackrock BCnevfile], 'file') & exist([blackrock DCnevfile], 'file') & exist([blackrock MC2nevfile], 'file')
		processLinearSuccessVelocity(conn, modelID, blackrock, nevfile, paramcode, threshold,  allunits, rerun);
		processLinearSuccessVelocity(conn, modelID, blackrock, BCnevfile, paramcode, threshold, allunits, rerun);
		processLinearSuccessVelocity(conn, modelID, blackrock, DCnevfile, paramcode, threshold, allunits, rerun);
		processLinearSuccessVelocity(conn, modelID, blackrock, MC2nevfile, paramcode, threshold, allunits, rerun);
	else
		display('Cannot find file, continuing')
	end
end
