%Script to run simple linear regression on manual control between 2013-09-20 and 2014-01-01
modelID = 1;
blackrock = './blackrock/';
labview = './labview/';
threshold = 5;
after = '2013-09-01';
before = '2014-09-20';
tasktype = 'dual';
duration = 240;

%Fetch paramcode to load
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
paramcode = exec(conn, ['SELECT `description` FROM models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};


%Fetch files to analyze
toprocess = exec(conn, ['SELECT `nev file`, `labview file` FROM recordings WHERE `nev date` BETWEEN "'...
 after '" AND "' before '" AND `tasktype` = "' tasktype '" AND `duration` > ' num2str(duration)]);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);
for idx = 604:nR
	nevfile = toprocess{idx, 1};
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
		processLinear(conn, modelID, blackrock, nevfile, paramcode, threshold, allunits);
	else
		display('Cannot find file, continuing')
	end
end