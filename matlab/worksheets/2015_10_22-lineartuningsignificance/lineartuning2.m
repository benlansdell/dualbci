%Script to run simple linear regression on manual control between 2013-09-20 and 2014-01-01
modelID = 1;
blackrock = './blackrock/';
labviewpath = './labview/';

%Fetch paramcode to load
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
paramcode = exec(conn, ['SELECT `description` FROM models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};

%Fetch each pair of nev files to run
tablename = '`experiment_tuning`';
colnames = '`1DBCrecording`, `manualrecording`, `experiment_id`';
toprocess = exec(conn, ['SELECT ' colnames ' FROM ' tablename]);
%toprocess = exec(conn, ['SELECT `1DBCrecording` FROM ' tablename]);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);

for idx = 1:nR
	BCnevfile = toprocess{idx,1};
	nevfile1 = toprocess{idx,2};
	expt_id = toprocess{idx,3};
	display(['Processing ' nevfile1])
	if exist([blackrock nevfile1], 'file') & exist([blackrock BCnevfile], 'file'))
		allunits = fetch(exec(conn, ['SELECT u1.unit FROM units u1 '...
		'INNER JOIN units u2 '...
		'ON u1.unit = u2.unit '...
		'WHERE u2.`nev file` = "' nevfile1 '" '...
		'AND u1.`nev file` = "' BCnevfile '" '...
		'AND u1.`firingrate` > 5 AND u2.`firingrate` > 5']));	
		allunits = allunits.Data;	
		processLinear(conn, modelID, blackrock, BCnevfile, paramcode, threshold, allunits);
		processLinear(conn, modelID, blackrock, nevfile1, paramcode, threshold, allunits);
	else
		display('Cannot find all files, continuing')
	end
end

closeJob(conn, id);
