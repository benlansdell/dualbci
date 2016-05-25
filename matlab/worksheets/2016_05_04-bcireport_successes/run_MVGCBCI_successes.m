%Script to run simple linear regression on manual control between 2013-09-20 and 2014-01-01
modelID = 42;
blackrock = './blackrock/';
labviewpath = './labview/';

%Fetch paramcode to load
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
paramcode = exec(conn, ['SELECT `description` FROM models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};

%Fetch each pair of nev files to run
toprocess = exec(conn, ['SELECT `1DBCrecording`, `dualrecording` FROM `experiment_tuning`']);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);

for idx = 2:nR
	BCnevfile = toprocess{idx,1};
	DCnevfile = toprocess{idx,2};
	%Check not labeled as waitfile:
	mapping2 = fetch(exec(conn, ['SELECT `labview task type` FROM recordings WHERE `nev file` = "' BCnevfile '"']));
	mapping2 = mapping2.Data{1};
	if strcmp(mapping2, '(waitfile)')
		display('Labeled as (waitfile), continuing')
		continue
	end
	if exist([blackrock BCnevfile], 'file') & exist([blackrock DCnevfile], 'file')
		display(['Processing ' BCnevfile])
		processMVGCBCISuccess(conn, modelID, blackrock, labviewpath, BCnevfile, paramcode);
		display(['Processing ' DCnevfile])
		processMVGCBCISuccess(conn, modelID, blackrock, labviewpath, DCnevfile, paramcode);
	else
		display('Cannot find file, continuing')
	end
end