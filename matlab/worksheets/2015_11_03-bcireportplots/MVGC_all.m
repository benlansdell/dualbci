%Script to run simple linear regression on manual control between 2013-09-20 and 2014-01-01
modelID = 3;
blackrock = './blackrock/';
labviewpath = './labview/';
threshold = 3;

%Fetch paramcode to load
conn = db_conn();
paramcode = exec(conn, ['SELECT `description` FROM models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};

%Fetch each pair of nev files to run
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
toprocess = exec(conn, ['SELECT `dualrecording` FROM experiment_tuning']);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);

for idx = 1:nR
	nevfile = toprocess{idx};
	display(['Processing ' nevfile])
	if exist([blackrock nevfile], 'file')
		processMVGCBCI(conn, modelID, blackrock, labviewpath, nevfile, paramcode, threshold);
	else
		display('Cannot find file, continuing')
	end
end

