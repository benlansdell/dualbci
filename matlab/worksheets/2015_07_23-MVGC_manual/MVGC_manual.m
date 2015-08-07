%Script to run simple linear regression on manual control between 2013-09-20 and 2014-01-01
modelID = 9;
blackrock = './blackrock/';
labviewpath = './labview/';
scriptname = './2015_07_23-MVGC_manual/MVGC_manual.m';
scriptdesc = 'Performing MVGC on cursor pos from manual control dataset to compare to MVGC on BC. Finish the job on the other files.';

%Fetch paramcode to load
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
paramcode = exec(conn, ['SELECT `description` FROM models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};
id = logJob(conn, scriptname, scriptdesc);

%Fetch each pair of nev files to run
tablename = '`experiment_tuning`';
colnames = '`1DBCrecording`, `manualrecording`, `manualrecordingafter`, `experiment_id`';
toprocess = exec(conn, ['SELECT ' colnames ' FROM ' tablename]);
%toprocess = exec(conn, ['SELECT `1DBCrecording` FROM ' tablename]);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);

for idx = 1:nR
	BCnevfile = toprocess{idx,1};
	nevfile1 = toprocess{idx,2};
	nevfile2 = toprocess{idx,3};
	expt_id = toprocess{idx,4};
	display(['Processing ' nevfile1])
	if exist([blackrock nevfile1], 'file')
		processMVGCmanual(conn, modelID, blackrock, labviewpath, nevfile1, BCnevfile, nevfile2, expt_id, paramcode);
	else
		display('Cannot find file, continuing')
	end
end

closeJob(conn, id);
