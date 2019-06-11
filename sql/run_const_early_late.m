%Script to run Transfer Entropy model, tracking units between recordings. 
%Model 45 runs the const function
modelID = 45;
blackrock = './blackrock/';
labviewpath = './labview/';

%Fetch paramcode to load
conn = database('',databaseuser,databasepwd,'com.mysql.jdbc.Driver', ...
	databaseurl);
paramcode = exec(conn, ['SELECT `description` FROM models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};

%Fetch each pair of nev files to run
toprocess = exec(conn, ['SELECT `1DBCrecording`, `manualrecording`, `manualrecordingafter`, `dualrecording`, `experiment_id` FROM experiment_tuning']);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);

for idx = 1:108
	BCnevfile1 = toprocess{idx,1};
	MCnevfile1 = toprocess{idx,2};
	MCnevfile2 = toprocess{idx,3};
	DCnevfile = toprocess{idx,4};
	expt_id = toprocess{idx,5};
	if ~strcmp(MCnevfile1, 'null') & ~strcmp(MCnevfile2, 'null') & ~strcmp(BCnevfile1, 'null') & ~strcmp(DCnevfile, 'null')
		display(['Processing ' BCnevfile1])
		processLinearConstEarlyLate(conn, modelID, blackrock, labviewpath, MCnevfile1, paramcode);
		processLinearConstEarlyLate(conn, modelID, blackrock, labviewpath, BCnevfile1, paramcode);
		processLinearConstEarlyLate(conn, modelID, blackrock, labviewpath, DCnevfile, paramcode);
		processLinearConstEarlyLate(conn, modelID, blackrock, labviewpath, MCnevfile2, paramcode);
	else
		display('Can''t find all files, continuing')
		continue
	end
end