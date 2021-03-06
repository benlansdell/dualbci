%Script to run Transfer Entropy model, tracking units between recordings. 
%Model 35 runs the covariance
modelID = 35;
blackrock = './blackrock/';
labviewpath = './labview/';

%Fetch paramcode to load
conn = db_conn();
paramcode = exec(conn, ['SELECT `description` FROM models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};

%Fetch each pair of nev files to run
conn = database('',databaseuser,databasepwd,'com.mysql.jdbc.Driver', ...
	databaseurl);
toprocess = exec(conn, ['SELECT `1DBCrecording`, `manualrecording`, `manualrecordingafter`, `dualrecording`, `experiment_id` FROM experiment_tuning']);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);

%rng('shuffle')
% 1-18 (americano) 19-45 (espresso) 46-67 (mocha) 68-89 (galao) 90-108 (mocha -- since latte sshfs not working...)
for idx = 105:108
	MCnevfile1 = toprocess{idx,2};
	MCnevfile2 = toprocess{idx,3};
	BCnevfile1 = toprocess{idx,1};
	DCnevfile = toprocess{idx,4};
	expt_id = toprocess{idx,5};
	if ~strcmp(MCnevfile1, 'null') & ~strcmp(MCnevfile2, 'null') & ~strcmp(BCnevfile1, 'null') & ~strcmp(DCnevfile, 'null')
		display(['Processing ' MCnevfile1])
		processCovariance(conn, modelID, blackrock, labviewpath, MCnevfile1, BCnevfile1, MCnevfile2, DCnevfile, expt_id, paramcode);
	else
		display('Can''t find all files, continuing')
		continue
	end
end