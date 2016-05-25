modelID = 34;
blackrock = './blackrock/';
labviewpath = './labview/';

%Fetch paramcode to load
conn = db_conn();
paramcode = exec(conn, ['SELECT `description` FROM models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};

%Fetch each pair of nev files to run
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
toprocess = exec(conn, ['SELECT `1DBCrecording`, `manualrecording`, `dualrecording`, `manualrecordingafter`, `experiment_id` FROM experiment_tuning']);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);

basename = './worksheets/2016_02_07-GPFApaired/run';

%rng('shuffle')
%for idx = 1:nR
for idx = 1:nR

	MCnevfile1 = toprocess{idx,2};
	BCnevfile1 = toprocess{idx,1};
	DCnevfile = toprocess{idx,3};
	MCnevfile2 = toprocess{idx,4};
	expt_id = toprocess{idx,5};

	%Get labview file from DB
	lbfile = exec(conn, ['SELECT `labview file` FROM recordings WHERE `nev file` = "' MCnevfile1 '"']);
	lbfile = fetch(lbfile);
	lbfile = lbfile.Data{1};

	if ~strcmp(MCnevfile1, 'null') & ~strcmp(MCnevfile2, 'null') & ~strcmp(BCnevfile1, 'null') & ~strcmp(DCnevfile, 'null')
		display(['Processing ' MCnevfile1])
		processGPFAPaired(conn, modelID, blackrock, labviewpath, lbfile, MCnevfile1, BCnevfile1, DCnevfile, MCnevfile2, expt_id, paramcode, basename);
	else
		display('Can''t find all files, continuing')
		continue
	end
end