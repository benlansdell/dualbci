%Script to run Granger causality model, tracking units between recordings. 
%Model 29 fits GC over a longer set of training data
modelID = 33;
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
toprocess = exec(conn, ['SELECT `1DBCrecording`, `manualrecording`, `dualrecording`, `experiment_id` FROM experiment_tuning']);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);

%rng('shuffle')
%for idx = 1:nR
a = {};
for idx = 1:5

	MCnevfile1 = toprocess{idx,2};
	BCnevfile1 = toprocess{idx,1};
	DCnevfile = toprocess{idx,3};
	expt_id = toprocess{idx,4};

	%Get labview file from DB
	lbfile = exec(conn, ['SELECT `labview file` FROM recordings WHERE `nev file` = "' MCnevfile1 '"']);
	lbfile = fetch(lbfile);
	lbfile = lbfile.Data{1};

	if ~strcmp(MCnevfile1, 'null') & ~strcmp(BCnevfile1, 'null') & ~strcmp(DCnevfile, 'null')
		display(['Processing ' MCnevfile1])
		[id, ll, D] = processGPFAIntrinsicDimDim(conn, modelID, blackrock, labviewpath, lbfile, MCnevfile1, BCnevfile1, DCnevfile, expt_id, paramcode);
		fn_out = ['./worksheets/2016_02_07-GPFApaired/IntDim_id_' num2str(id) '.mat']
		a{idx} = cumsum(D.^2)/sum(D.^2)
	else
		display('Can''t find all files, continuing')
		continue
	end
end

a_ave = [];
for idx = 1:5
	a_ave(idx) = a{idx}(5);
end

mean(a_ave)