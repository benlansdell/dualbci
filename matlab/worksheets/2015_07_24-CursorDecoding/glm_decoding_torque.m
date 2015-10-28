%Script to run simple linear regression on manual control between 2013-09-20 and 2014-01-01
modelID = 11;
verbose = 1;
blackrock = './blackrock/';
labview = './labview/';
scriptname = './2015_07_24-CursorDecoding/glm_decoding_torque.m';
scriptdesc = ['GLM decoding of torque values based on GLM w spike history filter '...
		' and cursor position fitlers. Applied to 1D brain control dataset.'];
%Fetch paramcode to load
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
paramcode = exec(conn, ['SELECT `description` FROM models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};
%Add job to things running
%id = logJob(conn, scriptname, scriptdesc);
%Fetch files to analyze
tablename = '`experiment_tuning`';
toprocess = exec(conn, ['SELECT `1DBCrecording`, `manualrecording`, `manualrecordingafter`, `experiment_id` FROM ' tablename]);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);
for idx = 1:nR
	nevfile2 = toprocess{idx, 3};
	nevfile = toprocess{idx, 2};
	BCnevfile = toprocess{idx, 1};
	expt_id = toprocess{idx,4};
	display(['Processing ' nevfile])
	if exist([blackrock nevfile], 'file')
		processGLMDecodingTorque(conn, modelID, blackrock, labview, nevfile, BCnevfile, nevfile2, expt_id, paramcode, verbose);
	else
		display('Cannot find file, continuing')
	end
end
%Finished successfully (?), close job
%closeJob(conn, id);