%Script to run simple linear and GLM decoding algorithms for comparison
modelID = 12;
verbose = 1;
blackrock = './blackrock/';
labview = './labview/';
scriptname = './2015_07_28-LinearVsGLMDecoding/glm_linear_decoding.m';
scriptdesc = ['Comparison between linear and GLM decoding performance on MC dataset.'];
%Fetch paramcode to load
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
paramcode = exec(conn, ['SELECT `description` FROM models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};
%Add job to things running
id = logJob(conn, scriptname, scriptdesc);
%Fetch files to analyze
tablename = '`experiment_tuning`';
toprocess = exec(conn, ['SELECT `manualrecording` FROM ' tablename]);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);
for idx = 2:nR
	nevfile = toprocess{idx};
	display(['Processing ' nevfile])
	if exist([blackrock nevfile], 'file')
		processGLMLinearDecoding(conn, modelID, blackrock, nevfile, paramcode, verbose);
	else
		display('Cannot find file, continuing')
	end
end
%Finished successfully (?), close job
closeJob(conn, id);