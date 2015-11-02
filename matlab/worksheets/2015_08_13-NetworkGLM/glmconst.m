%Script to run Granger causality model with 'causal' and 'acausal' stim filters. 
modelID = 31;
blackrock = './blackrock/';

%Fetch paramcode to load
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
paramcode = exec(conn, ['SELECT `description` FROM models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};

%Fetch each pair of nev files to run
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
toprocess = exec(conn, ['SELECT `manualrecording`, `1DBCrecording`, `manualrecordingafter`, `dualrecording`, `1DBCrecordingafter` FROM experiment_tuning']);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
toprocess = reshape(toprocess,1,[]);
toprocess = unique(toprocess);
nR = size(toprocess,2);

%rng('shuffle')
for idx = 1:nR
	nevfile = toprocess{idx};
	display(['Processing ' nevfile])
	if ~strcmp(nevfile, 'null')
		processGLM(conn, modelID, blackrock, labviewpath, nevfile, paramcode);
	end
end

