%Script to run simple linear regression on manual control between 2013-09-20 and 2014-01-01
modelID = 22;
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
toprocess = exec(conn, ['SELECT `manualrecording`, `1DBCrecording`, `manualrecordingafter` FROM experiment_tuning']);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
toprocess = reshape(toprocess,1,[]);
nR = size(toprocess,2);

for idx = 75:nR
	nevfile = toprocess{idx};
	display(['Processing ' nevfile])
	if ~strcmp(nevfile, 'null')
		processGrangerGLMShuffled(conn, modelID, blackrock, nevfile, paramcode);
	end
end

