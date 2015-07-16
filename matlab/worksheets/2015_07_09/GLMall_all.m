%Script to run simple linear regression on manual control between 2013-09-20 and 2014-01-01
modelID = 5;
blackrock = './blackrock/';
labviewpath = './labview/';

%Fetch paramcode to load
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/Spanky');
paramcode = exec(conn, ['SELECT `description` FROM Models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};

%Fetch each pair of nev files to run
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/Spanky');
tablename = '`AnalysisLinear`';
toprocess = exec(conn, ['SELECT `1DBCrecording`, `manualrecording` FROM ' tablename]);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
toprocess = reshape(toprocess, [], 1);
nR = size(toprocess,1);

for idx = 1:nR
	nevfile = toprocess{idx};
	display(['Processing ' nevfile])
	if exist([blackrock nevfile], 'file')
		processGLM(conn, modelID, blackrock, labviewpath, nevfile, paramcode);
	else
		display('Cannot find file, continuing')
	end
end

