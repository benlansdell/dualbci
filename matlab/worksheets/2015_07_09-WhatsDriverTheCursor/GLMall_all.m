%Script to run simple linear regression on manual control between 2013-09-20 and 2014-01-01
modelID = 14;
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
toprocess = exec(conn, ['SELECT `1DBCrecording`, `manualrecording`, `manualrecordingafter` FROM ' tablename]);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
toprocess = reshape(toprocess, [], 1);
nR = size(toprocess,1);
eval(paramcode);

for idx = 252:nR
	nevfile = toprocess{idx};
	%Figure out which units to process, that haven't already been processed by this model
	units = exec(conn, ['SELECT u.`unit` FROM Units u'...
		' WHERE u.`nev file` = "' nevfile '" AND u.`firingrate` > ' num2str(threshold)...
		' AND NOT EXISTS (SELECT * FROM Fits f WHERE f.`unit` = u.`unit` AND'...
		' f.modelID = ' num2str(modelID) ' AND f.`nev file` = u.`nev file`)'	]);
	units = fetch(units);
	units = units.Data;
	display(['Processing ' nevfile])
	if exist([blackrock nevfile], 'file')
		if ~strcmp(units, 'No Data')
			processGLM(conn, modelID, blackrock, labviewpath, nevfile, paramcode, units);
		else
			display('Already processed')
		end
	else
		display('Cannot find file, continuing')
	end
end

