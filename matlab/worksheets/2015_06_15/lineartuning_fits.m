%Script to run simple linear regression on manual control between 2013-09-20 and 2014-01-01
modelID = 1;
blackrock = './blackrock/';
threshold = 5;
after = '2014-01-01';
before = '2014-09-20';
tasktype = 'brain';
duration = 180;

%Fetch paramcode to load
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/Spanky');
paramcode = exec(conn, ['SELECT `description` FROM Models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};

%Fetch files to analyze
%Fetch all files 
toprocess = exec(conn, ['SELECT `nev file` FROM Recordings WHERE `nev date` BETWEEN "'...
 after '" AND "' before '" AND `tasktype` = "' tasktype '" AND `duration` > ' num2str(duration)]);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);
for idx = 1:nR
	nevfile = toprocess{idx, 1};
	display(['Processing ' nevfile])
	if exist([blackrock nevfile], 'file')
		processLinear(conn, modelID, blackrock, nevfile, paramcode, threshold);
	else
		display('Cannot find file, continuing')
	end
end