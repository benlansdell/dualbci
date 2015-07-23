%Script to run simple linear regression on manual control between 2013-09-20 and 2014-01-01
modelID = 9;
blackrock = './blackrock/';
labviewpath = './labview/';

%Fetch paramcode to load
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/Spanky');
paramcode = exec(conn, ['SELECT `description` FROM Models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};

%Fetch each pair of nev files to run
tablename = '`Analysis Feb-Aug 2014 Tuning Change`';
colnames = {'`File name MCP 1`', '`File name MCP 2`','`File name MCV`','`File name BC`'};
toprocess = exec(conn, ['SELECT `File name BC`,`File name MCP 1`, `File name MCP 2` FROM ' tablename]);
%toprocess = exec(conn, ['SELECT `1DBCrecording` FROM ' tablename]);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);

for idx = 2:nR
	BCnevfile = toprocess{idx,1};
	nevfile1 = toprocess{idx,2};
	nevfile2 = toprocess{idx,3};
	display(['Processing ' nevfile1])
	if exist([blackrock nevfile1], 'file')
		processMVGCmanual(conn, modelID, blackrock, labviewpath, nevfile1, BCnevfile, paramcode);
	else
		display('Cannot find file, continuing')
	end
	display(['Processing ' nevfile2])
	if exist([blackrock nevfile2], 'file')
		processMVGCmanual(conn, modelID, blackrock, labviewpath, nevfile2, BCnevfile, paramcode);
	else
		display('Cannot find file, continuing')
	end
end

