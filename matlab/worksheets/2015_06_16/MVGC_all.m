%Script to run simple linear regression on manual control between 2013-09-20 and 2014-01-01
modelID = 3;
blackrock = './blackrock/';
labviewpath = './labview/';
threshold = 3;
after = '2013-09-20';
before = '2014-01-01';

%Fetch paramcode to load
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/Spanky');
paramcode = exec(conn, ['SELECT `description` FROM Models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};

%Fetch each pair of nev files to run
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/Spanky');
tablename = '`Analysis Feb-Aug 2014 Tuning Change`';
colnames = {'`File name MCP 1`', '`File name MCP 2`','`File name MCV`','`File name BC`'};
toprocess = exec(conn, ['SELECT `File name BC` FROM ' tablename]);
%toprocess = exec(conn, ['SELECT `1DBCrecording` FROM ' tablename]);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);

rng('shuffle')
for idx = 34:nR
	nevfile = toprocess{idx};
	display(['Processing ' nevfile])
	if exist([blackrock nevfile], 'file')
		processMVGCBCI(conn, modelID, blackrock, labviewpath, nevfile, paramcode, threshold);
	else
		display('Cannot find file, continuing')
	end
end

