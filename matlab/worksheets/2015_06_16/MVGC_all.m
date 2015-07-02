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
tablename = 'AnalysisLinear';
tablename = 'AnalysisLinear2';
toprocess = exec(conn, ['SELECT `1DBCrecording` FROM ' tablename]);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
nR = size(toprocess,1);

rng('shuffle')
for idx = 74:nR
	nevfile = toprocess{idx};
	display(['Processing ' nevfile])
	processMVGCBCI(conn, modelID, blackrock, labviewpath, nevfile, paramcode, threshold);
end

