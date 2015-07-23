%Script to run simple linear regression on manual control between 2013-09-20 and 2014-01-01
modelID = 2;
blackrock = './blackrock/';
threshold = 5;
%after = '2013-12-05';
%before = '2014-01-01';
duration = 180;

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
toprocess = exec(conn, ['SELECT `File name MCP 1`, `File name MCP 2`,`File name MCV`,`File name BC` FROM ' tablename]);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
toprocess = reshape(toprocess,1,[]);
nR = size(toprocess,2);

%rng('shuffle')
for idx = 1:nR
	nevfile = toprocess{idx};
	display(['Processing ' nevfile])
	if ~strcmp(nevfile, 'null')
		processGrangerGLM(conn, modelID, blackrock, nevfile, paramcode, threshold);
	end
end

