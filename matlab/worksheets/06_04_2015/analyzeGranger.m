%Fetch each pair of nev files to run
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
tablename = 'experiment_tuning';
colnames = {'1DBCrecording', 'manualrecording'};
toprocess = exec(conn, ['SELECT `1DBCrecording`,`manualrecording` FROM ' tablename]);
toprocess = fetch(toprocess);
toprocess = toprocess.Data;
toprocess = reshape(toprocess,1,[]);
nR = size(toprocess,1);
files = toprocess;

performance = zeros(1,nR);
causaldensity = zeros(1,nR);

for idx = 1:nR
	conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
		'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
	nevfile = files{idx};
	perf = exec(conn, ['SELECT `successrate` FROM recordings WHERE `nev file` = "' nevfile '"']);
	perf = fetch(perf);
	performance(idx) = perf.Data{1};
	idx = 3; nevfile = files{idx};
	load(['./worksheets/06_04_2015/GC_' nevfile '.mat']);
	causaldensity(idx) = results.causaldensity;
end
