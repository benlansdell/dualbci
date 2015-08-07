
%%%%%%%%%%%%%
%Add Toffset%
%%%%%%%%%%%%%

%Add missing columns to recordings table
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', 'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db')
tablename = 'recordings';

%Toffset
colname = {'`Toffset`'};
toprocess = exec(conn,'select `nev file`, `labview file` from `recordings` WHERE `Toffset` IS NULL');
toprocess = fetch(toprocess);
nFiles = size(toprocess.Data);
if nFiles == 0
	display('No Toffset fields to update')
end

for idx = 1:nFiles
	%Find missing data in matfiles
	nevfile = toprocess.Data{idx,1};
	matfile = toprocess.Data{idx,2};
	if ~strcmp(matfile, 'null')
		display(['Processing ' nevfile])
		data = load(['./labview/' matfile]);
		for j = 1:length(data.data.nev)
			if strcmp(data.data.nev(j).nevfile, nevfile)
				Toffset = {data.data.nev(j).Toffset(1)};
				%Update
				whereclause = {['where `nev file` = ''' nevfile '''']};
				update(conn,tablename,colname,Toffset,whereclause);
				break
			end
		end
	else
		continue
	end
end

%Do the same with:

%%%%%%%%%%%%
%Add trials%
%%%%%%%%%%%%

%Performance. Trial info.
%Add missing columns to recordings table
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', 'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db')
tablename = 'trials';
colnames = {'`nev file`', '`labview file`', '`start`', '`end`', '`duration`', '`success`', '`valid`', '`startPosX`', '`startPosY`', '`targetPosX`', '`targetPosY`'};
nC = length(colnames);
matfiles = dir('./labview/*.mat');

%Iterate over all mat files
for i = 1:length(matfiles)
	%Add each trial info
	matfile = matfiles(i).name;

	%Determine if should add this mat file
	if findstr(matfile, 'Log')
		continue
	end

	%Determine if have already added this mat file
	toprocess = exec(conn,['SELECT * FROM `trials` WHERE `labview file` = ''' matfile '''']);
	toprocess = fetch(toprocess);
	nTrials = size(toprocess.Data,1);
	if nTrials == 0 | strcmp(toprocess.Data{1}, 'No Data')
		display([matfile ' trials have not been added, searching for trials to add.'])
	else
		display([matfile ' trials have already been added, skipping.'])
		continue
	end

	%Load nev file info
	nevinfo = exec(conn,['SELECT `nev file`, `starttime`, `endtime` FROM `recordings` WHERE `labview file` = ''' matfile '''']);
	nevinfo = fetch(nevinfo);

	nNev = size(nevinfo.Data,1);
	if nNev > 0 & ~strcmp(nevinfo.Data, 'No Data')
		nevrectimes = nevinfo.Data;
	else
		%If none then add no trials 
		display([matfile ' contains no recorded sessions'])
		continue
	end
	%nevrectimes

	load(matfile);
	nT = length(data.trials.time);
	newrows = {};
	%BCI units
	for j = 1:nT
		trialstart = data.trials.time(j);
		trialend = trialstart+data.trials.duration(j);
		%Figure out corresponding nev file, if any. 
		trialnev = findNev(nevrectimes, trialstart, trialend);
		%If none then skip
		if length(trialnev) > 0
			newrows = [newrows; {trialnev, matfile, trialstart, trialend, data.trials.duration(j), data.trials.success(j), data.trials.valid(j), data.trials.startPos(j,1), data.trials.startPos(j,2), data.trials.targetPos(j,1), data.trials.targetPos(j,2)}];
		end
	end

	datainsert(conn,tablename,colnames,newrows)
	display(['Added ' num2str(size(newrows, 1)) ' trials.'])
	clear data
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add BCI units to recordings%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', 'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db')

tablename = 'recordings';

%Toffset
colnames = {'`ch1`', '`ch2`', '`ch3`', '`ch4`'};
toprocess = exec(conn,'select `nev file`, `labview file` from `recordings` WHERE `ch1` IS NULL');
toprocess = fetch(toprocess);
nFiles = size(toprocess.Data,1);
if nFiles == 0
	display('No BCI channel fields to update')
end

for idx = 1:nFiles
	%Find missing data in matfiles
	nevfile = toprocess.Data{idx,1};
	matfile = toprocess.Data{idx,2};
	if ~strcmp(matfile, 'null')
		display(['Processing ' nevfile])
		data = load(['./labview/' matfile]);
		for j = 1:length(data.data.nev)
			if strcmp(data.data.nev(j).nevfile, nevfile)
				%Get channels
				chans = num2cell(data.data.nev(j).chans);
				%Update
				whereclause = {['where `nev file` = ''' nevfile '''']};
				update(conn,tablename,colnames,chans,whereclause);
				break
			end
		end
	else
		continue
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Add the number of units above 5Hz%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

threshold = 5;
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', 'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db')
tablename = 'recordings';

%Toffset
colname = {'`abovefive`'};
toprocess = exec(conn,'select `nev file` from `recordings` WHERE `abovefive` IS NULL');
toprocess = fetch(toprocess);
nFiles = size(toprocess.Data,1);
if nFiles == 0
	display('No abovefive fields to update')
end

for idx = 1:nFiles
	%Find missing data in matfiles
	nevfile = toprocess.Data{idx,1};
	nevpath = ['./blackrock/' nevfile];
	display(['Processing ' nevfile])
	NEV = openNEV(nevpath, 'nosave');
	nevsamplerate = NEV.MetaTags.TimeRes;
	display('Making structures')
	%Find the duration and sample rate of the nev file recording
	dur = NEV.MetaTags.DataDuration/nevsamplerate;
	%Convert spike times into array of binned spikes, one for each spike sorted channel
	nE = 128;
	nunits = 5; 
	nU = nE*nunits;
	spiketimes = double(NEV.Data.Spikes.TimeStamp)/nevsamplerate;
	elecs = cell(1,nU);
	spikemuas = struct('times', elecs);
	unitnames = cell(1,nU);
	averate = zeros(1,nU);
	isvalid = zeros(1,nU);
	for idx=1:nU
		spikemuas(idx).times = [0];    
	end
	display('Sorting into units')
	for i=1:length(spiketimes)
		E = NEV.Data.Spikes.Electrode(i);
		unit = NEV.Data.Spikes.Unit(i);
		U = single((E-1)*nunits)+single(unit)+1;
		spikemuas(U).times = [spikemuas(U).times; spiketimes(i)];
		unitnames{U} = [num2str(E) '.' num2str(unit)];
	end
	%Check which channels are doing stuff
	display('Checking units are not inactive for large parts of recording')
	for idx=1:nU
		averate(idx) = (length(spikemuas(idx).times)-1)/dur;
		if length(spikemuas(idx).times)>1
			if (spikemuas(idx).times(2)<20) & (spikemuas(idx).times(end)>(dur-20))
				isvalid(idx)=1;
			end
		end
	end
	%Set a threshold firing rate, below which we ignore that unit
	display('Counting spikes')
	abovethresh = (averate > threshold) & isvalid;
	nU = {sum(abovethresh)};
	display(['Found ' num2str(nU{1}) ' above ' num2str(threshold) 'Hz.'])
	%Write to database...
	whereclause = {['where `nev file` = ''' nevfile '''']};
	update(conn,tablename,colname,nU,whereclause);
end

%%%%%%%%%%%%%%%%%%%%%%%%
%Add intertrial periods%
%%%%%%%%%%%%%%%%%%%%%%%%

testnev = '20121205SpankyUtah001.nev';

conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', 'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db')
tablename = 'intertrials';
colnames = {'`nev file`', '`start`', '`end`', '`duration`'};

%Find all nev files without intertrials
nevfiles = exec(conn, ['SELECT `nev file`, `starttime`, `endtime` FROM `recordings` WHERE' ...
' not exists (select * from intertrials where `intertrials`.`nev file`=`recordings`.`nev file`)'...
' and duration is not null']);
nevfiles = fetch(nevfiles);
nN = size(nevfiles.Data, 1);

%Iterate over all nev files to process
for i = 1:nN
	%Add each trial info
	%nevfile = testnev;
	nevfile = nevfiles.Data{i};
	nevstart = nevfiles.Data{i, 2};
	nevend = nevfiles.Data{i, 3};
	display(['Processing ' nevfile]);
	%Fetch all trials within nev file
	trials = exec(conn, ['SELECT `start`, `end` FROM `trials` WHERE' ...
	' `nev file`="' nevfile '"']);
	trials = fetch(trials);
	trials = trials.Data;
	nT = size(trials, 1);
	if strcmp('No Data', trials)
		nT = 0;
	end
	newrows = {};
	%If no trials then add whole recording
	if nT == 0
		newrows = {nevfile, nevstart, nevend, nevend - nevstart};
	%Otherwise
	else
		%Add before first trial and after recording start (if trial starts after recording starts)
		trialstart = trials{1,1};
		if trialstart > nevstart
			newrows = vertcat(newrows, {nevfile, nevstart, trialstart, trialstart-nevstart});
		end
		for j = 2:nT
			%Iterate over remaining trials, adding period before current trial start and after previous trial end
			prevtrialend = trials{j-1,2};
			trialstart = trials{j,1};
			newrows = vertcat(newrows, {nevfile, prevtrialend, trialstart, trialstart-prevtrialend});
		end
		%If last trial ends before recording end, add period between last trial end and recording end
		trialend = trials{end,2};
		if trialend < nevend
			newrows = vertcat(newrows, {nevfile, trialend, nevend, nevend-trialend});
		end
	end
	datainsert(conn,tablename,colnames,newrows)
	display(['Added ' num2str(size(newrows, 1)) ' inter-trial periods.'])
end