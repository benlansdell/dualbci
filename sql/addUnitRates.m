conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', 'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db')
threshold = 0;
toprocess = exec(conn,'select `nev file` from `recordings` WHERE `nev date` > "2013-12-06"');
toprocess = fetch(toprocess);
nFiles = size(toprocess.Data,1);
tablename = 'units';
colnames = {'`nev file`', 'unit', 'firingrate'};

cannotfind = {};

for idx = 1:nFiles
	%Find missing data in matfiles
	nevfile = toprocess.Data{idx,1};
	nevpath = ['./blackrock/' nevfile];
	display(['Processing ' nevfile])
	if ~exist(nevpath, 'file')
		display(['Cannot find ' nevfile])
		cannotfind = {cannotfind{:}, nevfile};
		continue
	end
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
	abovethresh = (averate > threshold) & isvalid;
	nU = sum(abovethresh);
	%Add also the number of units above 5Hz
	abovefive = (averate > 5) & isvalid;
	nU5 = sum(abovefive);
	display([num2str(nU) ' units above ' num2str(threshold) 'Hz in ' nevfile])
	unitnames = unitnames(abovethresh);
	averate = averate(abovethresh);
	%Write to database...
	for idx = 1:length(unitnames)
		unit = unitnames{idx};
		firingrate = averate(idx);
		data = {nevfile, unit, firingrate};
		datainsert(conn,tablename,colnames,data);
	end
	whereclause = ['WHERE `nev file` = "' nevfile '"'];
	update(conn, 'recordings', {'abovefive'}, {nU5}, whereclause)	
	display([num2str(nU5) ' units above 5Hz in ' nevfile])
end