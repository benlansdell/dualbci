%Add missing columns to bci_units table
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
tablenameUnits = 'bci_units';
colnames = {'ID', 'unit', 'direction', 'angle'};

%Get a unique list of mat files
%toprocess = exec(conn,['select DISTINCT `labview file` from recordings where not exists'...
% ' (select * from bci_units where ID = `nev file`) and `labview file`="Spanky_2013-03-19-1408.mat"']);
toprocess = exec(conn,['select DISTINCT `labview file` from recordings where not exists'...
 ' (select * from bci_units where ID = `nev file`) AND `nev date` > "2013-12-05"']);
toprocess = fetch(toprocess);
nFiles = size(toprocess.Data,1);
if nFiles == 0
	display('No labview files to process')
end

for idx = 1:nFiles
	matfile = toprocess.Data{idx};
	if strcmp(matfile, 'null')
		continue
	end
	display(['Processing ' matfile])
	%Get the corresponding nev files, etc, associated with it
	nevfiles = exec(conn,['select `nev file`, `labview task type`, `starttime`, `endtime`, `ch1`,' ...
	' `ch2`, `ch3`, `ch4`, `duration` from recordings'...
	 ' where `labview file`="' matfile '"']);
	nevfiles = fetch(nevfiles);
	%Open the mat file
	load(['./labview/' matfile]);
	pvmap = data.pvmap;
	%For each nev, find all settings that overlap with it (including the one before it)
	for j = 1:size(nevfiles.Data,1)
		settings = {};
		nS = size(pvmap.time,1);
		nevfile = nevfiles.Data{j,1};

		%display(['Processing ' nevfile]);
		starttime = nevfiles.Data{j,3};
		endtime = nevfiles.Data{j,4};
		dur = nevfiles.Data{j,9};
		%Iterate over settings
		for k = 1:nS
			%If setting was made before end of nev recording, it was possibly 
			%in effect, otherwise it definitely wasn't
			if pvmap.time(k) < endtime;
				%If it was the only setting (nS == 1), or this is the last setting, then it was in effect
				if k == nS
					%Proportion of total time this setting was in effect
					proptime = (endtime - max(starttime, pvmap.time(k)))/dur;
					%Add settings
					settings = vertcat(settings, {nevfiles.Data{j,1}, pvmap.time(k), pvmap.angles(k,:), pvmap.enabled(k,:), proptime});
				%Otherwise we have to do more checking
				else
					%If the setting after this one was after the start of the recording, then
					%this setting was in place when the recording started...
					if pvmap.time(k+1) > starttime
						proptime = (min(endtime, pvmap.time(k+1)) - max(pvmap.time(k), starttime))/dur;
						settings = vertcat(settings, {nevfiles.Data{j,1}, pvmap.time(k), pvmap.angles(k,:), pvmap.enabled(k,:), proptime});
					end
				end
			end
		end
		%Go through settings and remove changes that don't change the mapping
		changes = {};
		for k = 1:size(settings,1)
			if isempty(changes)
				changes = settings(k, :);
			else
				%If the current settings are different in any way to the previous ones...
				currenabled = settings{k,4};
				prevenabled = changes{end,4};
				currdirs = settings{k,3}(currenabled == 1);
				prevdirs = changes{end,3}(prevenabled == 1);
				if any(currenabled ~= prevenabled) | any(currdirs ~= prevdirs)
					changes = vertcat(changes, settings(k,:));
				%otherwise, if they're the same, add the proptime to the last actual changes...
				else
					changes{end,5} = changes{end,5} + settings{k,5};
				end
			end
		end
		%If the settings haven't change in important ways within the recording, add the details to the 
		%database
		%More specifically, if there exists a mapping which took up more than 95% of the recording then use
		%that...
		unambigmapping = {};
		for k = 1:size(changes, 1)
			if changes{k,5}>0.95
				unambigmapping = changes(k,:);
				break
			end
		end
		%No unambiguous settings to use :(
		if isempty(unambigmapping)
			if dur > 20
				display(['No unambiguous settings found for ' nevfile '(dur=' num2str(dur) ...
					'). Charlie''s code labels file as ' nevfiles.Data{j,2}])
				changes	
			end
			%Add to conflict table, if not already in there
			isconflict = fetch(exec(conn, ['SELECT ID FROM `conflicted_recordings` WHERE ID = "' nevfile '"']));
			isconflict = isconflict.Data;
			if strcmp(isconflict, 'No Data')
				datainsert(conn, 'conflicted_recordings', {'ID'}, {nevfile});
			end
		else
			%Determine direction for BCI units, and add these...
			for k = 1:4
				chan = nevfiles.Data{j, 4+k};
				%Determine its direction
				theta = unambigmapping{3}(k);
				switch theta
				case 180
					direction = 'west';
					thta = pi;
				case 0
					direction = 'east';
					thta = 0;
				case 90
					direction = 'north';
					thta = pi/2;
				case 270
					direction = 'south';
					thta = 3*pi/2;
				end
				enabled = unambigmapping{4}(k);
				%If the channel is enabled add it
				if enabled & (chan > 0)
					mappedunit = {nevfile, chan, direction, thta};
					datainsert(conn, tablenameUnits, colnames, mappedunit);
				end
			end
		end
	end
end

%Were these added twice? Yes :(
%Processing Spanky_2013-12-04-1239.mat
%Processing Spanky_2013-12-05-1400.mat
