function weeklynevs = genweeklydataset(csv, condition, duration)
	%Produces a list of nev files of specified minimum duration and condition, one file for each week
	%of the year 2013.
	%
	%Usage:
	%	weeklynevs = genweeklydatasets(csv, condition, duration)
	%
	%Input:
	%	csv = filename to csv file output from nevmappings containing all nev and mapping info
	%	year = condition
	%	duration = number of seconds recording must be for. If nev file is below this amount,
	%		will move onto the next day until a file is found at least as long as this time.
	%
	%Output:
	%	weeklynevs = structure array containing a list of file names under requested conditions and duration
	%		.startdate
	%		.dur
	%		.date
	%		.nevfile
	%		.matfile
	%
	%Test code:
	%	csv = './nevmappings2.csv';
	%	condition = '2D Manual Position';
	%	duration = 360; %at least six minutes
	%	weeklynevs = genweeklydataset(csv, condition, duration);

	year = 2013;
	months = 1:12;
	%Every month
	%weekstarts = [1];
	%weekends = [31];
	%Every two weeks
	%weekstarts = [1 16];
	%weekends = [15 31];
	%Every week
	weekstarts = [1 8 16 23];
	weekends = [7 15 22 31];
	cnt = 1;
	for month = months
		for idx = 1:length(weekstarts)
			found = 0;
			days = weekstarts(idx):weekends(idx);
			date = [num2str(year) '/' num2str(month) '/' num2str(weekstarts(idx))];
			%Get all the nevs matching the condition for that week
			%This is terribly inefficient to sort through the csv file once for each week, rather
			%than once in total...
			nevs = retrievenev(csv, year, month, days, condition);
			%Check there is an nev file at the the duration needed
			for nev = nevs
				if nev.dur > duration
					found = 1;
					display([date ': ' nev.nevfile ' dur: ' num2str(nev.dur) ' is above min. dur'])
					break
				end
				display([date ': ' nev.nevfile ' dur: ' num2str(nev.dur)])
			end
			weeklynevs(cnt).startdate = date;
			if (found == 1)
				weeklynevs(cnt).date = nev.date;
				weeklynevs(cnt).nevfile = nev.nevfile;
				weeklynevs(cnt).matfile = nev.matfile;
				weeklynevs(cnt).dur = nev.dur;
			else
				display(['Cannot find nev file of min duration ' num2str(duration) '(s) for period starting ' date '.'])
			end
			cnt = cnt + 1;
		end
	end
end
