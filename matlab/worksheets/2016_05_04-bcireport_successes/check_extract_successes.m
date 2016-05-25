modelID = 40;
conn = database('','root','Fairbanks1!','com.mysql.jdbc.Driver', ...
	'jdbc:mysql://fairbanks.amath.washington.edu:3306/spanky_db');
paramcode = exec(conn, ['SELECT `description` FROM models WHERE modelID = ' num2str(modelID)]);
paramcode = fetch(paramcode);
paramcode = paramcode.Data{1};
eval(paramcode);

nevfile = '20130920SpankyUtah001.nev';
nevpath = ['./blackrock/' nevfile];
matfile = './labview/Spanky_2013-09-20-1356.mat';

fn_out = './worksheets/diagnostics/plots/test_extract_successes.eps';
verbose = 1;

processed = preprocess_spline_lv(nevpath, matfile, binsize, threshold, offset, fn_out, verbose);
processed_suc = extractSuccesses(processed, conn, nevfile, 0);

load(matfile);
labviewsamplerate = data.sampleRate;
foundnev = 0;
[pathstr,nevfilename,ext] = fileparts(nevfile);

for idx = 1:length(data.nev)
	%Find nev structure
	if length(data.nev(idx).nevfile)==0
		continue
	end
	[pathstr,curr_nevfilename,ext] = fileparts(data.nev(idx).nevfile);
	if strcmp(curr_nevfilename, nevfilename)
		%Find offset and duration in seconds
		nevoffset = single(data.nev(idx).Toffset(1))/60;
		nevdur = data.nev(idx).DurationSec;
		foundnev = 1;
		break
	end 
end
assert(foundnev == 1, 'Cannot find nev file in labview mat file, check your input files');

trials = import_trials(matfile);
%Find trials in nev file 
nevtrials = [];
for idx = 1:size(trials,2)
	if trials(idx).starttime > nevoffset & trials(idx).endtime < (nevoffset + nevdur)
		nevtrials = [nevtrials trials(idx)];
	end
end

%Plot cursor data 
nB = 1000;
dur = nB*binsize;
clf
subplot(2,1,1)
plot(processed.cursor(1:nB,1));
hold on 

%Label trial times 
for idx = 1:size(nevtrials,2)
	st = nevtrials(idx).starttime
	et = nevtrials(idx).endtime
	suc = nevtrials(idx).success;
	if (st > nevoffset) & (et < (nevoffset + dur))
		if suc == 1 
			c = 'r';
		else
			c = 'b';
		end
		targx = nevtrials(idx).target(1);
		plot([(st - nevoffset)/binsize, (st - nevoffset)/binsize], [-.5 .5], c, 'linewidth', 2)
		plot([(et - nevoffset)/binsize, (et - nevoffset)/binsize], [-.5 .5], c, 'linewidth', 2)
		plot([(st - nevoffset)/binsize, (et - nevoffset)/binsize], [targx targx], '--k')
	end
end
title('All cursor data. Red marks successful trial start/end points, dotted black marks target location')
ylabel('Cursor x')

%Plot processes_suc 
subplot(2,1,2)
plot(processed_suc.cursor(1:nB));
title('Successful trials extracted')
xlabel('time')
ylabel('Cursor x')
saveplot(gcf, './worksheets/2016_05_04-bcireport_successes/test_extract_success.eps')