conn = database('spanky_db',databaseuser,databasepwd,'com.mysql.jdbc.Driver', ...
	databaseurl);

MCDClist = ['MC'; 'DC'];
MCrec = '20130117SpankyUtah001.nev';
DCrec = '20130409SpankyUtah016.nev';
maxlines = 15;
figure

for idx1 = 1:2

	MCDC = MCDClist(idx1,:);
	if MCDC == 'DC'
		nevfile = DCrec;
	else
		nevfile = MCrec;
	end
	
	a = exec(conn, ['SELECT `start`, `end`, `success`, `valid`, `startPosX`, `startPosY`, `targetPosX`,' ...
									 '`targetPosY` FROM `trials` WHERE `nev file` = "' nevfile '"']);
	a = fetch(a);
	trials = a.Data;
	nT = size(trials, 1);
	
	startidx = round(cell2mat(trials(:,1))*60);
	endidx = round(cell2mat(trials(:,2))*60);
	success = cell2mat(trials(:,3));
	valid = cell2mat(trials(:,4));
	startposx = cell2mat(trials(:,5));
	startposy = cell2mat(trials(:,6));
	targetposx = cell2mat(trials(:,7));
	targetposy = cell2mat(trials(:,8));
	
	%Find tOffset for the file...
	a = exec(conn, ['SELECT `labview file` FROM `recordings` WHERE `nev file` = "' nevfile '"']);
	a = fetch(a);
	file = a.Data{1};
	
	labviewfile = load(file);
	successtrials = {};
	failtrials = {};
	trialtimes = [];
	
	maxplt = 20;
	plt = 0;
	
	successtrials_short = {};
	successtrials_med = {};
	successtrials_long = {};
	
	%Find the percentiles for times....
	for idx = 1:nT
		trialtime = (endidx(idx) - startidx(idx))/60;
		trialtimes = [trialtimes, trialtime];
	end
	sorted_times = sort(trialtimes);
	sorted_times = sorted_times(sorted_times~=0);
	n_times = length(sorted_times);
	shorttime = sorted_times(round(n_times*0.3));
	medtime = sorted_times(round(n_times*0.6));
	longtime = sorted_times(round(n_times*0.9));
	display(['Short time: ' num2str(shorttime) ' med time: ' num2str(medtime) ' long time: ' num2str(longtime)]);


	for idx = 1:nT
		trialtime = (endidx(idx) - startidx(idx))/60;
		if valid(idx)
			if success(idx)
				cursor = labviewfile.data.stateHist.cursor(startidx(idx):endidx(idx), :);
				%Scale and rotate by target position
				if length(cursor) > 0
					starttotarget = [startposx(idx) - targetposx(idx), startposy(idx) - targetposy(idx)];
					starttotarget_innerproduct = startposx(idx) - targetposx(idx);
					len = sqrt((startposx(idx) - targetposx(idx))^2 + (startposy(idx) - targetposy(idx))^2);
					cursor(:,1) = (cursor(:,1) - startposx(idx))/len;
					cursor(:,2) = (cursor(:,2) - startposy(idx))/len;
					%Rotate 
					theta = acos(starttotarget_innerproduct/len);
					if starttotarget(2) < 0
						theta = -theta;
					end
					R = [[cos(theta), -sin(theta)]; [sin(theta), cos(theta)]];
					rotated_cursor = cursor*R;
					if trialtime < shorttime
						successtrials_short{idx} = rotated_cursor;
					elseif trialtime < medtime
						successtrials_med{idx} = rotated_cursor;
					elseif trialtime < longtime
						successtrials_long{idx} = rotated_cursor;
					end
				end
			else
				cursor = labviewfile.data.stateHist.cursor(startidx(idx):endidx(idx), :);
				if length(cursor) > 0
					starttotarget = [startposx(idx) - targetposx(idx), startposy(idx) - targetposy(idx)];
					starttotarget_innerproduct = startposx(idx) - targetposx(idx);
					len = sqrt((startposx(idx) - targetposx(idx))^2 + (startposy(idx) - targetposy(idx))^2);
					cursor(:,1) = (cursor(:,1) - startposx(idx))/len;
					cursor(:,2) = (cursor(:,2) - startposy(idx))/len;
					%Rotate 
					theta = acos(starttotarget_innerproduct/len);
					if starttotarget(2) < 0
						theta = -theta;
					end
					R = [[cos(theta), -sin(theta)]; [sin(theta), cos(theta)]];
					rotated_cursor = cursor*R;
					failtrials{idx} = rotated_cursor;
				end
			end
		end
	end

	offset = idx1;
	subplot(3,2,offset)
	hold on
	l = 0;
	for idx = 1:length(successtrials_short)
		if length(successtrials_short{idx})
			l = l + 1;
			if l < maxlines
				plot(successtrials_short{idx}(:,1), successtrials_short{idx}(:,2), 'Color', [0.6 0.6 0.6])
				xlim([-2, 1])
				ylim([-2, 2])
				plot(successtrials_short{idx}(end,1), successtrials_short{idx}(end,2), 'r*')
				plot(successtrials_short{idx}(1,1), successtrials_short{idx}(1,2), 'k*')
			end
		end
	end
	if idx1 == 1
		ylabel('y')
	end
	
	subplot(3,2,offset+2)
	hold on
	l = 0;
	for idx = 1:length(successtrials_med)
		if length(successtrials_med{idx})
			l = l + 1;
			if l < maxlines
				plot(successtrials_med{idx}(:,1), successtrials_med{idx}(:,2), 'Color', [0.6 0.6 0.6])
				xlim([-2, 1])
				ylim([-2, 2])
				plot(successtrials_med{idx}(end,1), successtrials_med{idx}(end,2), 'r*')
				plot(successtrials_med{idx}(1,1), successtrials_med{idx}(1,2), 'k*')
			end
		end
	end
	if idx1 == 1
		ylabel('y')
	end
	
	subplot(3,2,offset+4)
	hold on
	l = 0;
	for idx = 1:length(successtrials_long)
		if length(successtrials_long{idx})
			l = l + 1;
			if l < maxlines
				plot(successtrials_long{idx}(:,1), successtrials_long{idx}(:,2), 'Color', [0.6 0.6 0.6])
				xlim([-2, 1])
				ylim([-2, 2])
				plot(successtrials_long{idx}(end,1), successtrials_long{idx}(end,2), 'r*')
				plot(successtrials_long{idx}(1,1), successtrials_long{idx}(1,2), 'k*')
			end
		end
	end
	xlabel('x')
	if idx1 == 1
		ylabel('y')
	end
end

saveplot(gcf, './figures/fig_1traj_examples.eps')