%State variables
RESTART_TRIAL = 0;
RUNNING = 1;
IN_TARGET = 2;
SUCCESS = 3;
FAIL = 4;
INTERTRIAL = 5;
SAMPLERATE = 60;

%Find DC recordings in recordings.csv and in experiments_tuning_dual.csv
recs_fn = './experiment_tuning_dual.csv';
recs_info_fn = './recordings.csv';

intertrial_time = 0.5;
success_time = 1.0;
timeout_time = 40.0;
new_target_radius = 0.12;
success_radius = 0.12;

intertrial_bin = round(intertrial_time*SAMPLERATE);
success_bin = round(success_time*SAMPLERATE);
timeout_bin = round(timeout_time*SAMPLERATE);

recs = readtable(recs_fn);
recs_info = readtable(recs_info_fn);
N = 100;

M = size(recs_info, 1);

null_success_rates = [];
actual_success_rates = [];

for idx = 1:N
	nev_file = recs{idx,1};
	mat_file = recs{idx,2};
	display(['Running new rec'])
	%Soooooooooo inefficient
	for rec_idx = 1:M
		nf = recs_info{rec_idx,4};
		mf = recs_info{rec_idx,1};
		if strcmp(nf, nev_file) & strcmp(mf, mat_file)
			break
		end
	end
	%Didn't find... continued to next nev file
	if rec_idx == M
		continue
	end
	%Find recording info for that nevfile
	start_time = str2num(cell2mat(recs_info{rec_idx, 6}));
	end_time = str2num(cell2mat(recs_info{rec_idx, 7}));
	n_act_success = str2num(cell2mat(recs_info{rec_idx, 15}));
	rec_time = str2num(cell2mat(recs_info{rec_idx, 8}));
	%Load the matfile
	data = load(cell2mat(mat_file));
	%Extract cursor information from the right time
	start_bin = round(start_time*SAMPLERATE);
	end_bin = round(end_time*SAMPLERATE);
	cursor = data.data.stateHist.cursor(start_bin:end_bin,:);

	%Run the trial logic for that whole recording
	state = RESTART_TRIAL;
	newstate = RUNNING;
	current_running_bins = 0;
	in_target_bins = 0;
	n_success = 0;
	n_trials = 0;
	success_times = [];
	target = [0,0];
	K = size(cursor, 1);
	for j = 1:K
		%Decide on the current state
		if (state == RESTART_TRIAL)
			%Restart counters
			current_running_bins = 0;
			in_target_bins = 0;
			%Pick a new target randomly outside radius of current position
			found_target = false;
			while (found_target == false)
				target = rand(1,2)-0.5;
				if (sqrt((cursor(j,1)-target(1))^2 + (cursor(j,2)-target(2))^2) > new_target_radius)
					found_target = true;
				end
			end
			newstate = RUNNING;
%		end
		elseif (state == RUNNING)
			in_target_bins = 0;
			current_running_bins = current_running_bins + 1;
			%If cursor is in target radius change state to IN_TARGET
			if (sqrt((cursor(j,1)-target(1))^2 + (cursor(j,2)-target(2))^2) < success_radius)
				newstate = IN_TARGET;
			end
			%If trial time is too long then fail
			if (current_running_bins > timeout_bin)
				newstate = FAIL;
			end
		elseif (state == IN_TARGET)
			in_target_bins = in_target_bins + 1;
			current_running_bins = current_running_bins + 1;
			%If in_target_bins is long enough then success
			if (in_target_bins > success_bin)
				%display('Success!')
				newstate = SUCCESS;
			end
			%If trial time is too long then fail
			if (current_running_bins > timeout_bin)
				%display('Timeout. Fail!')
				newstate = FAIL;
			end
			%If cursor is out of target radius change state to RUNNING
			if (sqrt((cursor(j,1)-target(1))^2 + (cursor(j,2)-target(2))^2) > success_radius)
				newstate = RUNNING;
			end
		elseif (state == SUCCESS)
			n_success = n_success + 1;
			n_trials = n_trials + 1;
			success_times(end+1) = current_running_bins/SAMPLERATE;
			current_running_bins = 0;
			newstate = INTERTRIAL;
		elseif (state == FAIL)
			n_trials = n_trials + 1;
			current_running_bins = 0;
			newstate = INTERTRIAL;
		elseif (state == INTERTRIAL)
			current_running_bins = current_running_bins + 1;			
			if (current_running_bins > intertrial_bin)
				newstate = RESTART_TRIAL;
			end
		end
		state = newstate;
	end
	%Finished... display the results
	null_success_rate = n_success/rec_time;
	actual_success_rate = n_act_success/rec_time;

	display(['Actual success rate ' num2str(actual_success_rate)])
	display(['Null success rate ' num2str(null_success_rate)])

	null_success_rates(end+1) = n_success/rec_time;
	actual_success_rates(end+1) = n_act_success/rec_time;
end

%null_success_rates = sort(null_success_rates);
alpha_null = max(null_success_rates);

%Number of actual success rates above alpha
sum(actual_success_rates>alpha_null)
sum(actual_success_rates>null_success_rates)