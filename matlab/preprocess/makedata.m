function dat = makedata(trials, spikes, cursor, target)

	nU = size(spikes, 2);
	dtrials = diff(trials);
	nztrials = find(trials==0);
	%dat = {};
	ntrials = 0;
	
	while length(trials > 0)
		startidx = find(dtrials==1,1);
		endidx = find(dtrials==-1,1);
		%start within trial
		if startidx > endidx
			startidx = 1;
		end
		%start within trial, no later trial start
		if isempty(startidx) & ~isempty(endidx)
			startidx = 1;
		end
		%end in trial
		if ~isempty(startidx) & isempty(endidx)
			endidx = length(trials);
		end
		if isempty(startidx) & isempty(endidx)
			%No changes, all within trial
			if trials(1) == 1
				startidx = 1;
				endidx = 1;
			%No changes, all outside trial
			else
				break
			end
		end
		ntrials = ntrials + 1;
		dat(ntrials).trialId = ntrials;
		indices = startidx+1:endidx;
		if nU > 1
			dat(ntrials).spikes = spikes(indices,:)';
			spikes = spikes(endidx+1:end,:);	
		else
			dat(ntrials).spikes = spikes(indices)';
			spikes = spikes(endidx+1:end);	
		end
		startPos = cursor(startidx+1,:);
		targetPos = target(startidx+1,:);
		dat(ntrials).quadrant = Targ2Quadrants(startPos, targetPos);
		cursor = cursor(endidx+1:end,:);
		target = target(endidx+1:end,:);
		trials = trials(endidx+1:end);
		dtrials = diff(trials);
	end
end