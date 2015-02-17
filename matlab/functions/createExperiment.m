function expts = createExperiment(matfile_in, ndays, condition)
	%Find all consecutive sets of recordings for a specified duration and condition
	% and extract the trial information for what happened in between the two recordings.
	%This can be used to fit a model for the start and the end of the expt and try to understand
	%the differences between the two as a function of the trials that took place during that time
    %
    %Usage:
    %   expts = createExperiment(matfile_in, ndays, condition)
    %
    %Input:
    %   matpath = filename to csv file output from nevmappings containing all nev and mapping info
    %   ndays = number of days between which to compare recordings
    %   condition = string of recording condition to find recordings for
    %
    %Output:
    %	expts = structure array listing experiment details:
    %		.nevfile_start = nevfile before
    %		.nevfile_end = nevfile after
    %		.matfile_start = matfile before
    %		.matfile_end = matfile after
    %		.cond_secs = [1xnC] array listing number of seconds of trials under each condition 
    %		.events = details about all the events that come in between start and end, including BCI units used
    %			between start and end recording -- what we'll be fitting a model to.
    %
    %Test code:
    %   matpath = './labview/';
    %   fn_out_event = './trialinfo_events.csv';
    %   fn_out_all = './trialinfo_all.csv';
    %   matfile_out = './trialinfo.mat';
    %   %Only run this if needed
    %   trialinfo(matpath, fn_out_event, fn_out_all, matfile_out);
    %   condition = '2D Manual Position';
    %   ndays = 1;
    %   expts = createExperiment(matfile_out, ndays, condition);
    %   save('./expts/2dmanualpos_1day.mat', 'expts');

    keys = {'1D Horiz Brain  Position',    '1D Horiz Brain  Velocity',    '1D Horiz Manual Position',    '1D Horiz Manual Velocity',    '1D Vert  Brain  Position',    '1D Vert  Brain  Velocity',    '1D Vert  Manual Position',    '1D Vert  Manual Velocity',    '2D Brain  Velocity',    '2D Manual Position',    '2D Manual Velocity', 'Dual Control'};
    vals = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    mappings = containers.Map(keys, vals);
    nC = 12;

	load(matfile_in);
	nexpts = 1;
	inexpt = false;
	%Iterate over events structure
	for idx = 1:length(allevents)
		event = allevents(idx);
        currdate = event.daten;
		%If in expt
		if inexpt == true
			%If within end of timeframe (<ndays after first day)
            if currdate < exptenddate
                %--Append info of event to expt struct
                expts(nexpts).cond_secs(mappings(event.mapping)) = expts(nexpts).cond_secs(mappings(event.mapping)) + event.totaldur;
                expts(nexpts).events = {expts(nexpts).events{:}, event};                    
            %If at end of timeframe (==ndays after first day)
            elseif currdate == exptenddate
                %--If come across condition of interest, use that as the terminal condition and finish expt struct,
                %  set inexpt = false
                if strcmp(event.mapping, condition)
                    expts(nexpts).nevfile_end = event.nevfile;
                    expts(nexpts).matfile_end = event.matfile;
                    nexpts = nexpts + 1;
                    inexpt = false;
                %--else append info of event to expt struct
                else
                    expts(nexpts).cond_secs(mappings(event.mapping)) = expts(nexpts).cond_secs(mappings(event.mapping)) + event.totaldur;
                    expts(nexpts).events = {expts(nexpts).events{:}, event};                    
                end
            %If outside of timeframe (>ndays after first day), 
            %cancel this one as timing conditions aren't met.
            else
				inexpt = false;
    	    end
        end
    	%If not in expt. Note that the same event can count as the end of an expt and the start of one...
		if inexpt == false
			%If find recording that matches condition, then we are in expt!
			if strcmp(event.mapping,condition)
				inexpt = true;
    			expts(nexpts).nevfile_start = event.nevfile;
    			expts(nexpts).matfile_start = event.matfile;
    			expts(nexpts).nevfile_end = '';
    			expts(nexpts).matfile_end = '';
    			expts(nexpts).cond_secs = zeros(1,nC);
    			expts(nexpts).events = {};
                expts(nexpts).startdate = event.daten;
                expts(nexpts).enddate = event.daten+ndays;
                exptenddate = event.daten+ndays;
    		end
		end
	end
end