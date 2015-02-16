function trialinfo(matpath, fn_out_event, fn_out_all, matfile_out)
    %Create a csv file which lists, for each date, the no. of hours of recording
    %that have taken place under each condition (2D manual, brain control, etc) upto that date.
    %And also extract, for each electrode, how many hours it has been used to drive a BCI
    %
    %Usage:
    %   trialinfo(matpath, fn_out_event, fn_out_all)
    %
    %Input:
    %   matpath = filename to csv file output from nevmappings containing all nev and mapping info
    %   fn_out_event = output csv file for only each date with recordings
    %   fn_out_all = output csv file for recording info for all days of the year 
    %   matfile_out = output mat file containing events structure
    %
    %Output:
    %   none
    %   Format of file written to, each column is:
    %       date number (using MATLAB's datenum function),
    %       1D Horiz Brain  Position rec. seconds,
    %       1D Horiz Brain  Velocity rec. seconds,
    %       1D Horiz Manual Position rec. seconds,
    %       1D Horiz Manual Velocity rec. seconds,
    %       1D Vert  Brain  Position rec. seconds,
    %       1D Vert  Brain  Velocity rec. seconds,
    %       1D Vert  Manual Position rec. seconds,
    %       1D Vert  Manual Velocity rec. seconds,
    %       2D Brain  Velocity rec. seconds,
    %       2D Manual Position rec. seconds,
    %       2D Manual Velocity rec. seconds,
    %       Dual Control rec. seconds,
    %   Following that are sets of 1x128 vectors containing how many seconds each electrode 
    %   has been used for different BCI conditions. These vectors are:
    %       BCIHpos, BCIVpos, BCIHvel, BCIVvel, BCI2vel, BCIDual  
    %
    %Test code:
    %   matpath = './labview/';
    %   fn_out_event = './trialinfo_events.csv';
    %   fn_out_all = './trialinfo_all.csv';
    %   trialinfo(matpath, fn_out_event, fn_out_all);

    if (nargin < 1) matpath = './labview/'; end
    if (nargin < 2) fn_out_event = './trialinfo_event.csv'; end
    if (nargin < 3) fn_out_all = './trialinfo_all.csv'; end

    %Find all .mat files in matpath directory
    matfiles = dir([matpath '/*.mat']);
    nE = 128;
    k = 1;
    keys = {'1D Horiz Brain  Position',    '1D Horiz Brain  Velocity',    '1D Horiz Manual Position',    '1D Horiz Manual Velocity',    '1D Vert  Brain  Position',    '1D Vert  Brain  Velocity',    '1D Vert  Manual Position',    '1D Vert  Manual Velocity',    '2D Brain  Velocity',    '2D Manual Position',    '2D Manual Velocity', 'Dual Control'};
    vals = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    mappings = containers.Map(keys, vals);
    for i = 1:length(matfiles)
        %For each mat file, load it, see if it has an nev structure, 
        %extract the recording info from it, write file, clear structure
        fn = matfiles(i).name
        clear data;
        load([matpath '/' fn]);
        if exist('data', 'var')
            if isfield(data, 'nev')
            for j = 1:length(data.nev)
                nevfile = data.nev(j).nevfile;
                matfile = fn;
                nevdur = data.nev(j).DurationSec;
                nevmap = data.nev(j).map;
                nevoffset = single(data.nev(j).Toffset(1))/60;
                activechannels = data.nev(j).chans;
                if strcmp(nevmap, 'Error') | strcmp(nevmap, '') | strcmp(nevmap, 'Multiple') | strcmp(nevmap, '(waitfile)')
                    continue
                end
                %Found an nev file whose info we want to record
                nevdate = strrep(data.fname, 'Spanky_', '');
                nevtime = nevdate(end-3:end);
                %Save info in an events structure
                events(k).date = nevdate;
                events(k).daten = datenum(nevdate(1:end-5), 'yyyy-mm-dd')
                events(k).mapping = nevmap;
                events(k).totaldur = nevdur;
                [bciunits, durs] = findBCIUnits(data, nevmap, activechannels, nevoffset, nevdur);
                events(k).bciunits = bciunits;
                events(k).dur = durs;
                events(k).matfile = matfile;
                events(k).nevfile = nevfile;
                events(k).nevtime = nevtime;
                k = k + 1;
            end
            end
        end
    end

    %Save events structure
    save(matfile_out, 'events');

    %%%%%%%%%%%%
    %Write file%
    %%%%%%%%%%%%
    fh_event = fopen(fn_out_event, 'w');
    fh_all = fopen(fn_out_all, 'w');
    %Create list of all dates to be written to file
    %date format: 2013-01-02-1318
    datestart = datenum('2013-01-01', 'yyyy-mm-dd');
    dateend = datenum('2013-12-31', 'yyyy-mm-dd');
    %Initialize counts of things...
    %       1D Horiz Brain  Position rec. seconds,
    %       1D Horiz Brain  Velocity rec. seconds,
    %       1D Horiz Manual Position rec. seconds,
    %       1D Horiz Manual Velocity rec. seconds,
    %       1D Vert  Brain  Position rec. seconds,
    %       1D Vert  Brain  Velocity rec. seconds,
    %       1D Vert  Manual Position rec. seconds,
    %       1D Vert  Manual Velocity rec. seconds,
    %       2D Brain  Velocity rec. seconds,
    %       2D Manual Position rec. seconds,
    %       2D Manual Velocity rec. seconds,
    %       Dual Control rec. seconds,
    cond_secs        = zeros(1,12);
    electrodeBCIHpos = zeros(1,nE);
    electrodeBCIVpos = zeros(1,nE);
    electrodeBCIHvel = zeros(1,nE);
    electrodeBCIVvel = zeros(1,nE);
    electrodeBCI2vel = zeros(1,nE);
    electrodeBCIDual = zeros(1,nE);    
    curr_date = datestart;
    %Iterate over events structure
    for idx = 1:length(events)
        %Write number of lines for intervening days between event day and 'curr_date' into all file
        writeLines(fh_all, curr_date, events(idx).daten, cond_secs, electrodeBCIHpos, electrodeBCIVpos, electrodeBCIHvel, electrodeBCIVvel, electrodeBCI2vel, electrodeBCIDual);
        %Update all info
        curr_date = events(idx).daten;
        map = events(idx).mapping;
        totaldur = events(idx).totaldur;
        for j = 1:length(events(idx).dur)
            units = events(idx).bciunits{j};
            units = units(units>0);
            dur = events(idx).dur(j);
            if strcmp(map, '1D Horiz Brain  Position')
                electrodeBCIHpos(units) = electrodeBCIHpos(units) + dur;
            elseif strcmp(map, '1D Horiz Brain  Velocity')
                electrodeBCIHvel(units) = electrodeBCIHvel(units) + dur;
            elseif strcmp(map, '1D Vert  Brain  Position')
                electrodeBCIVpos(units) = electrodeBCIVpos(units) + dur;
            elseif strcmp(map, '1D Vert  Brain  Velocity')
                electrodeBCIVvel(units) = electrodeBCIVvel(units) + dur;
            elseif strcmp(map, '2D Brain  Velocity')
                electrodeBCI2vel(units) = electrodeBCI2vel(units) + dur;
            elseif strcmp(map, 'Dual Control')
                electrodeBCIDual(units) = electrodeBCIDual(units) + dur;
            end
        end
        cond_secs(mappings(map)) = cond_secs(mappings(map)) + totaldur;
        %Write single line with updated info into event file
        writeLine(fh_event, curr_date, cond_secs, electrodeBCIHpos, electrodeBCIVpos, electrodeBCIHvel, electrodeBCIVvel, electrodeBCI2vel, electrodeBCIDual);
    end
    %Finish off all dates file
    writeLines(fh_all, curr_date, dateend, cond_secs, electrodeBCIHpos, electrodeBCIVpos, electrodeBCIHvel, electrodeBCIVvel, electrodeBCI2vel, electrodeBCIDual);
    %Close
    fclose(fh_event);
    fclose(fh_all);
end 

function writeLine(fh, curr_date, cond_secs, electrodeBCIHpos, electrodeBCIVpos, electrodeBCIHvel, electrodeBCIVvel, electrodeBCI2vel, electrodeBCIDual, nevfile, matfile)
    %Includes only current date
    %curr_date
    fprintf(fh, '%d,', curr_date);
    %cond_secs
    for idx = 1:length(cond_secs)
        fprintf(fh, '%f,', cond_secs(idx));
    end
    %electrodeBCIHpos, 
    for idx = 1:length(electrodeBCIHpos)
        fprintf(fh, '%f,', electrodeBCIHpos(idx));
    end
    %electrodeBCIVpos, 
    for idx = 1:length(electrodeBCIVpos)
        fprintf(fh, '%f,', electrodeBCIVpos(idx));
    end
    %electrodeBCIHvel, 
    for idx = 1:length(electrodeBCIHvel)
        fprintf(fh, '%f,', electrodeBCIHvel(idx));
    end
    %electrodeBCIVvel, 
    for idx = 1:length(electrodeBCIVvel)
        fprintf(fh, '%f,', electrodeBCIVvel(idx));
    end
    %electrodeBCI2vel, 
    for idx = 1:length(electrodeBCI2vel)
        fprintf(fh, '%f,', electrodeBCI2vel(idx));
    end
    %electrodeBCIDual
    for idx = 1:length(electrodeBCIDual)
        fprintf(fh, '%f,', electrodeBCIDual(idx));
    end
    fprintf(fh, '\n');
end

function writeLines(fh, startdate, enddate, cond_secs, electrodeBCIHpos, electrodeBCIVpos, electrodeBCIHvel, electrodeBCIVvel, electrodeBCI2vel, electrodeBCIDual)
    %Includes startdate but not enddate
    for curr_date = startdate:(enddate-1)
        writeLine(fh, curr_date, cond_secs, electrodeBCIHpos, electrodeBCIVpos, electrodeBCIHvel, electrodeBCIVvel, electrodeBCI2vel, electrodeBCIDual);
    end
end

function [units, durs] = findBCIUnits(data, mapping, activechannels, offset, dur)
    units = {};
    durs = [];
    %If manual just return an empty list, otherwise, find which were used based on data.pvmap
    if isempty(strfind(mapping, 'Manual'))
        %Find the options that were used during this nev file
        optidx1st = 0;
        optidx2nd = 0;
        for idx = 1:length(data.pvmap.time)
            if data.pvmap.time(idx)<offset
                optidx1st = idx;
            end
            if data.pvmap.time(idx)<offset+dur
                optidx2nd = idx;
            end
        end
        optidx1st = max(optidx1st,1);
        optidx2nd = max(optidx2nd,1);
        %If the options changed within the recording then we need to take more care...
        if optidx1st ~= optidx2nd
            for idx = optidx1st:optidx2nd
                enabledchans = data.pvmap.enabled(idx,1:4);
                units{idx-optidx1st+1} = floor(activechannels(enabledchans==1));
                %1st time
                if idx == optidx1st
                    segdur = data.pvmap.time(idx+1)-offset;
                %Last time
                elseif idx == optidx2nd
                    segdur = offset+dur-data.pvmap.time(idx);
                %Other times
                else
                    segdur = data.pvmap.time(idx+1)-data.pvmap.time(idx);
                end
                durs(idx-optidx1st+1) = segdur;
            end
        %Otherwise we don't
        else
            enabledchans = data.pvmap.enabled(optidx1st,1:4);
            units = {floor(activechannels(enabledchans==1))};
            durs = [dur];
        end
    end
end