%
% USAGE: LogRead(baseName, logpath, matpath)
%
% Looks for baseName.csv and baseName.log and parses logged data. See
% README.TXT in the LabVIEW Docs directory for information on log formats.
% Parsed data is returned in the 'data' struct variable, and 
% also saved to disk as baseName.mat containing that variable.

function data = LogRead(file, logpath, matpath, dbFlag) 

%% Log code definitions - see README.txt in LabVIEW code

statehistcode = 0;
messagecode = 1;
expcode = 2;
paramcode = 3;
trialcode = 4;
targetcode = 7;
pvmapcode = 5;
rewardcode = 6;
filtercode = 9;


%% Check to see if .mat file already exists & up to date
fileName = strcat(matpath,filesep,file,'.mat');

fnameBase =  regexprep(file,'.[^.]*','',2);     % strip extension

codeInfo = dir(which('LogRead'));

if (exist(fileName) && ~dbFlag)
    matInfo = dir(fileName);
    fprintf('Found mat file, checking date... ')
    if (datenum(codeInfo.date) < datenum(matInfo.date))
        fprintf('.mat file already up to date\n')
        if nargout == 1
            data = load([logpath filesep file '.mat'])
        end
        return
    end
end

%% Still executing? Parse the logfiles    

fprintf('\nNo .mat found, or out of date. Need to parse original logfiles.\n')
fprintf('Please wait, reading .csv file...');

csvFile = strcat(logpath,'/',fnameBase,'.csv');
rawdata = dlmread(csvFile);
csvInfo = dir(csvFile);

%%
disp('Separating state history from events....')
stateHist = rawdata(rawdata(:,2) == statehistcode,:);
events = rawdata(rawdata(:,2) ~= statehistcode,:);
clear data

tstart = events(1,1);   % time of very first event, used as offset

%stateHist(:,1) = (stateHist(:,1)-tstart)/1000;   % convert to seconds into exp

disp('Parsing events...')

%neurons = 4;    % TODO: get this from input map once implemented
%data.neurons = neurons;


% FILTER KERNELS
data.filters = events(events(:,2) == filtercode,:);
filters.time = (data.filters(:,1)-tstart)/1000;   % convert to seconds into exp;
filters.index = data.filters(:,3);
filters.kernel = data.filters(:,4:end);
% crop filter kernel by assuming only nonzero kernel sample values
%for i = 1:length(filters.kernel(:,1))
%    filters.kernel(i) = filters.kernel(i,filters.kernel(i,:) ~= 0); 
%end
data.filters = filters;

% BMI / POPULATION VECTOR MAPPING 
data.pvmap = events(events(:,2) == pvmapcode,1:end);

% pvmap params include 4 per neuron (offset, gain, angle, enable)
% + 2 output offset, 2 int disable, timestamp, decay value, log code

% Ok, silly ad-hoc logic
pvmapSize = find(data.pvmap(1,:), 1, 'last' );
switch (pvmapSize)  % last nonzero value equiv to length
    case (21)               % 4x 4 ch + 2 offset + time, decay, code = 21
        neurons = 4;     % prior to integrator disable feature
    case {29, 31}           % += 4x 2 ch, += 2 int disable
        neurons = 6;     % after int. disable & inclusion of 2 torque channels
    case {35,36}               % += 4x 1 ch
        neurons = 7;     % after inclusion of left torque channel (7ch total)
    otherwise
        neurons = 4;
        disp('Warning, guessed at 4 neurons');
        keyboard
end
data.neurons = neurons;

pvmap.time = (data.pvmap(:,1)-tstart)/1000;   % convert to seconds into exp
% data.pvmap(:,2) = log code
pvmap.offsets = data.pvmap(:,2+(1:neurons));
pvmap.gains = data.pvmap(:,(2+neurons)+(1:neurons));
pvmap.angles = data.pvmap(:,(2+2*neurons)+(1:neurons));
pvmap.outputOffsets = data.pvmap(:,(2+3*neurons)+(1:2));
pvmap.enabled = data.pvmap(:,(4+3*neurons)+(1:neurons));
pvmap.decay = data.pvmap(:,(4+4*neurons)+1);
if (pvmapSize > 29)
    pvmap.intDisable = data.pvmap(:,(5+4*neurons)+(1:2));
else  % int disable feature was not avaliable, so they are all effectively 0
    pvmap.intDisable = 0*pvmap.outputOffsets;   
end

% one-sided control feature added 2013-04-12
if (pvmapSize == 36)
    pvmap.oneSided = data.pvmap(:,end);
else
    pvmap.oneSided = zeros(length(data.pvmap),1);
end
data.pvmap = pvmap;

data.sampleRate = round(1000/mean(diff(stateHist(:,1))));
%data.neurons = max(data.filters(:,3))-1;    % 1st index 0, total # = neurons +2


% TELEMETRY
stateHist = stateHist(:,1:(12 + 2*data.neurons) );
%data.stateHist.raw = stateHist(:,4:end);     %% BLOATS FILES - for troubleshooting only 
data.stateHist.time = (stateHist(:,1)-tstart)/1000; % convert to seconds into exp;
data.stateHist.lag = stateHist(:,3);
data.stateHist.spikes = stateHist(:,(0:neurons-1)+4);

% in newer version of LabVIEW code, spike count per bin is replaced by
% instantaneous spike rate sample, hence it's greater by a factor of SR

if (max(max(data.stateHist.spikes)) > 10)
    data.stateHist.spikes = data.stateHist.spikes/data.sampleRate;
end

if (neurons < 7)
% This part works for dates < 9/21/2012
data.stateHist.rates = stateHist(:,(0:neurons-1)+4+neurons);
data.stateHist.targetVisible = stateHist(:,end);
data.stateHist.error = stateHist(:,end-2:end-1);
data.stateHist.cursor = stateHist(:,end-4:end-3);
data.stateHist.velocity = stateHist(:,end-6:end-5);

% NOTE this doesn't handle manual channels well
else
% The following is valid after 9/21/2012
data.stateHist.mapped = stateHist(:,(0:neurons-1)+4+neurons);
data.stateHist.rates = stateHist(:,(0:neurons-1)+4+2*neurons);
data.stateHist.velocity = stateHist(:,end-3:end-2);
data.stateHist.error = stateHist(:,end-1:end);
end

% EXPERIMENT EVENTS
data.expts = events(events(:,2) == expcode,1:6);
expts.time = (data.expts(:,1)-tstart)/1000;   % convert to seconds into exp;
expts.code = data.expts(:,3);
data.expts = expts;

% (COMPLETE) TRIAL DATA AT END OF TRIAL
data.trials = events(events(:,2) == trialcode,1:10);
trials.time = (data.trials(:,1)-tstart)/1000;   % convert to seconds into exp;
trials.index = data.trials(:,3);
trials.startPos = data.trials(:,4:5);
trials.duration = data.trials(:,6);
trials.targetPos = data.trials(:,7:8);
trials.threshold = data.trials(:,9);
trials.success = data.trials(:,10);
data.trials = trials;

% TEXT LOG MESSAGES
data.messages = events(events(:,2) == messagecode,:);
messages.time = (data.messages(:,1)-tstart)/1000;   % convert to seconds into exp;
messages.type = data.messages(:,3);
messages.data = data.messages(:,4:end);
data.messages = messages;

% REWARD EVENTS
data.rewards = events(events(:,2) == rewardcode,:);
rewards.time = (data.rewards(:,1)-tstart)/1000;   % convert to seconds into exp
rewards.prob = data.rewards(:,3);
rewards.given = data.rewards(:,4);
data.rewards = rewards;

% EXPERIMENT PARAMETERS
data.params = events(events(:,2) == paramcode,1:12);
params.raw = data.params(:,3:end);
params.time = (data.params(:,1)-tstart)/1000;   % convert to seconds into exp;
params.xrange = data.params(:,3:4);
params.yrange = data.params(:,5:6);
params.exclusionRadius = data.params(:,7);
params.successRadius = data.params(:,8);
params.targetRadius = data.params(:,9);
params.trialDelay = data.params(:,10);
params.timeLimit = data.params(:,11);
params.dwellTime = data.params(:,12);
data.params = params;

% apply the following correction only to logfiles prior to 12-23-2012
if (datenum([2012 12 22 00 00 00]) > datenum(csvInfo.date))
    data = FixCursorHistory(data, dbFlag);      % exactly what it says.
else
    if (dbFlag)
    fprintf('Logfile recorded after 2012-12-23, cursor valid as recorded\n')
    end
    % these offsets don't make a lot of sense, but work
    data.stateHist.cursor = stateHist(:,end-4:end-3);
    data = FixCursorHistory2(data, dbFlag);
end

disp(sprintf('Pulled data from %2.2f min long data file',(data.stateHist.time(end)-data.stateHist.time(1))/60))
% Pull out messages from human-readable file

data.messages.text = {};

disp('Reading .log file');
fid = fopen(strcat(logpath,filesep,fnameBase,'.log'));
%disp('System event record including user logs:')

 while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if ~isempty(strfind(tline, [' ',int2str(messagecode),': ']))
        tline = strtrim(tline);
        data.messages.text = [data.messages.text; {tline}];
        %messages = {messages{:} tline};
        %disp(tline)
        
        % fix short-sighted coding of pause events
    i = length(data.messages.text);
    if (strfind(cell2mat(data.messages.text(i)),'Run state = 1'))
        data.messages.type(i) = 1;
        data.messages.data(i,1) = 1;
    else if (strfind(cell2mat(data.messages.text(i)), 'Run state = 0'))
        data.messages.type(i) = 1;
        data.messages.data(i,1) = 0;
        else
            data.messages.type(i) = -1;
        end
    end


    end
 end
 
 if size(data.messages.time,1) ~= size(data.messages.type,1)
     fprintf('Entries in .log file inconsistent with .csv, double-check files\n')
     keyboard
 end

 fclose(fid);
 
data.fname = file;

save(fileName,'data');
fprintf('Done parsing logfiles for %s\n\n',file);

end
