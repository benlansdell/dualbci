csv = './nevmappings2.csv';
duration = 360; %at least six minutes
fn_out = './worksheets/01_13_2014/semiweeklynevs.mat';
conds = {'1D Horiz Brain  Position',
    '1D Horiz Brain  Velocity',
    '1D Horiz Manual Position',
    '1D Horiz Manual Velocity',
    '1D Vert  Brain  Position',
    '1D Vert  Brain  Velocity',
    '1D Vert  Manual Position',
    '1D Vert  Manual Velocity',
    '2D Brain  Velocity',
    '2D Manual Position',
    '2D Manual Velocity',
    'Dual Control'};
weeklynevs = containers.Map;
for idx = 1:length(conds)
	condition = conds{idx};
    weeklynevs(condition) = genweeklydataset(csv, condition, duration);
end
save(fn_out, 'weeklynevs');

recordingtimes = containers.Map;
for idx = 1:length(conds)
    condition = conds{idx};
    a.dates = [];
    a.durs = [];
    recordingtimes(condition) = a;
end
%Make plots of number of minutes of condition performed as a function of date
%Read in each line of nevmappings file and add recording duration to condition vector,
%along with time of recording
fh = fopen(csv, 'r');
tline = fgetl(fh);
while ischar(tline)
    %display(tline)
    words = strsplit(tline,',');
    if ~strcmp(words{1}, 'nev file')
        datetime = words{3};
        dur = str2num(words{4});
        mapping = words{5};
        nevfile = words{1};
        matfile = words{2};
        %From: 2013-01-08-1334
        curr_year = (datetime(1:4));
        curr_month = (datetime(6:7));
        curr_day = (datetime(9:10));
        curr_hr = (datetime(12:13));
        curr_min = (datetime(14:15));
        %To: 2000-03-01 15:45:17
        dateStr = [curr_year '-' curr_month '-' curr_day ' ' curr_hr ':' curr_min ':00'];
        curr_date = datenum(dateStr);
        try
            a = recordingtimes(mapping);
            a.durs = [a.durs; dur];
            a.dates = [a.dates; curr_date];
            recordingtimes(mapping) = a;
        end 
    end
    tline = fgetl(fh);
end
fclose(fh);

for idx = 1:length(conds)
    condition = conds{idx};
    a = recordingtimes(condition);
    totaldur = cumsum(a.durs)/3600;
    plot(a.dates, totaldur);
    title(condition);
    xlabel('Date')
    ylabel('Recording time (hrs)');
    datemin = min(a.dates);
    datemax = max(a.dates);
    if datemin ~= datemax
        set(gca, 'XTick', linspace(datemin, datemax, 6))
    end
    datetick('x','yyyy-mm-dd','keepticks')
    saveplot(gcf, ['./worksheets/01_13_2014/plots/' condition '_rectimes.eps'])
end