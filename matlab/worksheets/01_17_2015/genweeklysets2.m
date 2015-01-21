csv = './nevmappings2.csv';
duration = 360; %at least six minutes
fn_out = './worksheets/01_17_2015/weeklynevs.mat';
conds = {'2D Manual Position',
    'Dual Control'};
weeklynevs = containers.Map;
%Check inside genweeklydataset to see how often the function is collecting .nev files...
%Here check that it's collecting them weekly
for idx = 1:length(conds)
	condition = conds{idx};
    weeklynevs(condition) = genweeklydataset(csv, condition, duration);
end
save(fn_out, 'weeklynevs');

%Copy them to americano...
%Labview source dir
labviewsrc = '/home/lansdell/projects/bci/matlab/labview/';
%Blackrock source dir
blackrocksrc = '/home/lansdell/projects/bci/matlab/blackrock/';
%Labview target dir
labviewtgt = '/home/lansdell/americano/projects/bci/matlab/labview/';
%Blackrock target dir
blackrocktgt = '/home/lansdell/americano/projects/bci/matlab/blackrock/';

for idx = 1:length(conds)
	condition = conds{idx};
	nevs = weeklynevs(condition);
	for j = 1:length(nevs)
		nevfile = [blackrocksrc nevs(j).nevfile]
		matfile = [labviewsrc nevs(j).matfile];
		if (length(nevs(j).nevfile)>1) & (length(nevs(j).matfile)>1)
			system(['cp ' nevfile ' ' blackrocktgt])
			system(['cp ' matfile ' ' labviewtgt])
			%copyfile(nevfile,blackrocktgt, 'f');
			%copyfile(matfile,labviewtgt,'f');
		end
	end
end

%Now they've been copied across...
%Use grangerweekly.m script to run granger analysis on the files