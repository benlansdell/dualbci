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