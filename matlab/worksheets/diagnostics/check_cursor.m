fn = './worksheets/diagnostics/plots/Spanky_2013-01-17-1325_trial';
ntrials = 1;

trials = import_trials('Spanky_2013-01-17-1325.mat')

i = 1;
j = 0;

while j < ntrials
	t = trials(i);
	if t.valid & t.success
		plot_trial(t, [fn num2str(i) '.gif']);
		j = j + 1;
	end
	i = i + 1;
end
