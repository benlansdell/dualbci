%Script to check that firing rate function is working...

trials = import_trials('Spanky_2013-01-17-1325.mat');
t = trials(117);
t = import_spikes(t);
fn = './worksheets/diagnostics/plots/check_electrodes_trial_117.gif';
plot_electrodes(t, fn);
