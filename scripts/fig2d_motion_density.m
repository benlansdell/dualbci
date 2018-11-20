%n_files = 10;
%n_pts = 1000;

%conn = database('spanky_db',databaseuser,databasepwd,'com.mysql.jdbc.Driver', ...
%	databaseurl);
%
%%Get a bunch of files used in analysis
%a = exec(conn, ['SELECT et1.`manualrecording`, et1.`1DBCrecording`, et1.`manualrecordingafter`, et1.`dualrecording`, rec1.`labview file` FROM '...
%'`experiment_tuning` et1 '...
%'INNER JOIN `recordings` rec1 '...
%'ON rec1.`nev file` = et1.`manualrecording` '...
%'WHERE et1.`manualrecording` IS NOT NULL AND et1.`1DBCrecording` IS NOT NULL AND et1.`manualrecordingafter` IS NOT NULL AND et1.`dualrecording` IS NOT NULL']);
%bci_data = fetch(a);
%bci_data = bci_data.Data;
%
%mcfiles = bci_data(1:n_files,1);
%bcfiles = bci_data(1:n_files,2);
%dcfiles = bci_data(1:n_files,4);
%labviewfiles = bci_data(1:n_files,5);
%
%%Extract wrist torque info from nev files
%pos_mctorque = zeros(n_files, n_pts, 2);
%pos_bctorque = zeros(n_files, n_pts, 2);
%pos_dctorque = zeros(n_files, n_pts, 2);
%
%vel_mctorque = zeros(n_files, n_pts, 2);
%vel_bctorque = zeros(n_files, n_pts, 2);
%vel_dctorque = zeros(n_files, n_pts, 2);
%
%%Load MC nevfile
%for idx = 1:n_files
%	nevfile = ['./blackrock/' mcfiles{idx}];
%	display(['Loading ' nevfile]);
%	matfile = ['./labview/' labviewfiles{idx}];
%	binsize = 0.002;
%	offset = 0.0;
%	threshold = 5;
%	processed = preprocess_spline_lv(nevfile, matfile, binsize, threshold, offset);
%	%Get n_pts samples from processed wrist torque
%	pos_mctorque(idx,:,:) = datasample(processed.torque, n_pts);
%	vel_mctorque(idx,:,:) = datasample(processed.dtorque, n_pts);
%end
%
%%Load BC nevfile
%for idx = 1:n_files
%	nevfile = ['./blackrock/' bcfiles{idx}];
%	display(['Loading ' nevfile]);
%	matfile = ['./labview/' labviewfiles{idx}];
%	binsize = 0.002;
%	offset = 0.0;
%	threshold = 5;
%	processed = preprocess_spline_lv(nevfile, matfile, binsize, threshold, offset);
%	%Get n_pts samples from processed wrist torque
%	pos_bctorque(idx,:,:) = datasample(processed.torque, n_pts);
%	vel_bctorque(idx,:,:) = datasample(processed.dtorque, n_pts);
%end
%
%%Load DC nevfile
%for idx = 1:n_files
%	nevfile = ['./blackrock/' dcfiles{idx}];
%	display(['Loading ' nevfile]);
%	matfile = ['./labview/' labviewfiles{idx}];
%	binsize = 0.002;
%	offset = 0.0;
%	threshold = 5;
%	processed = preprocess_spline_lv(nevfile, matfile, binsize, threshold, offset);
%	%Get n_pts samples from processed wrist torque
%	pos_dctorque(idx,:,:) = datasample(processed.torque, n_pts);
%	vel_dctorque(idx,:,:) = datasample(processed.dtorque, n_pts);
%end
%
%save('scripts/fig2d.mat')
load('scripts/fig2d.mat')

%Density plots
smoothness = 3;
nbins = 200;
vellim = 1;
poslim = 0.5;
cutoff = 0.01;

%MC density
pos_mctorque = reshape(pos_mctorque, [], 2);
vel_mctorque = reshape(vel_mctorque, [], 2);
subplot(3,2,1)
smoothhist2D(pos_mctorque, smoothness, [nbins, nbins], cutoff);
title('Manual control position')
xlim([-0.5, 0.5])
ylim([-0.5, 0.5])
subplot(3,2,2)
smoothhist2D(vel_mctorque, smoothness, [nbins, nbins], cutoff);
title('Manual control velocity')
xlim([-vellim, vellim])
ylim([-vellim, vellim])

%BC density plots
pos_bctorque = reshape(pos_bctorque, [], 2);
vel_bctorque = reshape(vel_bctorque, [], 2);
subplot(3,2,3)
smoothhist2D(pos_bctorque, smoothness, [nbins, nbins], cutoff);
title('Brain control position')
xlim([-0.5, 0.5])
ylim([-0.5, 0.5])
subplot(3,2,4)
smoothhist2D(vel_bctorque, smoothness, [nbins, nbins], cutoff);
title('Brain control velocity')
xlim([-vellim, vellim])
ylim([-vellim, vellim])

%DC density plots
pos_dctorque = reshape(pos_dctorque, [], 2);
vel_dctorque = reshape(vel_dctorque, [], 2);
subplot(3,2,5)
smoothhist2D(pos_dctorque, smoothness, [nbins, nbins], cutoff);
title('Dual control position')
xlim([-0.5, 0.5])
ylim([-0.5, 0.5])
subplot(3,2,6)
smoothhist2D(vel_dctorque, smoothness, [nbins, nbins], cutoff);
title('Dual control velocity')
xlim([-vellim, vellim])
ylim([-vellim, vellim])

saveplot(gcf, './figures/fig2d_motion_densities.eps', 'eps', [3 6]);