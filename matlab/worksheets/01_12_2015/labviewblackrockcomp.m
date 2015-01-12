nevfile = './testdata/20130117SpankyUtah001.nev';
matfile = './testdata/Spanky_2013-01-17-1325.mat';
binsize = 0.002;
offset = 0.0;
threshold = 5;
fn_out = './worksheets/01_12_2015/test_spline_pre_lv.eps';
processed = preprocess_spline_lv(nevfile, matfile, binsize, threshold, offset, fn_out);

%Also also a file that uses a BCI:
nevfile = './testdata/20130117SpankyUtah005.nev';
matfile = './testdata/Spanky_2013-01-17-1325.mat';
binsize = 0.002;
offset = 0.0;
threshold = 5;
fn_out = './worksheets/01_12_2015/test_spline_pre_BCI_lv.eps';
processed = preprocess_spline_lv(nevfile, matfile, binsize, threshold, offset, fn_out);
