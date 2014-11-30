fn_out = './worksheets/11_30_2014/plots/timescale';
nevfile = './testdata/20130117SpankyUtah001.nev';
const = 'on';
nK_sp = 100; 
nK_pos = 5;
threshold = 5;
offset = 0;
binsize = 0.002;
dt_sp = binsize;
dt_pos = 0.2;
alpha = 0.05;

processed = preprocess_spline(nevfile, binsize, threshold, offset);
data = filters_sp_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
model = MLE_glmfit(data, const);
gof_timerescale_ISI(model, data, processed, alpha, fn_out);