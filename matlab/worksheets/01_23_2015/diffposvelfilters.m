nevfile = './testdata/20130117SpankyUtah001.nev';
binsize = 0.002;
dt_sp = binsize;
dt_pos = 0.1;
dt_vel = 0.1;
offset = 0.0;
threshold = 5;
processed = preprocess_spline(nevfile, binsize, threshold, offset);
nU = length(processed.unitnames);
N = size(processed.rates, 1);
const = 'on';
%Prepare data for fitting GLM
%Prepare:
nK_sp = 50; 
nK_pos = 10;

%Position filter
data = filters_sp_pos(processed, nK_sp, nK_pos);
%Fit each of the above GLMs
model1 = MLE_glmfit(data, const);
%Plot filters
fn_out1 = './worksheets/01_23_2015/plots/pos';
plot_filters(model1, data, processed, fn_out1);

%Velocity filter
data = filters_sp_vel(processed, nK_sp, nK_pos);
%Fit each of the above GLMs
model2 = MLE_glmfit(data, const);
%Plot filters
fn_out2 = './worksheets/01_23_2015/plots/vel_inst';
plot_filters(model2, data, processed, fn_out2);

%Velocity filter based on difference between adjacent columns of data matrix
%data = filters_sp_pos(processed, nK_sp, nK_pos);
dataveldisc = data;
dataveldisc.X = data.X(:,:,[1:59,61:69]);
dataveldisc.X(:,:,51:59)=diff(data.X(:,:,[51:60]),1,3);
dataveldisc.X(:,:,60:68)=diff(data.X(:,:,[61:70]),1,3);
dataveldisc.k{2,2} = 51:59;
dataveldisc.k{3,2} = 60:68;
clear data
%Fit each of the above GLMs
model3 = MLE_glmfit(dataveldisc, const);
%Plot filters
fn_out3 = './worksheets/01_23_2015/plots/vel_disc';
plot_filters(model3, dataveldisc, processed, fn_out3);

%Check if filter for unit 55.1 disc is equal once converted