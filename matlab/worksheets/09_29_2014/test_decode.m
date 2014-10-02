%%%%%%%%%%%%%%%
%test_decode.m%
%%%%%%%%%%%%%%%

%Code to test my glm_decode function works

sigma = 0.25;
%sigma = 0.001;
sigma = sigma/processed.binsize;
sz = sigma*3*3;
x = linspace(-sz/2, sz/2, sz);
gaussFilter = exp(-x.^2/(2*sigma^2));
gaussFilter = gaussFilter/sum(gaussFilter);

%%%%%%%%%%%%
%First part%
%%%%%%%%%%%%

%Generate filters and data from GLM with those filters to test decoder works...

freqlow = 20;
freqhigh = 60;
%freqhigh = 2000;
N = 100000;
binsize = 0.002;
dur = N*binsize;
dt_sp = binsize;
dt_pos = 0.05;
seed = 1000000;
fn_out = './testdata/glm_gen_testdecode.mat'

%Filters
k_const = [0];
nK_sp = 50;
nK_pos = 1;

t_sp = linspace(0,1,nK_sp);
%k_sp = 1.5*exp(-4*t_sp)-2*exp(-6*t_sp);
k_sp = -2*exp(-6*t_sp);
k_sp = fliplr(k_sp);

t_pos = 0; %linspace(0,1,nK_pos);

%Filter 1
k_RU1 = -4*25*exp(-6*t_pos);%.*sin(t_pos*2*pi);
k_FE1 = 6*25*exp(-5*t_pos);%.*sin(t_pos*2*pi);

%Filter 2
k_RU2 = 1*25*exp(-6*t_pos);%.*sin(t_pos*2*pi);
k_FE2 = 3*25*exp(-5*t_pos);%.*sin(t_pos*2*pi);

%Filter 3
k_RU3 = -1*25*exp(-6*t_pos);%.*sin(t_pos*2*pi);
k_FE3 = 3*25*exp(-5*t_pos);%.*sin(t_pos*2*pi);

%Simulate GLM
processed1 = generate_glm_data(freqlow, freqhigh, k_const, k_sp, k_RU1, k_FE1, N, binsize, seed);
processed2 = generate_glm_data(freqlow, freqhigh, k_const, k_sp, k_RU2, k_FE2, N, binsize, seed);
processed3 = generate_glm_data(freqlow, freqhigh, k_const, k_sp, k_RU3, k_FE3, N, binsize, seed);

%Combine into one structure
processed = processed2;
processed.binnedspikes = [processed1.binnedspikes, processed2.binnedspikes, processed3.binnedspikes];
processed.rates = [processed1.rates, processed2.rates, processed3.rates]; 
processed.unitnames = {'unit 1', 'unit 2', 'unit 3'};

%Save data
save(fn_out, 'processed');

%Reestimate filters
fn_out = './worksheets/09_29_2014/plots/testdecode_fakedata';
data = filters_sp_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
k = data.k;
nK = size(k,1);
const = 'on';
order = 1;
model = MLE_glmfit(data, const);
[F, Q, mu] = fit_AR_LS(data.torque, order);

%load('./testdata/testglm.mat')
%Decode based on filters
decoded_torque = glm_decode(processed, data, model, F, Q, mu, fn_out);

%Seems to work well on fabricated data (that probably has a very high firing rate...)

%%%%%%%%%%%%%
%Second part%
%%%%%%%%%%%%%

%Read in recording where the cursor position was infact mapped from firing rates and repeat

nevfile = './testdata/20130117SpankyUtah005.nev';
fn_out = './worksheets/09_29_2014/testdecode_braincontrol.eps';
fn_out_filters = './worksheets/09_29_2014/testdecode_filters';
%Load the mat file with some info on the BCI
load('./testdata/Spanky_2013-01-17-1325.mat')
binsize = 0.002;
offset = 0;
threshold = 5;
const = 'on';
nK_sp = 1;
nK_pos = 5;
dt_sp = 0.002;
dt_pos = 0.05;
order = 1;

processed = preprocess(nevfile, binsize, threshold, offset);
data = filters_sp_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
%Save a short version of the data for later...
%Load the short version of data
smthratesU = conv(processed.binnedspikes(:,12), gaussFilter, 'same');
smthRU = conv(processed.torque(:,1), gaussFilter, 'same');
smthFE = conv(processed.torque(:,2), gaussFilter, 'same');
t_i = 60;
%t_f = 33;
t_f = 200; 
ii = 1:size(processed.torque,1);
tt = ii*processed.binsize;
clf
plot(tt, smthratesU, tt, smthRU, tt, smthFE);
xlim([t_i t_f])

model = MLE_glmfit(data, const);
%Plot filters fitted... we should hopefully observe something significant! Or else the fitting must be wrong
plot_filters(model, data, processed, fn_out_filters);
plot_predictions(model, data, processed, fn_out);

[F, Q, mu] = fit_AR_LS(data.torque, order);
decoded_torque = glm_decode(processed, data, model, F, Q, mu, fn_out);

%%%%%%%%%%%%
%Third part%
%%%%%%%%%%%%

%Run the decoder on a manual control dataset...

nevfile = './testdata/20130117SpankyUtah001.nev';
fn_out = './worksheets/09_29_2014/testdecode_braincontrol.eps';
%fn_out_filters = './worksheets/09_29_2014/testdecode_filters';
%Load the mat file with some info on the BCI
%load('./testdata/Spanky_2013-01-17-1325.mat')
binsize = 0.002;
offset = -0.150;
threshold = 5;
const = 'on';
nK_sp = 1;
nK_pos = 1;
dt_sp = 0.002;
dt_pos = 0.05;
order = 1;

processed = preprocess(nevfile, binsize, threshold, offset);
data = filters_sp_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
%Save a short version of the data for later...
%Load the short version of data
%smthratesU = conv(processed.binnedspikes(:,12), gaussFilter, 'same');
%smthRU = conv(processed.torque(:,1), gaussFilter, 'same');
%smthFE = conv(processed.torque(:,2), gaussFilter, 'same');
%t_i = 60;
%t_f = 33;
%t_f = 200; 
%ii = 1:size(processed.torque,1);
%tt = ii*processed.binsize;
%clf
%plot(tt, smthratesU, tt, smthRU, tt, smthFE);
%xlim([t_i t_f])

model = MLE_glmfit(data, const);
%Plot filters fitted... we should hopefully observe something significant! Or else the fitting must be wrong
%plot_filters(model, data, processed, fn_out_filters);
%plot_predictions(model, data, processed, fn_out);
[F, Q, mu] = fit_AR_LS(data.torque, order);
figure;
decoded_torque = glm_decode(processed, data, model, F, Q, mu, fn_out);