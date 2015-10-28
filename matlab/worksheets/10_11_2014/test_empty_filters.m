%%%%%%%%%%%%%%%%%%%%%%
%test_empty_filters.m%
%%%%%%%%%%%%%%%%%%%%%%

%Code to generate data in which the spike history filter is significant, but in which the stimulus filters are not important at all 
%in determining the spike rate

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

%Generate filters and data from GLM with those filters

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
nK_sp = 10;
nK_pos = 5;

t_sp = linspace(0,1,nK_sp);
%k_sp = 1.5*exp(-4*t_sp)-2*exp(-6*t_sp);
k_sp = -2*exp(-6*t_sp);
k_sp = fliplr(k_sp);

t_pos = linspace(0,1,nK_pos);

%Filter 1
k_RU1 = -0*25*exp(-6*t_pos);%.*sin(t_pos*2*pi);
k_FE1 = 0*25*exp(-5*t_pos);%.*sin(t_pos*2*pi);

%Filter 2
k_RU2 = 0*exp(-6*t_pos);%.*sin(t_pos*2*pi);
k_FE2 = 0*exp(-5*t_pos);%.*sin(t_pos*2*pi);

%Filter 3
k_RU3 = 0*exp(-6*t_pos);%.*sin(t_pos*2*pi);
k_FE3 = 0*exp(-5*t_pos);%.*sin(t_pos*2*pi);

%Simulate GLM
processed1 = generate_glm_data(freqlow, freqhigh, k_const, k_sp, k_RU1, k_FE1, N, binsize, seed);
processed2 = generate_glm_data(freqlow, freqhigh, k_const, k_sp, k_RU2, k_FE2, N, binsize, seed);
processed3 = generate_glm_data(freqlow, freqhigh, k_const, k_sp, k_RU3, k_FE3, N, binsize, seed);

%Combine into one structure
processed = processed2;
processed.binnedspikes = [processed1.binnedspikes, processed2.binnedspikes, processed3.binnedspikes];
processed.rates = [processed1.rates, processed2.rates, processed3.rates]; 
processed.unitnames = {'unit 1', 'unit 2', 'unit 3'};

%Reestimate filters
fn_out = './worksheets/10_11_2014/plots/testempty_fakedata';
data = filters_sp_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
k = data.k;
nK = size(k,1);
const = 'on';
order = 1;
model = MLE_glmfit(data, const);
plot_filters(model, data, processed, fn_out);

%Thus we see that in a case with no stimulus filters, trying to estimate the stim filters 
%gives something that looks a lot like what we see with the actual data... random oscillations
%This is for stimulus that is fairly slowly varying, though, what happens when we try the fitting 
%in the case when the stimulus is basically white noise?

%%%%%%%%%%%%%
%Second part%
%%%%%%%%%%%%%

%Redo fitting of 'empty' stim filters for the case when the input data is much noisier...

freqlow = 20;
%freqhigh = 60;
freqhigh = 2000;
N = 100000;
binsize = 0.002;
dur = N*binsize;
dt_sp = binsize;
dt_pos = 0.05;
seed = 1000000;
fn_out = './testdata/glm_gen_testdecode.mat'

%Filters
k_const = [0];
nK_sp = 10;
nK_pos = 5;

t_sp = linspace(0,1,nK_sp);
%k_sp = 1.5*exp(-4*t_sp)-2*exp(-6*t_sp);
k_sp = -2*exp(-6*t_sp);
k_sp = fliplr(k_sp);

t_pos = linspace(0,1,nK_pos);

%Filter 1
k_RU1 = -0*25*exp(-6*t_pos);%.*sin(t_pos*2*pi);
k_FE1 = 0*25*exp(-5*t_pos);%.*sin(t_pos*2*pi);

%Filter 2
k_RU2 = 0*exp(-6*t_pos);%.*sin(t_pos*2*pi);
k_FE2 = 0*exp(-5*t_pos);%.*sin(t_pos*2*pi);

%Filter 3
k_RU3 = 0*exp(-6*t_pos);%.*sin(t_pos*2*pi);
k_FE3 = 0*exp(-5*t_pos);%.*sin(t_pos*2*pi);

%Simulate GLM
processed1 = generate_glm_data(freqlow, freqhigh, k_const, k_sp, k_RU1, k_FE1, N, binsize, seed);
processed2 = generate_glm_data(freqlow, freqhigh, k_const, k_sp, k_RU2, k_FE2, N, binsize, seed);
processed3 = generate_glm_data(freqlow, freqhigh, k_const, k_sp, k_RU3, k_FE3, N, binsize, seed);

%Combine into one structure
processed = processed2;
processed.binnedspikes = [processed1.binnedspikes, processed2.binnedspikes, processed3.binnedspikes];
processed.rates = [processed1.rates, processed2.rates, processed3.rates]; 
processed.unitnames = {'unit 1', 'unit 2', 'unit 3'};

%Reestimate filters
fn_out = './worksheets/10_11_2014/plots/testempty_fakedata_highfreq';
data = filters_sp_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
k = data.k;
nK = size(k,1);
const = 'on';
order = 1;
model = MLE_glmfit(data, const);
plot_filters(model, data, processed, fn_out);

%Doesn't really look any different... of course not... if it's ignoring completely the stimulus then changing
%the stimulus properties won't have any effect!

%Try something slightly more nuanced

%%%%%%%%%%%%
%Third part%
%%%%%%%%%%%%

%Have the filters depend only a small amount on the stimulus, and vary the size of the noise...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%First make filters with small dependence and no high frequency noise%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
nK_sp = 10;
nK_pos = 5;

t_sp = linspace(0,1,nK_sp);
%k_sp = 1.5*exp(-4*t_sp)-2*exp(-6*t_sp);
k_sp = -2*exp(-6*t_sp);
k_sp = fliplr(k_sp);

t_pos = linspace(0,1,nK_pos);

%Filter 1
k_RU1 = 25*exp(-6*t_pos);%.*sin(t_pos*2*pi);
k_FE1 = -25*exp(-5*t_pos);%.*sin(t_pos*2*pi);

%Filter 2
k_RU2 = 20*exp(-6*t_pos);%.*sin(t_pos*2*pi);
k_FE2 = -5*exp(-5*t_pos);%.*sin(t_pos*2*pi);

%Filter 3
k_RU3 = 1*exp(-6*t_pos);%.*sin(t_pos*2*pi);
k_FE3 = -3*exp(-5*t_pos);%.*sin(t_pos*2*pi);

%Simulate GLM
processed1 = generate_glm_data(freqlow, freqhigh, k_const, k_sp, k_RU1, k_FE1, N, binsize, seed);
processed2 = generate_glm_data(freqlow, freqhigh, k_const, k_sp, k_RU2, k_FE2, N, binsize, seed);
processed3 = generate_glm_data(freqlow, freqhigh, k_const, k_sp, k_RU3, k_FE3, N, binsize, seed);

%Combine into one structure
processed = processed2;
processed.binnedspikes = [processed1.binnedspikes, processed2.binnedspikes, processed3.binnedspikes];
processed.rates = [processed1.rates, processed2.rates, processed3.rates]; 
processed.unitnames = {'unit 1', 'unit 2', 'unit 3'};

%Reestimate filters
fn_out = './worksheets/10_11_2014/plots/testsmall_fakedata';
data = filters_sp_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
k = data.k;
nK = size(k,1);
const = 'on';
order = 1;
model = MLE_glmfit(data, const);
plot_filters(model, data, processed, fn_out);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Then make filters with small dependence and with high frequency noise%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

freqlow = 20;
freqhigh = 2000;
N = 100000;
binsize = 0.002;
dur = N*binsize;
dt_sp = binsize;
dt_pos = 0.05;
seed = 1000000;
fn_out = './testdata/glm_gen_testdecode.mat'

%Filters
k_const = [0];
nK_sp = 10;
nK_pos = 5;

t_sp = linspace(0,1,nK_sp);
%k_sp = 1.5*exp(-4*t_sp)-2*exp(-6*t_sp);
k_sp = -2*exp(-6*t_sp);
k_sp = fliplr(k_sp);

t_pos = linspace(0,1,nK_pos);

%Filter 1
k_RU1 = 25*exp(-6*t_pos);%.*sin(t_pos*2*pi);
k_FE1 = -25*exp(-5*t_pos);%.*sin(t_pos*2*pi);

%Filter 2
k_RU2 = 20*exp(-6*t_pos);%.*sin(t_pos*2*pi);
k_FE2 = -5*exp(-5*t_pos);%.*sin(t_pos*2*pi);

%Filter 3
k_RU3 = 1*exp(-6*t_pos);%.*sin(t_pos*2*pi);
k_FE3 = -3*exp(-5*t_pos);%.*sin(t_pos*2*pi);

%Simulate GLM
processed1 = generate_glm_data(freqlow, freqhigh, k_const, k_sp, k_RU1, k_FE1, N, binsize, seed);
processed2 = generate_glm_data(freqlow, freqhigh, k_const, k_sp, k_RU2, k_FE2, N, binsize, seed);
processed3 = generate_glm_data(freqlow, freqhigh, k_const, k_sp, k_RU3, k_FE3, N, binsize, seed);

%Combine into one structure
processed = processed2;
processed.binnedspikes = [processed1.binnedspikes, processed2.binnedspikes, processed3.binnedspikes];
processed.rates = [processed1.rates, processed2.rates, processed3.rates]; 
processed.unitnames = {'unit 1', 'unit 2', 'unit 3'};

%Reestimate filters
fn_out = './worksheets/10_11_2014/plots/testsmall_fakedata_highfreq';
data = filters_sp_pos(processed, nK_sp, nK_pos, dt_sp, dt_pos);
k = data.k;
nK = size(k,1);
const = 'on';
order = 1;
model = MLE_glmfit(data, const);
plot_filters(model, data, processed, fn_out);

%Don't notice any effect of different stimulus statistics....Need to do more searching to nail down
%what causes this effect I suppose

save('./worksheets/10_11_2014/empty_filters_data.mat');