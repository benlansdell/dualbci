%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%filter_strengths_fminunc.m%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Code to see how fitting performs on GLMs generated from known filters of different
%strengths, with actual stimulus data used and generated spike trains from these GLMs
%We compare the quality of the fit to the filter strength and the #of parameters of the fit model

clear

fn_out = './worksheets/11_12_2014/data.mat';
N = 200000;
binsize = 0.002;
dur = N*binsize;
dt_sp = binsize;
dt_pos = binsize;
seed = 1000000;
const = 'on';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Make up 10 filters. They'll have a common spike history filter, and varying stimulus strengths
%Change the constant term such that they each have the same average firing rate

%Filters
nK_sp = 100;
nK_pos = 100;
%Factors to scale stim filters by
k_const = zeros(nAlpha,1);
target_nsp = 10000;
k_const_guess = -5;

%Spike history filter
t_sp = linspace(0,1,nK_sp);
k_sp = -0.1*exp(-50*t_sp);
k_sp = fliplr(k_sp);

t_pos = linspace(0,1,nK_pos);
k_RU = -1.0*exp(-6*t_pos);
k_FE = 0.5*exp(-5*t_pos);


nevfile = './testdata/20130117SpankyUtah001.nev';
threshold = 5; offset = 0;
pre = preprocess(nevfile, binsize, threshold, offset);
%Truncate to only one unit
idx = 9;
pre.binnedspikes = pre.binnedspikes(:,idx);
pre.rates = pre.rates(:,idx);              
pre.unitnames = pre.unitnames(idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Generate spike train data with these filters and actual stimulus data as input
%Tweak the constant term so that the average firing rate is about what is seen in the actual data sets
%Fit filters of different lengths/resolution to this data

nK_poss = [2 5 10 20 50 100];
dt_poss = binsize*100./nK_poss;

processed = {};
K = {};
%Determine how to set constant rate from train with the same filters
for j = 1:length(nK_poss)
	nK_pos = nK_poss(j);
	dt_pos = dt_poss(j);
	%Base stim filter that will be scaled below
	t_pos = linspace(0,1,nK_pos);
	k_RU_sc = -1.0*exp(-6*t_pos);
	k_FE_sc = 0.5*exp(-5*t_pos);
	%Rescale since binsize is changing
	k_RU_sc = k_RU_sc*dt_pos/dt_sp;
	k_FE_sc = k_FE_sc*dt_pos/dt_sp;
	K{1,j} = k_sp;
	K{2,j} = k_RU_sc;
	K{3,j} = k_FE_sc;
	%Simulate GLM with filters wanted but constant term left at a guess value
	p = generate_glm_data_torque(pre, k_const_guess, k_sp, k_RU_sc, k_FE_sc, dt_sp, dt_pos, N, binsize);
	%Count number of spikes
	nspikes1 = sum(p.binnedspikes);
	%Estimate constant to produce wanted average firing rates
	k_const(j) = k_const_guess+log(target_nsp/nspikes1);
	%Redo with estimated constant values
	processed{j} = generate_glm_data_torque(pre, k_const(idx), k_sp, k_RU_sc, k_FE_sc, dt_sp, dt_pos, N, binsize);
	%Check that total number of spikes is around the same...
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Can start from here by loading fn_out, as just saved above
fn_out = './worksheets/11_12_2014/data.mat';
save(fn_out);

fn_out = './worksheets/11_12_2014/data.mat';
load(fn_out);

%Compute auto-correlation for the two torque channels (within 3s)
maxlag = 4/dt_sp;
autocorrRU = xcorr(pre.torque(:,1), maxlag);
autocorrFE = xcorr(pre.torque(:,2), maxlag);

devs_IRLS = zeros(length(nK_poss));
filt_IRLS = {};

g = figure;
h = figure;
for i = 1:length(nK_poss)
	i
	figure(g);
	clf;
	dt_pos_act = dt_poss(i);
	for j = 1:length(nK_poss);
		%set dt_pos accordingly
		dt_pos = dt_poss(j);
		nK_pos = nK_poss(j);
		%for each resolution generate the data structure
		data = filters_sp_pos(processed{i}, nK_sp, nK_pos, dt_sp, dt_pos);
		%Takes a long time... will run later
		model_IRLS = MLE_glmfit(data, const);
		%Compute the deviance
		devs_IRLS(i,j) = deviance(model_IRLS, data);
		%Plot the filters estimated for each
		fn_out2 = ['./worksheets/11_12_2014/plots/filterstrengths_dtposact_' num2str(dt_pos_act) '_dtpos_' num2str(dt_pos)];
		filt_IRLS{i,j} = model_IRLS.b_hat;
		figure(g);
	end
end
%Save the outcome of all the above
save(fn_out);
fn_out = './worksheets/11_12_2014/data.mat';
load(fn_out);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make a plot of the deviance as a function of params

nK = repmat([nK_poss]+nK_sp + 1, nAlpha, 1);
AIC_IRLS = devs_IRLS + 2*nK;
AICc_IRLS = AIC_IRLS + 2*nK.*(nK+1)./(N-1-nK);
BIC_IRLS = devs_IRLS + nK.*log(N);

clf
imagesc(devs_IRLS)
title('Deviance')
xlabel('resolution (s)')
ylabel('alpha')
%set(gca,'XTick',1:size(S,1));
set(gca,'XTickLabel',dt_poss);
%set(gca,'YTick',1:size(S,2));
set(gca,'YTickLabel',alphas);
colorbar
saveplot(gcf, './worksheets/11_12_2014/plots/filterstrengths_dev.eps')

clf
imagesc(AIC_IRLS)
title('AIC')
xlabel('resolution (s)')
ylabel('alpha')
%set(gca,'XTick',1:size(S,1));
set(gca,'XTickLabel',dt_poss);
%set(gca,'YTick',1:size(S,2));
set(gca,'YTickLabel',alphas);
colorbar
saveplot(gcf, './worksheets/11_12_2014/plots/filterstrengths_AIC.eps')

clf
imagesc(BIC_IRLS)
title('BIC')
xlabel('resolution (s)')
ylabel('alpha')
%set(gca,'XTick',1:size(S,1));
set(gca,'XTickLabel',dt_poss);
%set(gca,'YTick',1:size(S,2));
set(gca,'YTickLabel',alphas);
colorbar
saveplot(gcf, './worksheets/11_12_2014/plots/filterstrengths_BIC.eps')

%Compute dot product of fitted filters and actual filters and plot result
%K{1,idx,j} = k_sp;
%K{2,idx,j} = k_RU_sc;
%K{3,idx,j} = k_FE_sc;
%k_const(idx,j) = k_const_guess+log(target_nsp/nspikes1);

%Must interpolate so that all filters are the same length

%Fitted 
%filt_IRLS
dp = [];
for i = 1:nAlpha
	dt_pos_act = dt_poss(i);
	act_filt = [k_const(i), K{1,i},K{2,i},K{3,i}];
	%Fit spline to cursor filters and resample at binsize to get high res filter

	for j = 1:length(nK_poss)
			dt_pos = dt_poss(j);
			fit_filt = filt_IRLS{i,j};
			asdf
			dp(i,j) = act_filt*filt_IRLS{i,j}'/norm(act_filt)/norm(filt_IRLS{i,j});
	end
end
clf
imagesc(dp)
title('k.k_{IRLS}/|k||k_{IRLS}|')
xlabel('resolution fitted (s)')
ylabel('resolution actual (s)')
%set(gca,'XTick',1:size(S,1));
set(gca,'XTickLabel',dt_poss);
%set(gca,'YTick',1:size(S,2));
set(gca,'YTickLabel',alphas);
colorbar
saveplot(gcf, './worksheets/11_12_2014/plots/filterstrengths_DP.eps')
