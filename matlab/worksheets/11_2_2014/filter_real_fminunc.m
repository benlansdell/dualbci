%%%%%%%%%%%%%%%
%filter_real.m%
%%%%%%%%%%%%%%%

%See how filters look on actual dataset
%Plot spatial filters as a function of resolution (for IRLS and fminunc)
%Plot likelihood (and AIC and BIC) as a function of resolution

fn_out = './worksheets/11_2_2014/data_real.mat';
N = 200000;
binsize = 0.002;
dur = N*binsize;
dt_sp = binsize;
nK_sp = 100;
filterlen = 500;
seed = 1000000;
const = 'on';

%Fit filters of different lengths/resolution to this data
nK_poss = [2 5 10 20 50 100];
dt_poss = binsize*filterlen./nK_poss;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nevfile = './testdata/20130117SpankyUtah001.nev';
threshold = 5; offset = 0;
pre = preprocess(nevfile, binsize, threshold, offset);
%Truncate to only one unit
idx = 9;
pre.binnedspikes = pre.binnedspikes(:,idx);
pre.rates = pre.rates(:,idx);              
pre.unitnames = pre.unitnames(idx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Can start from here by loading fn_out, as just saved above
fn_out = './worksheets/11_2_2014/data.mat';
save(fn_out);

fn_out = './worksheets/11_2_2014/data.mat';
load(fn_out);

%Compute auto-correlation for the two torque channels (within 3s)
maxlag = 4/dt_sp;

devs_IRLS = zeros(size(nK_poss));
devs_SD = zeros(size(nK_poss));
filt_IRLS = {};
filt_SD = {};
g = figure;clf
h = figure;clf
%Fit different filters for different timebin sizes
for j = 1:length(nK_poss);
	%set dt_pos accordingly
	dt_pos = dt_poss(j);
	nK_pos = nK_poss(j);
	%for each resolution generate the data structure
	data = filters_sp_pos(pre, nK_sp, nK_pos, dt_sp, dt_pos);
	%Takes a long time... will run later
	%model_SD = MLE_SD(data, const);
	model_IRLS = MLE_glmfit(data, const);
	%Compute the deviance
	%devs_SD(j) = deviance(model_SD, data);
	devs_IRLS(j) = deviance(model_IRLS, data);
	%Plot the filters estimated for each
	fn_out2 = ['./worksheets/11_2_2014/plots/filters_dtpos_' num2str(dt_pos)];
	%plot_filters(model_SD, data, pre, [fn_out2 '_fminunc'])
	plot_filters(model_IRLS, data, pre, [fn_out2 '_IRLS'])
	filt_IRLS{j} = model_IRLS.b_hat;
	%filt_SD{j} = model_SD.b_hat;
	for k = 1:size(data.k,1)
		%Plot the estimated filters
		name = data.k{k,1};
		filt = model_IRLS.b_hat(data.k{k,2}+1);
		filt_fmin = model_SD.b_hat(data.k{k,2}+1);
		if k == 1
			dt = dt_sp;
		else
			%rescale filter because of different time scales
			dt = dt_pos;
			filt = filt*dt_sp/dt_pos;
			filt_fmin = filt_fmin*dt_sp/dt_pos;
		end
		figure(g);
		subplot(3,1,k)
		hold on
		plot((0:length(filt)-1)*dt, filt, 'Color', [1 0.5 0.5])  
		plot((0:length(filt)-1)*dt, filt, 'r.')  
		plot((0:length(filt_fmin)-1)*dt, filt_fmin, 'Color', [0.5 0.5 1])  
		plot((0:length(filt_fmin)-1)*dt, filt_fmin, 'b.')  
		title(name);
		figure(h);
		subplot(3,4,(4*(k-1)+(1:3)))
		hold on
		NFFT = 2^nextpow2(length(filt)); % Next power of 2 from length of y
		Fs = 1/dt;
		Y = fft(filt,NFFT)/length(filt);
		f = Fs/2*linspace(0,1,NFFT/2+1);
		% Plot single-sided amplitude spectrum.
		plot(f,2*abs(Y(1:NFFT/2+1)), 'Color', [1 0.5 0.5])
		plot(f,2*abs(Y(1:NFFT/2+1)), 'r.')
		NFFT = 2^nextpow2(length(filt_fmin)); % Next power of 2 from length of y
		Fs = 1/dt;
		Y = fft(filt_fmin,NFFT)/length(filt_fmin);
		f = Fs/2*linspace(0,1,NFFT/2+1);
		% Plot single-sided amplitude spectrum.
		plot(f,2*abs(Y(1:NFFT/2+1)), 'Color', [0.5 0.5 1])
		plot(f,2*abs(Y(1:NFFT/2+1)), 'b.')
		title(name);
	end
end
autocorrRU = xcorr(pre.torque(:,1), maxlag);
autocorrFE = xcorr(pre.torque(:,2), maxlag);
autocorrdRU = xcorr(pre.dtorque(:,1), maxlag);
autocorrdFE = xcorr(pre.dtorque(:,2), maxlag);
%Plot the actual filters
figure(g);
subplot(3,1,1)
title({['Fit filters.'];['spike history filter']})
hold on
subplot(3,1,2)
hold on
subplot(3,1,3)
hold on
xlabel('time (s)')
fn_out3 = ['./worksheets/11_2_2014/plots/filters_fminunc.eps'];
saveplot(gcf, fn_out3, 'eps', [9 6]);
%Plot the actual filters in fourier domain
figure(h);
subplot(3,4,[1 2 3])
title({['Fourier transform of fit filters.'];['spike history filter']})
hold on
ylabel('|\beta (k)|')
subplot(3,4,[5 6 7])
hold on
ylabel('|\beta (k)|')
subplot(3,4,[9 10 11])
hold on
ylabel('|\beta (k)|')
xlabel('Freq (Hz)')
subplot(3,4,8)
tt = ((1:length(autocorrRU))-length(autocorrRU)/2)*dt_sp;
plot(tt,autocorrRU)
title('Auto-correlation RU')
subplot(3,4,12)
plot(tt,autocorrFE)
title('Auto-correlation FE')
xlabel('time(s)')
fn_out3 = ['./worksheets/11_2_2014/plots/filters_fminunc_fourier.eps'];
saveplot(gcf, fn_out3, 'eps', [9 6]);

%Save the outcome of all the above
save(fn_out);
%fn_out = './worksheets/10_27_2014/data.mat';
%load(fn_out);

%Compute AIC and BIC
nK = [nK_poss]+nK_sp + 1;
AIC_IRLS = devs_IRLS + 2*nK;
AICc_IRLS = AIC_IRLS + 2*nK.*(nK+1)./(N-1-nK);
BIC_IRLS = devs_IRLS + ([nK_poss]+nK_sp+1)*log(N);

%AIC_SD = devs_SD + 2*nK;
%AICc_SD = AIC_SD + 2*nK.*(nK+1)./(N-1-nK);
%BIC_SD = devs_SD + ([nK_poss]+nK_sp+1)*log(N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make a plot of the deviance as a function of params
clf
plot(dt_poss, devs_IRLS, dt_poss, AIC_IRLS, dt_poss, AICc_IRLS, dt_poss, BIC_IRLS);
legend('Deviance', 'AIC', 'AICc', 'BIC')
%hold on
%plot(dt_poss, devs_SD, '--', dt_poss, AIC_SD, '--', dt_poss, AICc_SD, '--', dt_poss, BIC_SD, '--');
xlabel('resolution (s)')
ylabel('likelihood')
%title('Dotted line is fminunc, solid is IRLS')
saveplot(gcf, './worksheets/11_2_2014/plots/filters_real_dev.eps')

%Plot auto correlation in separate plot
figure 
subplot(2,2,1)
tt = ((1:length(autocorrRU))-length(autocorrRU)/2)*dt_sp;
plot(tt,autocorrRU)
xlabel('time(s)')
title('Auto-correlation RU')
subplot(2,2,2)
tt = ((1:length(autocorrFE))-length(autocorrFE)/2)*dt_sp;
plot(tt,autocorrFE)
xlabel('time(s)')
title('Auto-correlation FE')
subplot(2,2,3)
plot(tt,autocorrdFE)
title('Auto-correlation vel RU')
xlabel('time(s)')
subplot(2,2,4)
tt = ((1:length(autocorrdRU))-length(autocorrdRU)/2)*dt_sp;
plot(tt,autocorrdRU)
xlabel('time(s)')
title('Auto-correlation vel FE')
fn_out4 = ['./worksheets/11_2_2014/plots/filters_fminunc_autocorr.eps'];
saveplot(gcf, fn_out4, 'eps', [6 4]);

%Compute norms of fitted filters to each other and to actual
inf_filts = [];
rel_filts = [];
max_filts_SD = [];
max_filts_IRLS = [];
FI_SD = {};
FI_IRLS = {};
kappa_SD = [];
kappa_IRLS = [];
for idx = 1:length(nK_poss)
	idx
	diff_filts = filt_IRLS{idx}-filt_SD{idx};
	inf_filts(idx) = norm(diff_filts, Inf);
	max_filts_SD(idx) = norm(filt_SD{idx}, Inf);
	max_filts_IRLS(idx) = norm(filt_IRLS{idx}, Inf);
	rel_filts(idx) = inf_filts(idx)/max_filts_IRLS(idx);
	%Fisher information
	%set dt_pos accordingly
	dt_pos = dt_poss(idx);
	nK_pos = nK_poss(idx);
	%recompute data matrix structure for each resolution generate the data structure
	data = filters_sp_pos(pre, nK_sp, nK_pos, dt_sp, dt_pos);
	N = size(data.X,2);
	X = [ones(N,1),squeeze(data.X)];
	beta_IRLS = filt_IRLS{idx}';
	beta_SD = filt_SD{idx}';
	mu_IRLS = spdiags(exp(X*beta_IRLS), [0], N,N);
	mu_SD = spdiags(exp(X*beta_SD), [0], N,N);
	FI_SD{idx} = X'*mu_IRLS*X;
	FI_IRLS{idx} = X'*mu_SD*X;
	kappa_IRLS(idx) = cond(FI_IRLS{idx});
	kappa_SD(idx) = cond(FI_SD{idx});
end

%Compute condition number of Fisher information matrices
save(fn_out);